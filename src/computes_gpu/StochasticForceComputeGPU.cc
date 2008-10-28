/*
Highly Optimized Object-Oriented Molecular Dynamics (HOOMD) Open
Source Software License
Copyright (c) 2008 Ames Laboratory Iowa State University
All rights reserved.

Redistribution and use of HOOMD, in source and binary forms, with or
without modification, are permitted, provided that the following
conditions are met:

* Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names HOOMD's
contributors may be used to endorse or promote products derived from this
software without specific prior written permission.

Disclaimer

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND
CONTRIBUTORS ``AS IS''  AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 

IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS  BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.
*/

// $Id: StochasticForceComputeGPU.cc 1234 2008-09-11 16:29:13Z joaander $
// $URL: http://svn2.assembla.com/svn/hoomd/trunk/src/computes_gpu/StochasticForceComputeGPU.cc $

/*! \file StochasticForceComputeGPU.cc
	\brief Defines the StochasticForceComputeGPU class
*/

#ifdef WIN32
#pragma warning( push )
#pragma warning( disable : 4103 4244 )
#endif

#include "StochasticForceComputeGPU.h"
#include "cuda_runtime.h"

#include <stdexcept>
#include <stdlib.h>
#include <math.h>

#include <boost/python.hpp>
using namespace boost::python;

#include <boost/bind.hpp>

using namespace boost;
using namespace std;

/*! \param pdata ParticleData to compute forces on
	\param Temp Temperature of the bath of random particles
	\param deltaT Length of the computation timestep
	\param seed	Seed for initializing the RNG
*/
StochasticForceComputeGPU::StochasticForceComputeGPU(boost::shared_ptr<ParticleData> pdata, Scalar deltaT, Scalar Temp, unsigned int seed) 
	: StochasticForceCompute(pdata, deltaT, Temp, seed)
	{
	const ExecutionConfiguration& exec_conf = m_pdata->getExecConf();

//  I DO NOT KNOW IF THIS APPLIES TO THIS FORCE		
	// default block size is the highest performance in testing on different hardware
	// choose based on compute capability of the device
	cudaDeviceProp deviceProp;
	int dev;
	exec_conf.gpu[0]->call(bind(cudaGetDevice, &dev));	
	exec_conf.gpu[0]->call(bind(cudaGetDeviceProperties, &deviceProp, dev));
	if (deviceProp.major == 1 && deviceProp.minor < 2)
		m_block_size = 128;
	else if (deviceProp.major == 1 && deviceProp.minor < 4)
		m_block_size = 96;
	else
		{
		cout << "***Warning! Unknown compute " << deviceProp.major << "." << deviceProp.minor << " when tuning block size for StochasticForceComputeGPU" << endl;
		m_block_size = 96;
		}


	// allocate the gamma data on the GPU
	int nbytes = sizeof(float)*m_pdata->getNTypes();
	// allocate the coeff data on the CPU
	h_gammas = new float[m_pdata->getNTypes()];
	//All gamma coefficients initialized to 1.0
	for (unsigned int j = 0; j < m_pdata->getNTypes(); j++) h_gammas[j] = 1.0;  
	d_gammas.resize(exec_conf.gpu.size());
	
	for (unsigned int cur_gpu = 0; cur_gpu < exec_conf.gpu.size(); cur_gpu++)
		{
		exec_conf.gpu[cur_gpu]->setTag(__FILE__, __LINE__);
		exec_conf.gpu[cur_gpu]->call(bind(cudaMalloc, (void **)((void *)&d_gammas[cur_gpu]), nbytes));
		assert(d_gammas[cur_gpu]);
		exec_conf.gpu[cur_gpu]->call(bind(cudaMemcpy, (void **)(void *) d_gammas[cur_gpu], h_gammas, nbytes, cudaMemcpyHostToDevice));
		}

	}
	

StochasticForceComputeGPU::~StochasticForceComputeGPU()
	{
	// deallocate our memory
	const ExecutionConfiguration& exec_conf = m_pdata->getExecConf();
	for (unsigned int cur_gpu = 0; cur_gpu < exec_conf.gpu.size(); cur_gpu++)
		{
		assert(d_gammas[cur_gpu]);
		exec_conf.gpu[cur_gpu]->setTag(__FILE__, __LINE__);
		exec_conf.gpu[cur_gpu]->call(bind(cudaFree, (void *)d_gammas[cur_gpu]));
		}
	delete[] h_gammas;
	}	
		
/*! \param block_size Size of the block to run on the device
	Performance of the code may be dependant on the block size run
	on the GPU. \a block_size should be set to be a multiple of 32.
	\todo error check value
*/
void StochasticForceComputeGPU::setBlockSize(int block_size)
	{
	m_block_size = block_size;
	}

/*! \post The parameter \a gamma is set for \a typ, 
	\note \a gamma is a low level parameters used in the calculation. 
	
	\param typ Specifies the particle type
	\param gamma Parameter used to calcluate forces
*/
void StochasticForceComputeGPU::setParams(unsigned int typ, Scalar gamma)
	{
	assert(h_gammas);
	if (typ >= m_ntypes)
		{
		cerr << endl << "***Error! Trying to set Stochastic Force param Gamma for a non existant type! " << typ << endl << endl;
		throw runtime_error("StochasticForceComputeGpu::setParams argument error");
		}
	
	// set gamma coeffs 
	h_gammas[typ] = gamma;
	
	int nbytes = sizeof(float)*m_pdata->getNTypes();
	const ExecutionConfiguration& exec_conf = m_pdata->getExecConf();
	for (unsigned int cur_gpu = 0; cur_gpu < exec_conf.gpu.size(); cur_gpu++)
		exec_conf.gpu[cur_gpu]->call(bind(cudaMemcpy, d_gammas[cur_gpu], h_gammas, nbytes, cudaMemcpyHostToDevice));
	}

/*! \post The Temeperature of the Stochastic Bath \a T
	\note \a T is a low level parameter used in the calculation. 
	
	\param T Temperature of Stochastic Bath
*/

void StochasticForceComputeGPU::setT(Scalar T)
	{
	if (T <= 0)
		{
		cerr << endl << "***Error! Trying to set a Temperature <= 0 " << endl << endl;
		throw runtime_error("StochasticForceComputeGpu::setT argument error");
		}
	
	// set Temperature
	m_T = T;	
	cout << "Set T to " << m_T << endl;	
	}	

/*! \post The timestep of the Stochastic Bath \a T
	\note \a deltaT is a low level parameter used in the calculation. 
	
	\param deltaT timestep of Stochastic Bath
*/

void StochasticForceComputeGPU::setDeltaT(Scalar deltaT)
	{
	if (deltaT <= 0)
		{
		cerr << endl << "***Error! Trying to set a timestep <= 0 " << endl << endl;
		throw runtime_error("StochasticForceComputeGpu::setDeltaT argument error");
		}
	
	m_dt=deltaT;
	
		
	}	

		
/*! \post The stochastic forces are computed for the given timestep on the GPU. 
 	\param timestep Current time step of the simulation
*/
void StochasticForceComputeGPU::computeForces(unsigned int timestep)
	{
	// check the execution configuration
	const ExecutionConfiguration& exec_conf = m_pdata->getExecConf();
	
	// start the profile
	if (m_prof) m_prof->push(exec_conf, "Stochastic Baths");

	// access the particle data
	vector<gpu_pdata_arrays>& pdata = m_pdata->acquireReadOnlyGPU();
	
	for (unsigned int cur_gpu = 0; cur_gpu < exec_conf.gpu.size(); cur_gpu++)
		{
		exec_conf.gpu[cur_gpu]->setTag(__FILE__, __LINE__);
		exec_conf.gpu[cur_gpu]->callAsync(bind(gpu_stochasticforce, m_gpu_forces[cur_gpu].d_data, &pdata[cur_gpu], m_dt, m_T, d_gammas[cur_gpu], m_seed, timestep, m_pdata->getNTypes(), m_block_size));
		}
	for (unsigned int cur_gpu = 0; cur_gpu < exec_conf.gpu.size(); cur_gpu++)
		exec_conf.gpu[cur_gpu]->sync();
	
	m_pdata->release();
	
	// the force data is now only up to date on the gpu
	m_data_location = gpu;

//	int64_t mem_transfer = m_pdata->getN() * (4 + 16 + 16) + n_calc * (4 + 16);
//	int64_t flops = n_calc * (3+12+5+2+2+6+3+7);
//	if (m_prof) m_prof->pop(exec_conf, flops, mem_transfer);
	
	}

#ifdef WIN32
#pragma warning( pop )
#endif


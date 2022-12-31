// Copyright (c) 2009-2022 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include "TwoStepNVTMTK.h"

#ifndef __TWO_STEP_NVT_MTK_GPU_H__
#define __TWO_STEP_NVT_MTK_GPU_H__

/*! \file TwoStepNVTMTKGPU.h
    \brief Declares the TwoStepNVTMTKGPU class
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#include <pybind11/pybind11.h>

#include "hoomd/Autotuner.h"

namespace hoomd
    {
namespace md
    {
//! Integrates part of the system forward in two steps in the NVT ensemble on the GPU
/*! Implements Nose-Hoover NVT integration through the IntegrationMethodTwoStep interface, runs on
   the GPU

    In order to compute efficiently and limit the number of kernel launches integrateStepOne()
   performs a first pass reduction on the sum of m*v^2 and stores the partial reductions. A second
   kernel is then launched to reduce those to a final \a sum2K, which is a scalar but stored in a
   GPUArray for convenience.

    \ingroup updaters
*/
class PYBIND11_EXPORT TwoStepNVTMTKGPU : public TwoStepNVTMTK
    {
    public:
    //! Constructs the integration method and associates it with the system
    TwoStepNVTMTKGPU(std::shared_ptr<SystemDefinition> sysdef,
                     std::shared_ptr<ParticleGroup> group,
                     std::shared_ptr<ComputeThermo> thermo,
                     Scalar tau,
                     std::shared_ptr<Variant> T);
    virtual ~TwoStepNVTMTKGPU() {};

    //! Performs the first step of the integration
    virtual void integrateStepOne(uint64_t timestep);

    //! Performs the second step of the integration
    virtual void integrateStepTwo(uint64_t timestep);

    protected:
    /// Autotuner for block size (step one kernel).
    std::shared_ptr<Autotuner<1>> m_tuner_one;

    /// Autotuner for block size (step two kernel).
    std::shared_ptr<Autotuner<1>> m_tuner_two;

    /// Autotuner_angular for block size (angular step one kernel).
    std::shared_ptr<Autotuner<1>> m_tuner_angular_one;

    /// Autotuner_angular for block size (angular step two kernel).
    std::shared_ptr<Autotuner<1>> m_tuner_angular_two;
    };

    } // end namespace md
    } // end namespace hoomd

#endif // #ifndef __TWO_STEP_NVT_MTK_GPU_H__

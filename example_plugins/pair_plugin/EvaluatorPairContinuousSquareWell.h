// Copyright (c) 2009-2023 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#ifndef __PAIR_EVALUATOR_CONTINUOUSSQUAREWELL_H__
#define __PAIR_EVALUATOR_CONTINUOUSSQUAREWELL_H__

#ifndef __HIPCC__
#include <string>
#endif


#include "hoomd/HOOMDMath.h"

/*! \file EvaluatorPairContinuousSquareWell.h
    \brief Defines the pair evaluator class for the ContinuousSquareWell potential
*/

// need to declare these class methods with __device__ qualifiers when building in nvcc
// DEVICE is __host__ __device__ when included in nvcc and blank when included into the host
// compiler
#ifdef __HIPCC__
#define DEVICE __device__
#define HOSTDEVICE __host__ __device__
#else
#define DEVICE
#define HOSTDEVICE
#endif

namespace hoomd
    {
namespace md
    {

class EvaluatorPairContinuousSquareWell
    {
    public:
    //! Define the parameter type used by this pair potential evaluator
    struct param_type
        {

        //TODO: add all variable names needed for potential
        Scalar k;     //!< Spring constant
        Scalar sigma; //!< Minima of the spring

        DEVICE void load_shared(char*& ptr, unsigned int& available_bytes) { }

        HOSTDEVICE void allocate_shared(char*& ptr, unsigned int& available_bytes) const { }

#ifdef ENABLE_HIP
        //! Set CUDA memory hints
        void set_memory_hint() const
            {
            // default implementation does nothing
            }
#endif

#ifndef __HIPCC__
        param_type() : k(0), sigma(0) { }

        param_type(pybind11::dict v, bool managed = false)
            {
            //TODO: add all variables needed - this communicates the variables between python and c++
            k = v["k"].cast<Scalar>();
            sigma = v["sigma"].cast<Scalar>();
            }

        pybind11::dict asDict()
            {
            pybind11::dict v;
            //TODO: add all variables needed
            v["k"] = k;
            v["sigma"] = sigma;
            return v;
            }
#endif
        }
#if HOOMD_LONGREAL_SIZE == 32
        __attribute__((aligned(8)));
#else
        __attribute__((aligned(16)));
#endif

    //! Constructs the pair potential evaluator
    /*! \param _rsq Squared distance between the particles
        \param _rcutsq Squared distance at which the potential goes to 0
        \param _params Per type pair parameters of this potential
    */
    DEVICE EvaluatorPairExample(Scalar _rsq, Scalar _rcutsq, const param_type& _params)
    //TODO: add all variable assignements here  
        : rsq(_rsq), rcutsq(_rcutsq), k(_params.k), sigma(_params.sigma)
        {
        }

    //! Example doesn't use charge
    DEVICE static bool needsCharge()
        {
        return false;
        }
    //! Accept the optional charge value
    /*! \param qi Charge of particle i
        \param qj Charge of particle j
    */
    DEVICE void setCharge(Scalar qi, Scalar qj) { }

    //! Evaluate the force and energy
    /*! \param force_divr Output parameter to write the computed force divided by r.
        \param pair_eng Output parameter to write the computed pair energy
        \param energy_shift If true, the potential must be shifted so that
        V(r) is continuous at the cutoff
        \note There is no need to check if rsq < rcutsq in this method.
        Cutoff tests are performed in PotentialPair.

        \return True if they are evaluated or false if they are not because
        we are beyond the cutoff
    */
    DEVICE bool evalForceAndEnergy(Scalar& force_divr, Scalar& pair_eng, bool energy_shift)
        {
        // compute the force divided by r in force_divr
        if (rsq < rcutsq)
            {
            Scalar r = fast::sqrt(rsq);
            Scalar rinv = 1 / r;
           
            //TODO: change this to be force divided by r
            force_divr = k * overlap * rinv;
            //TODO: change this to be energy 
            pair_eng = Scalar(0.5) * k * overlap * overlap;

            //TODO: this is related to 'none', 'xplor', and 'shift' - look into hoomd documentation
            //TODO: to see which mode makes sense for this potential. 
            //TODO: Might only be 'none', i.e. energy_shift = false 
            if (energy_shift)
                {
                Scalar rcut = fast::sqrt(rcutsq);
                Scalar cut_overlap = sigma - rcut;
                pair_eng -= Scalar(0.5) * k * cut_overlap * cut_overlap;
                }
            return true;
            }
        else
            return false;
        }

    //! Example doesn't eval LRC integrals
    DEVICE Scalar evalPressureLRCIntegral()
        {
        return 0;
        }

    //! Example doesn't eval LRC integrals
    DEVICE Scalar evalEnergyLRCIntegral()
        {
        return 0;
        }

#ifndef __HIPCC__
    //! Get the name of this potential
    /*! \returns The potential name.
     */
    static std::string getName()
        {
        return std::string("squarewell_pair");
        }

    std::string getShapeSpec() const
        {
        throw std::runtime_error("Shape definition not supported for this pair potential.");
        }
#endif

    protected:
    Scalar rsq;    //!< Stored rsq from the constructor
    Scalar rcutsq; //!< Stored rcutsq from the constructor
    //TODO: add parameters needed for potential here 
    Scalar k;      //!< Stored k from the constructor
    Scalar sigma;  //!< Stored sigma from the constructor
    };

    }  // end namespace md
    }  // end namespace hoomd

#endif // __PAIR_EVALUATOR_CONTINUOUSSQUAREWELL_H__

// Copyright (c) 2009-2023 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#ifndef __PAIR_EVALUATOR_GCMSADJ_H__
#define __PAIR_EVALUATOR_GCMSADJ_H__

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

class EvaluatorPairGCMSAdj
    {
    public:
    //! Define the parameter type used by this pair potential evaluator
    struct param_type
        {

        //TODO: add all variable names needed for potential
        Scalar w;     //!< width of well
        Scalar sigma; //!< core size
        Scalar a; //!< well depth
        Scalar q; //!< exponent to determine well steepness


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
        param_type() : w(0), sigma(0), a(0), q(0) { }

        param_type(pybind11::dict v, bool managed = false)
            {
            //TODO: add all variables needed - this communicates the variables between python and c++
            w = v["w"].cast<Scalar>();
            sigma = v["sigma"].cast<Scalar>();
            a = v["a"].cast<Scalar>();
            q = v["q"].cast<Scalar>();
            }

        pybind11::dict asDict()
            {
            pybind11::dict v;
            //TODO: add all variables needed
            v["w"] = w;
            v["sigma"] = sigma;
            v["a"] = a;
            v["q"] = q;
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
    DEVICE EvaluatorPairGCMSAdj(Scalar _rsq, Scalar _rcutsq, const param_type& _params)
    //TODO: add all variable assignments here
        : rsq(_rsq), rcutsq(_rcutsq), w(_params.w), sigma(_params.sigma), a(_params.a), q(_params.q)
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
        if (rsq < rcutsq && a != 0 && w != 0)
            {
            Scalar r = fast::sqrt(rsq);
            Scalar rinv = 1 / r;
            Scalar w_rinv_shifted = w / (r - sigma + w);
            Scalar core_repuls = pow(w_rinv_shifted, q);
            Scalar exponent = exp(q*(r - sigma - w)/w);

            //TODO: change this to be force divided by r
            force_divr = a*rinv*(q*core_repuls*w_rinv_shifted - q*exponent/(w*(pow(1 + exponent, 2.0))));
            //TODO: change this to be energy
            pair_eng = a*(1/(1 + exponent) - core_repuls);

            //TODO: this is related to 'none', 'xplor', and 'shift' - look into hoomd documentation
            //TODO: to see which mode makes sense for this potential.
            //TODO: Might only be 'none', i.e. energy_shift = false
            if (energy_shift)
                {
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
        return std::string("gcms_adj_pair");
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
    Scalar w;      //!< Stored width of well
    Scalar sigma;  //!< Stored core size from the constructor
    Scalar a; //!< Stored depth from the constructor
    Scalar q; //?< Stored exponent from the constructor
    };

    }  // end namespace md
    }  // end namespace hoomd

#endif // __PAIR_EVALUATOR_CONTINUOUSSQUAREWELL_H__

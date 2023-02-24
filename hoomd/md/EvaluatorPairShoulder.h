// Copyright (c) 2009-2023 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#ifndef __PAIR_EVALUATOR_SHOULDER_H__
#define __PAIR_EVALUATOR_SHOULDER_H__

#ifndef __HIPCC__
#include <string>
#endif

#include "hoomd/HOOMDMath.h"

/*! \file EvaluatorPairShoulder.h
    \brief Defines the pair evaluator class for Shoulder potentials
    \details As the prototypical example of a MD pair potential, this also serves as the primary
   documentation and base reference for the implementation of pair evaluators.
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
//! Class for evaluating the Shoulder pair potential
/*! <b>General Overview</b>

    EvaluatorPairLJ is a low level computation class that computes the LJ pair potential V(r). As
   the standard MD potential, it also serves as a well documented example of how to write additional
   pair potentials. "Standard" pair potentials in hoomd are all handled via the template class
   PotentialPair. PotentialPair takes a potential evaluator as a template argument. In this way, all
   the complicated data management and other details of computing the pair force and potential on
   every single particle is only written once in the template class and the difference V(r)
   potentials that can be calculated are simply handled with various evaluator classes. Template
   instantiation is equivalent to inlining code, so there is no performance loss.

    In hoomd, a "standard" pair potential is defined as V(rsq, rcutsq, params, di, dj, qi, qj),
   where rsq is the squared distance between the two particles, rcutsq is the cutoff radius at which
   the potential goes to 0, params is any number of per type-pair parameters, di, dj are the
   diameters of particles i and j, and qi, qj are the charges of particles i and j respectively.

    Diameter and charge are not always needed by a given pair evaluator, so it must provide the
   functions needsDiameter() and needsCharge() which return boolean values signifying if they need
   those quantities or not. A false return value notifies PotentialPair that it need not even load
   those values from memory, boosting performance.

    If needsDiameter() returns true, a setDiameter(Scalar di, Scalar dj) method will be called to
   set the two diameters. Similarly, if needsCharge() returns true, a setCharge(Scalar qi, Scalar
   qj) method will be called to set the two charges.

    All other arguments are common among all pair potentials and passed into the constructor.
   Coefficients are handled in a special way: the pair evaluator class (and PotentialPair) manage
   only a single parameter variable for each type pair. Pair potentials that need more than 1
   parameter can specify that their param_type be a compound structure and reference that.

    The program flow will proceed like this: When a potential between a pair of particles is to be
   evaluated, a PairEvaluator is instantiated, passing the common parameters to the constructor and
   calling setDiameter() and/or setCharge() if need be. Then, the evalForceAndEnergy() method is
   called to evaluate the force and energy (more on that later). Thus, the evaluator must save all
   of the values it needs to compute the force and energy in member variables.

    evalForceAndEnergy() makes the necessary computations and sets the out parameters with the
   computed values. Specifically after the method completes, \a force_divr must be set to the value
    \f$ -\frac{1}{r}\frac{\partial V}{\partial r}\f$ and \a pair_eng must be set to the value \f$
   V(r) \f$ if \a energy_shift is false or \f$ V(r) - V(r_{\mathrm{cut}}) \f$ if \a energy_shift is
   true.

    A pair potential evaluator class is also used on the GPU. So all of its members must be declared
   with the DEVICE keyword before them to mark them __device__ when compiling in nvcc and blank
   otherwise. If any other code needs to diverge between the host and device (i.e., to use a special
   math function like __powf on the device), it can similarly be put inside an ifdef __HIPCC__
   block.

    <b>LJ specifics</b>

    EvaluatorPairLJ evaluates the function:
    \f[ V_{\mathrm{LJ}}(r) = 4 \varepsilon \left[ \left( \frac{\sigma}{r} \right)^{12} -
                                            \left( \frac{\sigma}{r} \right)^{6} \right] \f]
    broken up as follows for efficiency
    \f[ V_{\mathrm{LJ}}(r) = r^{-6} \cdot \left( 4 \varepsilon \sigma^{12} \cdot r^{-6} -
                                            4 \varepsilon \sigma^{6} \right) \f]
    . Similarly,
    \f[ -\frac{1}{r} \frac{\partial V_{\mathrm{LJ}}}{\partial r} = r^{-2} \cdot r^{-6} \cdot
            \left( 12 \cdot 4 \varepsilon \sigma^{12} \cdot r^{-6} - 6 \cdot 4 \varepsilon
   \sigma^{6} \right) \f]

    The LJ potential does not need diameter or charge. Two parameters are specified and stored in
    the parameter structure. It stores precomputed 4 * epsilon and sigma**6 which can be converted
    back to epsilon and sigma for the user.

    The force computation later precomputes:
    - \a lj1 = 4.0 * epsilon * pow(sigma,12.0)
    - \a lj2 = 4.0 * epsilon * pow(sigma,6.0);

*/
class EvaluatorPairShoulder
    {
    public:
    //! Define the parameter type used by this pair potential evaluator
    struct param_type
        {
        Scalar A;
        Scalar m;
        Scalar l;

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
        param_type() : A(0), m(0), l(0) { }

        param_type(pybind11::dict v, bool managed = false)
            {
            A = v["A"].cast<Scalar>();
            m = v["m"].cast<Scalar>();
            l = v["l"].cast<Scalar>();
            }

        pybind11::dict asDict()
            {
            pybind11::dict v;
            v["A"] = A;
            v["m"] = m;
            v["l"] = l;
            return v;
            }
#endif
        }
#ifdef SINGLE_PRECISION
        __attribute__((aligned(8)));
#else
        __attribute__((aligned(16)));
#endif

    //! Constructs the pair potential evaluator
    /*! \param _rsq Squared distance between the particles
        \param _rcutsq Squared distance at which the potential goes to 0
        \param _params Per type pair parameters of this potential
    */
    DEVICE EvaluatorPairShoulder(Scalar _rsq, Scalar _rcutsq, const param_type& _params)
        : rsq(_rsq), rcutsq(_rcutsq), A(_params.A), m(_params.m), l(_params.l)
        {
        }

    //! LJ doesn't use diameter
    DEVICE static bool needsDiameter()
        {
        return false;
        }
    //! Accept the optional diameter values
    /*! \param di Diameter of particle i
        \param dj Diameter of particle j
    */
    DEVICE void setDiameter(Scalar di, Scalar dj) { }

    //! LJ doesn't use charge
    DEVICE static bool needsCharge()
        {
        return false;
        }
    //! Accept the optional diameter values
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
            Scalar Exp_factor = fast::exp(-m * (r - l));
            Scalar Exp_inv = 1.0 / (1.0 + Exp_factor);

            force_divr = A * m * Exp_factor * Exp_inv * Exp_inv / r;

            pair_eng = 0.5 * A * (1 - (1.0 - Exp_factor) / (1.0 + Exp_factor));

            return true;
            }
        else
            return false;
        }

    DEVICE Scalar evalPressureLRCIntegral()
        {
        return 0;
        }

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
        return std::string("shoulder");
        }

    std::string getShapeSpec() const
        {
        throw std::runtime_error("Shape definition not supported for this pair potential.");
        }
#endif

    protected:
    Scalar rsq;    //!< Stored rsq from the constructor
    Scalar rcutsq; //!< Stored rcutsq from the constructor
    Scalar A;      //!< lj1 parameter extracted from the params passed to the constructor
    Scalar m;      //!< lj2 parameter extracted from the params passed to the constructor
    Scalar l;      //!< lj2 parameter extracted from the params passed to the constructor
    };

    } // end namespace md
    } // end namespace hoomd

#endif // __PAIR_EVALUATOR_SHOULDER_H__

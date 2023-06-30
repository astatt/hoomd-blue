// Copyright (c) 2009-2022 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#ifndef __BOND_EVALUATOR_QUARTIC_H__
#define __BOND_EVALUATOR_QUARTIC_H__

#ifndef __HIPCC__
#include <string>
#endif

#include "hoomd/HOOMDMath.h"

/*! \file EvaluatorBondQuartic.h
    \brief Defines the bond evaluator class for Quartic potentials
*/

// need to declare these class methods with __device__ qualifiers when building in nvcc
// DEVICE is __host__ __device__ when included in nvcc and blank when included into the host
// compiler
#ifdef __HIPCC__
#define DEVICE __device__
#else
#define DEVICE
#endif

namespace hoomd
    {
namespace md
    {
struct quartic_params
    {
    Scalar k;
    Scalar r_cut;
    Scalar r_1;


#ifndef __HIPCC__
    quartic_params()
        {
        k = 0;
        r_cut = 0;
        r_1 = 0;
        }

    quartic_params(pybind11::dict v)
        {
        k = v["k"].cast<Scalar>();
        r_cut = v["rcut"].cast<Scalar>();
        r_1 = v["r1"].cast<Scalar>();

        }

    pybind11::dict asDict()
        {
        pybind11::dict v;
        v["k"] = k;
        v["rcut"] = r_cut;
        v["r1"] = r_1;
        return v;
        }
#endif
    } __attribute__((aligned(16)));

//! Class for evaluating the Quartic bond potential
/*! The parameters are:
    - \a K (params.x) Stiffness parameter for the force computation
    - \a r_cut (params.y) maximum bond length for the force computation
    - \a r_1 (params.z) Value of lj1 = 4.0*epsilon*pow(sigma,12.0)
*/
class EvaluatorBondQuartic
    {
    public:
    //! Define the parameter type used by this pair potential evaluator
    typedef quartic_params param_type;

    //! Constructs the pair potential evaluator
    /*! \param _rsq Squared distance between the particles
        \param _params Per type pair parameters of this potential
    */
    DEVICE EvaluatorBondQuartic(Scalar _rsq, const param_type& _params)
        : rsq(_rsq), k(_params.k), r_cut(_params.r_cut), r_1(_params.r_1)
        {
        }

    //! This evaluator uses diameter information
    DEVICE static bool needsDiameter()
        {
        return false;
        }

    //! Accept the optional diameter values
    /*! \param da Diameter of particle a
        \param db Diameter of particle b
    */
    DEVICE void setDiameter(Scalar da, Scalar db) { }

    //! Quartic  doesn't use charge
    DEVICE static bool needsCharge()
        {
        return false;
        }

    //! Accept the optional charge values
    /*! \param qa Charge of particle a
        \param qb Charge of particle b
    */
    DEVICE void setCharge(Scalar qa, Scalar qb) { }

    //! Evaluate the force and energy
    /*! \param force_divr Output parameter to write the computed force divided by r.
        \param bond_eng Output parameter to write the computed bond energy

        \return True if they are evaluated or false if the bond
                energy is not defined
    */
    DEVICE bool evalForceAndEnergy(Scalar& force_divr, Scalar& bond_eng)
        {

        Scalar r = sqrt(rsq);

        // cutoff force and energy after r_cut to be zero
        if (rsq >= r_cut * r_cut)
        {
            force_divr = 0;
            bond_eng = 0;
            return true;
        }

        // force_divr = - 1/r *  d/dr bond_eng 
        force_divr = k*(r-r_cut)*(r-r_cut)*(3*r_1+r_cut-4*r)/r;
  
// if the result is not finite, it is likely because of a division by 0, setting force_divr to 0
// will correctly result in a 0 force in this case
#ifdef __HIPCC__
        if (!isfinite(force_divr))
#else
        if (!std::isfinite(force_divr))
#endif
            {
            force_divr = Scalar(0);
            }

        bond_eng = k*(r-r_cut)*(r-r_cut)*(r-r_cut)*(r-r_1);

        return true;
        }

#ifndef __HIPCC__
    //! Get the name of this potential
    /*! \returns The potential name.
     */
    static std::string getName()
        {
        return std::string("quartic");
        }
#endif

    protected:
    Scalar rsq;   //!< Stored rsq from the constructor
    Scalar k;     //!< K parameter
    Scalar r_cut;   //!< r_0 parameter
    Scalar r_1;   //!< r_1 parameter
    };

    } // end namespace md
    } // end namespace hoomd

#endif // __BOND_EVALUATOR_QUARTIC_H__

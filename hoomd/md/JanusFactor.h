// Copyright (c) 2009-2022 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.


#ifndef __JANUS_FACTOR_H__
#define __JANUS_FACTOR_H__

#include "hoomd/HOOMDMath.h"
#include "hoomd/VectorMath.h"

#ifdef NVCC
#define DEVICE __device__
#else
#define DEVICE
#endif

namespace hoomd
    {
namespace md
    {
/*! \file JanusFactor.h
*/
class JanusFactor
{
public:

    struct param_type
    {
        param_type()
            {
            }

        param_type(pybind11::dict params)
            : cosalpha( fast::cos(params["alpha"].cast<Scalar>()) ),
              omega(params["omega"].cast<Scalar>())
            {
            }

        pybind11::dict asDict()
            {
                pybind11::dict v;

                v["alpha"] = fast::acos(cosalpha);
                v["omega"] = omega;

                return v;
            }

        Scalar cosalpha;
        Scalar omega;
    }
#ifdef SINGLE_PRECISION
        __attribute__((aligned(8)));
#else
        __attribute__((aligned(16)));
#endif

    //!
    /*! \param

     */
    DEVICE JanusFactor(const Scalar3& _dr,
                       const Scalar4& _qi,
                       const Scalar4& _qj,
                       const Scalar& _rcutsq,
                       const param_type& _params)
        : dr(_dr), qi(_qi), qj(_qj), oi(_params.patch_orientation), oj(_q_patchj), params(_params) // TODO: this is wrong
        {
            // compute current janus direction vectors
            vec3<Scalar> ex { make_scalar3(1, 0, 0) };
            vec3<Scalar> ey { make_scalar3(0, 1, 0) };
            vec3<Scalar> ez { make_scalar3(0, 0, 1) };

            vec3<Scalar> ei, ej;
            vec3<Scalar> a1, a2, a3, b1, b2, b3;

            ei = rotate(quat<Scalar>(qi), ex);
            a1 = rotate(quat<Scalar>(qi), ex);
            a2 = rotate(quat<Scalar>(qi), ey);
            a3 = rotate(quat<Scalar>(qi), ez);
            // patch points relative to x (a1) direction of particle
            ni = rotate(quat<Scalar>(a1), oi);

            ej = rotate(quat<Scalar>(qj), ex);
            b1 = rotate(quat<Scalar>(qj), ex);
            b2 = rotate(quat<Scalar>(qj), ey);
            b3 = rotate(quat<Scalar>(qj), ez);

            // patch points relative to x (b1) direction of particle
            ni = rotate(quat<Scalar>(b1), oj);

            // compute distance
            drsq = dot(dr, dr);
            magdr = fast::sqrt(drsq);

            // cos(angle between dr and pointing vector)
            // which as first implemented is the same as the angle between the patch and pointing director
            doti = -dot(vec3<Scalar>(dr), ei) / magdr; // negative because dr = dx = pi - pj
            dotj = dot(vec3<Scalar>(dr), ej) / magdr;

            costhetai = -dot(vec3<Scalar>(dr), ni) / magdr;
            costhetaj = dot(vec3<Scalar>(dr), nj) / magdr;
        }

    DEVICE inline Scalar Modulatori()
        {
            return Scalar(1.0) / ( Scalar(1.0) + fast::exp(-params.omega*(costhetai-params.cosalpha)) );
        }

    DEVICE inline Scalar Modulatorj()
        {
            return Scalar(1.0) / ( Scalar(1.0) + fast::exp(-params.omega*(costhetaj-params.cosalpha)) );
        }

    DEVICE Scalar ModulatorPrimei() // D[f[\[Theta], \[Alpha], \[Omega]], Cos[\[Theta]]]
        {
            Scalar fact = Modulatori();
            // weird way of writing out the derivative of f with respect to doti = Cos[theta] =
            return params.omega * fast::exp(-params.omega*(costhetai-params.cosalpha)) * fact * fact;
        }

    DEVICE Scalar ModulatorPrimej()
        {
            Scalar fact = Modulatorj();
            return params.omega * fast::exp(-params.omega*(costhetaj-params.cosalpha)) * fact * fact;
        }


    Scalar3 dr;
    Scalar4 qi;
    Scalar4 qj;

    Scalar4 oi; // quaternion representing orientation of patch with respect to particle x direction (a1)
    Scalar4 oj; // quaternion representing orientation of patch with respect to particle x direction (b1)
    param_type params;

    Scalar3 ni; // pointing vector for patch on particle i
    Scalar3 nj; // pointing vector for patch on particle j

    Scalar3 a1;
    Scalar3 a2;
    Scalar3 a3;
    Scalar3 b1;
    Scalar3 b2;
    Scalar3 b3;
    Scalar drsq;
    Scalar magdr;
    Scalar doti;
    Scalar dotj;
};



    } // end namespace md
    } // end namespace hoomd

#endif // __JANUS_FACTOR_H__
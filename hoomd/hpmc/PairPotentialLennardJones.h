// Copyright (c) 2009-2023 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#pragma once

#include "hoomd/hpmc/IntegratorHPMC.h"

namespace hoomd
    {
namespace hpmc
    {

/*** Compute Lennard-Jones energy between two particles.

For use with HPMC simulations.
*/
class PairPotentialLennardJones : public hpmc::PairPotential
    {
    public:
    PairPotentialLennardJones(std::shared_ptr<SystemDefinition> sysdef);
    virtual ~PairPotentialLennardJones() { }

    virtual LongReal getRCut();
    virtual LongReal energy(const vec3<LongReal>& r_ij,
                             unsigned int type_i,
                             const quat<LongReal>& q_i,
                             LongReal charge_i,
                             unsigned int type_j,
                             const quat<LongReal>& q_j,
                             LongReal charge_j);

    virtual void setParamsPython(pybind11::tuple typ, pybind11::dict params);
    virtual pybind11::dict getParamsPython(pybind11::tuple typ);

    protected:
    /// Shifting modes that can be applied to the energy
    enum EnergyShiftMode
        {
        no_shift = 0,
        shift,
        xplor
        };

    /// Type pair parameters of LJ potential
    struct ParamType
        {
        ParamType()
            {
            sigma_6 = 0;
            epsilon_x_4 = 0;
            r_cut_squared = 0;
            r_on_squared = 0;
            mode = no_shift;
            }

        ParamType(pybind11::dict v)
            {
            auto sigma(v["sigma"].cast<LongReal>());
            auto epsilon(v["epsilon"].cast<LongReal>());
            auto r_cut(v["r_cut"].cast<LongReal>());
            auto r_on(v["r_on"].cast<LongReal>());
            auto mode_str(v["mode"].cast<std::string>());

            sigma_6 = sigma * sigma * sigma * sigma * sigma * sigma;
            epsilon_x_4 = LongReal(4.0) * epsilon;
            r_cut_squared = r_cut * r_cut;
            r_on_squared = r_on * r_on;

            if (mode_str == "none")
                {
                mode = no_shift;
                }
            else if (mode_str == "shift")
                {
                mode = shift;
                }
            else if (mode_str == "xplor")
                {
                mode = xplor;
                }
            else
                {
                throw std::domain_error("Invalid mode " + mode_str);
                }
            }

        pybind11::dict asDict()
            {
            pybind11::dict result;

            result["sigma"] = pow(sigma_6, 1. / 6.);
            result["epsilon"] = epsilon_x_4 / 4.0;
            result["r_cut"] = slow::sqrt(r_cut_squared);
            result["r_on"] = slow::sqrt(r_on_squared);
            result["mode"] = "none";

            if (mode == shift)
                {
                result["mode"] = "shift";
                }
            if (mode == xplor)
                {
                result["mode"] = "xplor";
                }

            return result;
            }

        LongReal sigma_6;
        LongReal epsilon_x_4;
        LongReal r_cut_squared;
        LongReal r_on_squared;
        EnergyShiftMode mode;
        };

    Index2DUpperTriangular m_type_param_index;
    std::vector<ParamType> m_params;
    };

namespace detail
    {
void export_PairPotentialLennardJones(pybind11::module& m);

    } // end namespace detail
    } // end namespace hpmc
    } // end namespace hoomd

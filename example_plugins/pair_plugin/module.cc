// Copyright (c) 2009-2022 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// Include the defined classes that are to be exported to python
#include "EvaluatorPairExample.h"
#include "EvaluatorBondQuartic.h"

#include "hoomd/md/PotentialPair.h"
#include "hoomd/md/PotentialBond.h"
#include <pybind11/pybind11.h>
#ifdef ENABLE_HIP
#include "hoomd/md/PotentialPairGPU.h"
#include "hoomd/md/PotentialBondGPU.h"
#endif

namespace hoomd
    {
namespace md
    {

// specify the python module. Note that the name must explicitly match the PROJECT() name provided
// in CMakeLists (with an underscore in front)
PYBIND11_MODULE(_pair_plugin, m)
    {
    detail::export_PotentialPair<EvaluatorPairExample>(m, "PotentialPairExample");
    detail::export_PotentialBond<EvaluatorBondQuartic>(m, "PotentialBondQuartic");
#ifdef ENABLE_HIP
    detail::export_PotentialPairGPU<EvaluatorPairExample>(m, "PotentialPairExampleGPU");
    detail::export_PotentialBondGPU<EvaluatorBondQuartic>(m, "PotentialBondQuarticGPU");

#endif
    }

    } // end namespace md
    } // end namespace hoomd

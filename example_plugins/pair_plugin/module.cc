// Copyright (c) 2009-2023 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// Include the defined classes that are to be exported to python
#include "EvaluatorPairExample.h"
#include "EvaluatorPairContinuousSquareWell.h"
#include "EvaluatorPairExpandedYukawa.h"
#include "EvaluatorPairGCMSAdj.h"

#include "hoomd/md/PotentialPair.h"
#include <pybind11/pybind11.h>
#ifdef ENABLE_HIP
#include "hoomd/md/PotentialPairGPU.h"
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
    detail::export_PotentialPair<EvaluatorPairContinuousSquareWell>(m, "PotentialPairContinuousSquareWell");
    detail::export_PotentialPair<EvaluatorPairExpandedYukawa>(m, "PotentialPairExpandedYukawa");
    detail::export_PotentialPair<EvaluatorPairGCMSAdj>(m, "PotentialPairGCMSAdj");

#ifdef ENABLE_HIP
    detail::export_PotentialPairGPU<EvaluatorPairExample>(m, "PotentialPairExampleGPU");
    detail::export_PotentialPairGPU<EvaluatorPairContinuousSquareWell>(m, "PotentialPairContinuousSquareWellGPU");
    detail::export_PotentialPairGPU<EvaluatorPairExpandedYukawa>(m, "PotentialPairExpandedYukawaGPU");
    detail::export_PotentialPairGPU<EvaluatorPairGCMSAdj>(m, "PotentialPairGCMSAdjGPU");

#endif
    }

    } // end namespace md
    } // end namespace hoomd

// Copyright (c) 2009-2022 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// See md/CMakeLists.txt for the source of these variables to be processed by CMake's
// configure_file().

// clang-format off
#include "hoomd/md/PotentialBondGPU.h"
#include "EvaluatorBondQuartic.h"


// clang-format on

namespace hoomd
    {
namespace md
    {

namespace detail
    {

void export_PotentialBondQuarticGPU(pybind11::module& m)
    {
    export_PotentialBondGPU<EvaluatorBondQuartic>(m, "PotentialBondQuarticGPU");
    }

    } // end namespace detail
    } // end namespace md
    } // end namespace hoomd

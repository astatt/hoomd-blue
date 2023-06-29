// Copyright (c) 2009-2022 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// See md/CMakeLists.txt for the source of these variables to be processed by CMake's
// configure_file().

// clang-format off
#include "hoomd/md/PotentialBondGPU.cuh"
#include "EvaluatorBondQuartic.h"


namespace hoomd
    {
namespace md
    {
namespace kernel
    {
template __attribute__((visibility("default"))) hipError_t
gpu_compute_bond_forces<EvaluatorBondQuartic, 2>(const kernel::bond_args_t<2>& bond_args,
                                            const typename EvaluatorBondQuartic::param_type* d_params,
                                            unsigned int* d_flags);
    } // end namespace kernel
    } // end namespace md
    } // end namespace hoomd

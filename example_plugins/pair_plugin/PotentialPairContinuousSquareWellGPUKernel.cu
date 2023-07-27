// Copyright (c) 2009-2023 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include "EvaluatorPairContinuousSquareWell.h"
#include "hoomd/md/PotentialPairGPU.cuh"

namespace hoomd
    {
namespace md
    {
namespace kernel
    {
template __attribute__((visibility("default"))) hipError_t
gpu_compute_pair_forces<EvaluatorPairExample>(const pair_args_t& pair_args,
                                              const EvaluatorPairContinuousSquareWell::param_type* d_params);
    } // end namespace kernel
    } // end namespace md
    } // end namespace hoomd

# Copyright (c) 2009-2023 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Example pair potential."""

# Import the C++ module.
from hoomd.pair_plugin import _pair_plugin

# Impot the hoomd Python package and other necessary components.
from hoomd.md import pair
from hoomd.data.parameterdicts import TypeParameterDict
from hoomd.data.typeparam import TypeParameter


class ExamplePair(pair.Pair):
    """Example pair potential."""

    # set static class data
    _ext_module = _pair_plugin
    _cpp_class_name = "PotentialPairExample"
    _accepted_modes = ("none", "shift", "xplor")

    def __init__(self, nlist, default_r_cut=None, default_r_on=0., mode='none'):
        super().__init__(nlist, default_r_cut, default_r_on, mode)
        params = TypeParameter(
            'params', 'particle_types',
            TypeParameterDict(k=float, sigma=float, len_keys=2))
        self._add_typeparam(params)

class ContinuousSquareWellPair(pair.Pair):
    """ContinuousSquareWell pair potential.
    
    #TODO: write proper documentation for this - see hoomd documentation for examples 

    """

    # set static class data
    _ext_module = _pair_plugin
    _cpp_class_name = "PotentialPairContinuousSquareWell"

    #TODO: check that these modes make sense for potential - I think only none might make sense?
    #TODO: check in hoomd documentation for what these do

    _accepted_modes = ("none", "shift")

    def __init__(self, nlist, default_r_cut=None, default_r_on=0., mode='none'):
        super().__init__(nlist, default_r_cut, default_r_on, mode)
        params = TypeParameter(
            'params', 'particle_types',
            TypeParameterDict(n=float, m=float, lambda_val=float, A=float, len_keys=2))
        self._add_typeparam(params)


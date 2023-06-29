# Copyright (c) 2009-2022 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Example pair potential."""

# Import the C++ module.
from hoomd.pair_plugin import _pair_plugin

# Import the hoomd Python package and other necessary components.
from hoomd.md import bond
from hoomd.data.parameterdicts import TypeParameterDict
from hoomd.data.typeparam import TypeParameter


class Quartic(bond.Bond):
    """Example bond potential."""

    # set static class data
    _ext_module = _pair_plugin
    _cpp_class_name = "PotentialBondQuartic"
  
    def __init__(self):
        super().__init__()
        params = TypeParameter("params", "bond_types",
                               TypeParameterDict(k=float, rcut=float, r1=float, len_keys=1))
        self._add_typeparam(params)



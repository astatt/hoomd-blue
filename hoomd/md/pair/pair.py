# Copyright (c) 2009-2021 The Regents of the University of Michigan
# This file is part of the HOOMD-blue project, released under the BSD 3-Clause
# License.

"""Pair potentials."""

import hoomd
from hoomd.md import _md
from hoomd.md import force
from hoomd.md.nlist import NList
from hoomd.data.parameterdicts import ParameterDict, TypeParameterDict
from hoomd.data.typeparam import TypeParameter
from hoomd.data.typeconverter import (OnlyFrom, OnlyTypes, positive_real,
                                      nonnegative_real)

validate_nlist = OnlyTypes(NList)


class Pair(force.Force):
    r"""Common pair potential documentation.

    Users should not invoke `Pair` directly. It is a base command
    that provides common features to all standard pair forces. Common
    documentation for all pair potentials is documented here.

    All pair force commands specify that a given potential energy and force be
    computed on all non-excluded particle pairs in the system within a short
    range cutoff distance :math:`r_{\mathrm{cut}}`.

    The force :math:`\vec{F}` applied between each pair of particles is:

    .. math::
        :nowrap:

        \begin{eqnarray*}
        \vec{F}  = & -\nabla V(r) & r < r_{\mathrm{cut}} \\
                  = & 0           & r \ge r_{\mathrm{cut}} \\
        \end{eqnarray*}

    where :math:`\vec{r}` is the vector pointing from one particle to the other
    in the pair, and :math:`V(r)` is chosen by a mode switch:

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V(r)  = & V_{\mathrm{pair}}(r) & \mathrm{mode\ is\ no\_shift} \\
              = & V_{\mathrm{pair}}(r) - V_{\mathrm{pair}}(r_{\mathrm{cut}})
              & \mathrm{mode\ is\ shift} \\
              = & S(r) \cdot V_{\mathrm{pair}}(r) & \mathrm{mode\ is\
              xplor\ and\ } r_{\mathrm{on}} < r_{\mathrm{cut}} \\
              = & V_{\mathrm{pair}}(r) - V_{\mathrm{pair}}(r_{\mathrm{cut}})
              & \mathrm{mode\ is\ xplor\ and\ } r_{\mathrm{on}} \ge
              r_{\mathrm{cut}}
        \end{eqnarray*}

    :math:`S(r)` is the XPLOR smoothing function:

    .. math::
        :nowrap:

        \begin{eqnarray*}
        S(r) = & 1 & r < r_{\mathrm{on}} \\
             = & \frac{(r_{\mathrm{cut}}^2 - r^2)^2 \cdot
             (r_{\mathrm{cut}}^2 + 2r^2 -
             3r_{\mathrm{on}}^2)}{(r_{\mathrm{cut}}^2 -
             r_{\mathrm{on}}^2)^3}
               & r_{\mathrm{on}} \le r \le r_{\mathrm{cut}} \\
             = & 0 & r > r_{\mathrm{cut}} \\
         \end{eqnarray*}

    and :math:`V_{\mathrm{pair}}(r)` is the specific pair potential chosen by
    the respective command.

    Enabling the XPLOR smoothing function :math:`S(r)` results in both the
    potential energy and the force going smoothly to 0 at :math:`r =
    r_{\mathrm{cut}}`, reducing the rate of energy drift in long simulations.
    :math:`r_{\mathrm{on}}` controls the point at which the smoothing starts,
    so it can be set to only slightly modify the tail of the potential. It is
    suggested that you plot your potentials with various values of
    :math:`r_{\mathrm{on}}` in order to find a good balance between a smooth
    potential function and minimal modification of the original
    :math:`V_{\mathrm{pair}}(r)`. A good value for the LJ potential is
    :math:`r_{\mathrm{on}} = 2 \cdot \sigma`.

    The split smoothing / shifting of the potential when the mode is ``xplor``
    is designed for use in mixed WCA / LJ systems. The WCA potential and it's
    first derivative already go smoothly to 0 at the cutoff, so there is no need
    to apply the smoothing function. In such mixed systems, set
    :math:`r_{\mathrm{on}}` to a value greater than :math:`r_{\mathrm{cut}}`
    for those pairs that interact via WCA in order to enable shifting of the WCA
    potential to 0 at the cutoff.

    The following coefficients must be set per unique pair of particle types.
    See `hoomd.md.pair` for information on how to set coefficients.


    .. py:attribute:: r_cut

        *r_cut* (in distance units), *optional*: defaults to the value ``r_cut``
        specified on construction.

        Type: `TypeParameter` [`tuple` [``particle_type``, ``particle_type``],
        `float`])

    .. py:attribute:: r_on

        *r_on* (in distance units),  *optional*: defaults to the value ``r_on``
        specified on construction

        Type: `TypeParameter` [`tuple` [``particle_type``, ``particle_type``],
        `float`])
    """

    def __init__(self, nlist, r_cut=None, r_on=0., mode='none'):
        self._nlist = validate_nlist(nlist)
        tp_r_cut = TypeParameter('r_cut', 'particle_types',
                                 TypeParameterDict(positive_real, len_keys=2))
        if r_cut is not None:
            tp_r_cut.default = r_cut
        tp_r_on = TypeParameter('r_on', 'particle_types',
                                TypeParameterDict(nonnegative_real, len_keys=2))
        if r_on is not None:
            tp_r_on.default = r_on
        self._extend_typeparam([tp_r_cut, tp_r_on])
        self._param_dict.update(
            ParameterDict(mode=OnlyFrom(['none', 'shift', 'xplor'])))
        self.mode = mode

    def compute_energy(self, tags1, tags2):
        r"""Compute the energy between two sets of particles.

        Args:
            tags1 (``ndarray<int32>``): a numpy array of particle tags in the
                first group
            tags2 (``ndarray<int32>``): a numpy array of particle tags in the
                second group

        .. math::

            U = \sum_{i \in \mathrm{tags1}, j \in \mathrm{tags2}} V_{ij}(r)

        where :math:`V_{ij}(r)` is the pairwise energy between two particles
        :math:`i` and :math:`j`.

        Assumed properties of the sets *tags1* and *tags2* are:

        - *tags1* and *tags2* are disjoint
        - all elements in *tags1* and *tags2* are unique
        - *tags1* and *tags2* are contiguous numpy arrays of dtype int32

        None of these properties are validated.

        Examples::

            tags=numpy.linspace(0,N-1,1, dtype=numpy.int32)
            # computes the energy between even and odd particles
            U = mypair.compute_energy(tags1=numpy.array(tags[0:N:2]),
                                      tags2=numpy.array(tags[1:N:2]))

        """
        # TODO future versions could use np functions to test the assumptions
        # above and raise an error if they occur.
        return self._cpp_obj.computeEnergyBetweenSets(tags1, tags2)

    def _attach(self):
        # create the c++ mirror class
        if not self._nlist._added:
            self._nlist._add(self._simulation)
        else:
            if self._simulation != self._nlist._simulation:
                raise RuntimeError("{} object's neighbor list is used in a "
                                   "different simulation.".format(type(self)))
        if not self.nlist._attached:
            self.nlist._attach()
        if isinstance(self._simulation.device, hoomd.device.CPU):
            cls = getattr(_md, self._cpp_class_name)
            self.nlist._cpp_obj.setStorageMode(
                _md.NeighborList.storageMode.half)
        else:
            cls = getattr(_md, self._cpp_class_name + "GPU")
            self.nlist._cpp_obj.setStorageMode(
                _md.NeighborList.storageMode.full)
        self._cpp_obj = cls(self._simulation.state._cpp_sys_def,
                            self.nlist._cpp_obj)

        super()._attach()

    @property
    def nlist(self):
        """Neighbor list used to compute the pair potential."""
        return self._nlist

    @nlist.setter
    def nlist(self, value):
        if self._attached:
            raise RuntimeError("nlist cannot be set after scheduling.")
        else:
            self._nlist = validate_nlist(value)

    @property
    def _children(self):
        return [self.nlist]


class LJ(Pair):
    r"""Lennard-Jones pair potential.

    Args:
        nlist (`hoomd.md.nlist.NList`): Neighbor list
        r_cut (float): Default cutoff radius (in distance units).
        r_on (float): Default turn-on radius (in distance units).
        mode (str): energy shifting/smoothing mode

    `LJ` specifies that a Lennard-Jones pair potential should be
    applied between every non-excluded particle pair in the simulation.

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{LJ}}(r)  = & 4 \varepsilon \left[ \left(
        \frac{\sigma}{r} \right)^{12} - \left( \frac{\sigma}{r}
        \right)^{6} \right] & r < r_{\mathrm{cut}} \\
        = & 0 & r \ge r_{\mathrm{cut}} \\
        \end{eqnarray*}

    See `Pair` for details on how forces are calculated and the
    available energy shifting and smoothing modes.  Use `params` dictionary
    to set potential coefficients. The coefficients must be set per
    unique pair of particle types.

    .. py:attribute:: params

        The LJ potential parameters. The dictionary has the following keys:

        * ``epsilon`` (`float`, **required**) -
          energy parameter :math:`\varepsilon` (in energy units)
        * ``sigma`` (`float`, **required**) -
          particle size :math:`\sigma` (in distance units)

        Type: `TypeParameter` [`tuple` [``particle_type``, ``particle_type``],
        `dict`]

    Example::

        nl = nlist.Cell()
        lj = pair.LJ(nl, r_cut=3.0)
        lj.params[('A', 'A')] = {'sigma': 1.0, 'epsilon': 1.0}
        lj.r_cut[('A', 'B')] = 3.0
    """
    _cpp_class_name = "PotentialPairLJ"

    def __init__(self, nlist, r_cut=None, r_on=0., mode='none'):
        super().__init__(nlist, r_cut, r_on, mode)
        params = TypeParameter(
            'params', 'particle_types',
            TypeParameterDict(epsilon=float, sigma=float, len_keys=2))
        self._add_typeparam(params)


class Gauss(Pair):
    r"""Gaussian pair potential.

    Args:
        nlist (`hoomd.md.nlist.NList`): Neighbor list
        r_cut (float): Default cutoff radius (in distance units).
        r_on (float): Default turn-on radius (in distance units).
        mode (str): energy shifting/smoothing mode.

    `Gauss` specifies that a Gaussian pair potential should be applied
    between every non-excluded particle pair in the simulation.

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{gauss}}(r)  = & \varepsilon \exp \left[ -\frac{1}{2}
                                  \left( \frac{r}{\sigma} \right)^2 \right]
                                  & r < r_{\mathrm{cut}} \\
                                 = & 0 & r \ge r_{\mathrm{cut}} \\
        \end{eqnarray*}

    See `Pair` for details on how forces are calculated and the
    available energy shifting and smoothing modes. Use `params` dictionary to
    set potential coefficients. The coefficients must be set per unique pair of
    particle types.

    .. py:attribute:: params

        The Gauss potential parameters. The dictionary has the following
        keys:

        * ``epsilon`` (`float`, **required**) - energy parameter
          :math:`\varepsilon` (in energy units)
        * ``sigma`` (`float`, **required**) - particle size :math:`\sigma`
          (in distance units)

        Type: `TypeParameter` [`tuple` [``particle_type``, ``particle_type``],
        `dict`]

    Example::

        nl = nlist.Cell()
        gauss = pair.Gauss(r_cut=3.0, nlist=nl)
        gauss.params[('A', 'A')] = dict(epsilon=1.0, sigma=1.0)
        gauss.r_cut[('A', 'B')] = 3.0
    """
    _cpp_class_name = "PotentialPairGauss"

    def __init__(self, nlist, r_cut=None, r_on=0., mode='none'):
        super().__init__(nlist, r_cut, r_on, mode)
        params = TypeParameter(
            'params', 'particle_types',
            TypeParameterDict(epsilon=float, sigma=float, len_keys=2))
        self._add_typeparam(params)


class SLJ(Pair):
    r"""Shifted Lennard-Jones pair potential.

    Args:
        nlist (`hoomd.md.nlist.NList`): Neighbor list
        r_cut (float): Default cutoff radius (in distance units).
        r_on (float): Default turn-on radius (in distance units).
        mode (str): Energy shifting/smoothing mode

    `SLJ` specifies that a shifted Lennard-Jones type pair potential
    should be applied between every non-excluded particle pair in the
    simulation.

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{SLJ}}(r)  = & 4 \varepsilon \left[ \left(
                                \frac{\sigma}{r - \Delta} \right)^{12} -
                                \left( \frac{\sigma}{r - \Delta}
                                \right)^{6} \right] & r < (r_{\mathrm{cut}}
                                + \Delta) \\
                             = & 0 & r \ge (r_{\mathrm{cut}} + \Delta) \\
        \end{eqnarray*}

    where :math:`\Delta = (d_i + d_j)/2 - 1` and :math:`d_i` is the diameter of
    particle :math:`i`.

    See `Pair` for details on how forces are calculated and the
    available energy shifting and smoothing modes. Use `params` dictionary to
    set potential coefficients. The coefficients must be set per unique pair of
    particle types.

    Attention:
        Due to the way that `SLJ` modifies the cutoff criteria, a shift_mode
        of *xplor* is not supported.

    Set the ``max_diameter`` property of the neighbor list object to the largest
    particle diameter in the system (where **diameter** is a per-particle
    property of the same name in `hoomd.State`).

    Warning:
        Failure to set ``max_diameter`` will result in missing pair
        interactions.

    .. py:attribute:: params

        The potential parameters. The dictionary has the following keys:

        * ``epsilon`` (`float`, **required**) - energy parameter
          :math:`\varepsilon` (in energy units)
        * ``sigma`` (`float`, **required**) - particle size :math:`\sigma`
          (in distance units)

        Type: `TypeParameter` [`tuple` [``particle_type``, ``particle_type``],
        `dict`]

    Example::

        nl = nlist.Cell()
        nl.max_diameter = 2.0
        slj = pair.SLJ(r_cut=3.0, nlist=nl)
        slj.params[('A', 'B')] = dict(epsilon=2.0, r_cut=3.0)
        slj.r_cut[('B', 'B')] = 2**(1.0/6.0)
    """
    _cpp_class_name = 'PotentialPairSLJ'

    def __init__(self, nlist, r_cut=None, r_on=0., mode='none'):
        if mode == 'xplor':
            raise ValueError("xplor is not a valid mode for SLJ potential")

        super().__init__(nlist, r_cut, r_on, mode)
        params = TypeParameter(
            'params', 'particle_types',
            TypeParameterDict(epsilon=float, sigma=float, len_keys=2))
        self._add_typeparam(params)

        # mode not allowed to be xplor, so re-do param dict entry without that
        # option
        param_dict = ParameterDict(mode=OnlyFrom(['none', 'shift']))
        self._param_dict.update(param_dict)
        self.mode = mode

        # this potential needs diameter shifting on
        self._nlist.diameter_shift = True

        # NOTE do we need something to automatically set the max_diameter
        # correctly?


class Yukawa(Pair):
    r"""Yukawa pair potential.

    Args:
        nlist (`hoomd.md.nlist.NList`): Neighbor list
        r_cut (float): Default cutoff radius (in distance units).
        r_on (float): Default turn-on radius (in distance units).
        mode (str): Energy shifting mode.

    `Yukawa` specifies that a Yukawa pair potential should be applied between
    every non-excluded particle pair in the simulation.

    .. math::
        :nowrap:

        \begin{eqnarray*}
          V_{\mathrm{yukawa}}(r) = & \varepsilon \frac{ \exp \left(
          -\kappa r \right) }{r} & r < r_{\mathrm{cut}} \\
                                  = & 0 & r \ge r_{\mathrm{cut}} \\
        \end{eqnarray*}

    See `Pair` for details on how forces are calculated and the available
    energy shifting and smoothing modes. Use `params` dictionary to set
    potential coefficients. The coefficients must be set per unique pair of
    particle types.

    .. py:attribute:: params

        The Yukawa potential parameters. The dictionary has the following
        keys:

        * ``epsilon`` (`float`, **required**) - energy parameter
          :math:`\varepsilon` (in energy units)
        * ``kappa`` (`float`, **required**) - scaling parameter
          :math:`\kappa` (in units of 1/distance)

        Type: `TypeParameter` [`tuple` [``particle_type``, ``particle_type``],
        `dict`]

    Example::

        nl = nlist.Cell()
        yukawa = pair.Yukawa(r_cut=3.0, nlist=nl)
        yukawa.params[('A', 'A')] = dict(epsilon=1.0, kappa=1.0)
        yukawa.r_cut[('A', 'B')] = 3.0
    """
    _cpp_class_name = "PotentialPairYukawa"

    def __init__(self, nlist, r_cut=None, r_on=0., mode='none'):
        super().__init__(nlist, r_cut, r_on, mode)
        params = TypeParameter(
            'params', 'particle_types',
            TypeParameterDict(kappa=float, epsilon=float, len_keys=2))
        self._add_typeparam(params)


class Ewald(Pair):
    r"""Ewald pair potential.

    Args:
        nlist (`hoomd.md.nlist.NList`): Neighbor list
        r_cut (float): Default cutoff radius (in distance units).
        r_on (float): Default turn-on radius (in distance units).
        mode (str): Energy shifting mode.

    `Ewald` specifies that a Ewald pair potential should be applied between
    every non-excluded particle pair in the simulation.

    .. math::
        :nowrap:

        \begin{eqnarray*}
         V_{\mathrm{ewald}}(r)  = & q_i q_j \left[\mathrm{erfc}\left(\kappa
                                    r + \frac{\alpha}{2\kappa}\right)
                                    \exp(\alpha r) \\
                                    + \mathrm{erfc}\left(\kappa r -
                                    \frac{\alpha}{2 \kappa}\right)
                                    \exp(-\alpha r)\right]
                                    & r < r_{\mathrm{cut}} \\
                            = & 0 & r \ge r_{\mathrm{cut}} \\
        \end{eqnarray*}

    The Ewald potential is designed to be used in conjunction with PPPM.

    See `Pair` for details on how forces are calculated and the available
    energy shifting and smoothing modes. Use the `params` dictionary to set
    potential coefficients. The coefficients must be set per unique pair of
    particle types.

    .. py:attribute:: params

        The Ewald potential parameters. The dictionary has the following keys:

        * ``kappa`` (`float`, **required**) - Splitting parameter
          :math:`\kappa` (in units of 1/distance)
        * ``alpha`` (`float`, **required**) - Debye screening length
          :math:`\alpha` (in units of 1/distance)

        Type: `TypeParameter` [`tuple` [``particle_type``, ``particle_type``],
        `dict`]

    Example::

        nl = nlist.Cell()
        ewald = pair.Ewald(r_cut=3.0, nlist=nl)
        ewald.params[('A', 'A')] = dict(kappa=1.0, alpha=1.5)
        ewald.r_cut[('A', 'B')] = 3.0
    """
    _cpp_class_name = "PotentialPairEwald"

    def __init__(self, nlist, r_cut=None, r_on=0., mode='none'):
        super().__init__(nlist, r_cut, r_on, mode)
        params = TypeParameter(
            'params', 'particle_types',
            TypeParameterDict(kappa=float, alpha=0.0, len_keys=2))
        self._add_typeparam(params)


class Morse(Pair):
    r"""Morse pair potential.

    Args:
        nlist (`hoomd.md.nlist.NList`): Neighbor list
        r_cut (float): Default cutoff radius (in distance units).
        r_on (float): Default turn-on radius (in distance units).
        mode (str): energy shifting/smoothing mode.

    `Morse` specifies that a Morse pair potential should be applied between
    every non-excluded particle pair in the simulation.

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{morse}}(r) = & D_0 \left[ \exp \left(-2\alpha\left(
            r-r_0\right)\right) -2\exp \left(-\alpha\left(r-r_0\right)
            \right) \right] & r < r_{\mathrm{cut}} \\
            = & 0 & r \ge r_{\mathrm{cut}} \\
        \end{eqnarray*}

    See `Pair` for details on how forces are calculated and the available
    energy shifting and smoothing modes. Use `params` dictionary to set
    potential coefficients. The coefficients must be set per unique pair of
    particle types.

    .. py:attribute:: params

        The potential parameters. The dictionary has the following keys:

        * ``D0`` (`float`, **required**) - depth of the potential at its
          minimum :math:`D_0` (in energy units)
        * ``alpha`` (`float`, **required**) - the width of the potential well
          :math:`\alpha` (in units of 1/distance)
        * ``r0`` (`float`, **required**) - position of the minimum
          :math:`r_0` (in distance units)

        Type: `TypeParameter` [`tuple` [``particle_type``, ``particle_type``],
        `dict`]

    Example::

        nl = nlist.Cell()
        morse = pair.Morse(r_cut=3.0, nlist=nl)
        morse.params[('A', 'A')] = dict(D0=1.0, alpha=3.0, r0=1.0)
        morse.r_cut[('A', 'B')] = 3.0
    """

    _cpp_class_name = "PotentialPairMorse"

    def __init__(self, nlist, r_cut=None, r_on=0., mode='none'):
        super().__init__(nlist, r_cut, r_on, mode)
        params = TypeParameter(
            'params', 'particle_types',
            TypeParameterDict(D0=float, alpha=float, r0=float, len_keys=2))
        self._add_typeparam(params)


class DPD(Pair):
    r"""Dissipative Particle Dynamics.

    Args:
        nlist (`hoomd.md.nlist.NList`): Neighbor list
        kT (`hoomd.variant` or `float`): Temperature of
          thermostat (in energy units).
        r_cut (float): Default cutoff radius (in distance units).
        r_on (float): Default turn-on radius (in distance units).

    `DPD` specifies that a DPD pair force should be applied between every
    non-excluded particle pair in the simulation, including an interaction
    potential, pairwise drag force, and pairwise random force. See `Groot and
    Warren 1997 <http://dx.doi.org/10.1063/1.474784>`_.

    .. math::
        :nowrap:

        \begin{eqnarray*}
        F = F_{\mathrm{C}}(r) + F_{\mathrm{R,ij}}(r_{ij}) +
        F_{\mathrm{D,ij}}(v_{ij}) \\
        \end{eqnarray*}

    .. math::
        :nowrap:

        \begin{eqnarray*}
        F_{\mathrm{C}}(r) = & A \cdot  w(r_{ij}) \\
        F_{\mathrm{R, ij}}(r_{ij}) = & - \theta_{ij}\sqrt{3}
        \sqrt{\frac{2k_b\gamma T}{\Delta t}}\cdot w(r_{ij})  \\
        F_{\mathrm{D, ij}}(r_{ij}) = & - \gamma w^2(r_{ij})\left(
        \hat r_{ij} \circ v_{ij} \right)  \\
        \end{eqnarray*}

    .. math::
        :nowrap:

        \begin{eqnarray*}
        w(r_{ij}) = &\left( 1 - r/r_{\mathrm{cut}} \right)
        & r < r_{\mathrm{cut}} \\
                  = & 0 & r \ge r_{\mathrm{cut}} \\
        \end{eqnarray*}

    where :math:`\hat r_{ij}` is a normalized vector from particle i to
    particle j, :math:`v_{ij} = v_i - v_j`, and :math:`\theta_{ij}` is a
    uniformly distributed random number in the range [-1, 1].

    `C. L. Phillips et. al. 2011 <http://dx.doi.org/10.1016/j.jcp.2011.05.021>`_
    describes the DPD implementation details in HOOMD-blue. Cite it if you
    utilize the DPD functionality in your work.

    `DPD` does not implement and energy shift / smoothing modes due to the
    function of the force. Use `params` dictionary to set potential
    coefficients. The coefficients must be set per unique pair of particle
    types.

    To use the DPD thermostat, an `hoomd.md.methods.NVE` integrator
    must be applied to the system and the user must specify a temperature.  Use
    of the dpd thermostat pair force with other integrators will result in
    unphysical behavior. To use pair.dpd with a different conservative potential
    than :math:`F_C`, set A to zero and define the conservative pair potential
    separately.  Note that DPD thermostats are often defined in terms of
    :math:`\sigma` where :math:`\sigma = \sqrt{2k_b\gamma T}`.

    .. py:attribute:: params

        The force parameters. The dictionary has the following keys:

        * ``A`` (`float`, **required**) - :math:`A` (in force units)
        * ``gamma`` (`float`, **required**) - :math:`\gamma` (in units of
          force/velocity)

        Type: `TypeParameter` [`tuple` [``particle_type``, ``particle_type``],
        `dict`]

    Example::

        nl = nlist.Cell()
        dpd = pair.DPD(nlist=nl, kT=1.0, r_cut=1.0)
        dpd.params[('A', 'A')] = dict(A=25.0, gamma=4.5)
        dpd.params[('A', 'B')] = dict(A=40.0, gamma=4.5)
        dpd.params[('B', 'B')] = dict(A=25.0, gamma=4.5)
        dpd.params[(['A', 'B'], ['C', 'D'])] = dict(A=40.0, gamma=4.5)
    """
    _cpp_class_name = "PotentialPairDPDThermoDPD"

    def __init__(self, nlist, kT, r_cut=None, r_on=0., mode='none'):
        super().__init__(nlist, r_cut, r_on, mode)
        params = TypeParameter(
            'params', 'particle_types',
            TypeParameterDict(A=float, gamma=float, len_keys=2))
        self._add_typeparam(params)

        d = ParameterDict(kT=hoomd.variant.Variant)
        self._param_dict.update(d)

        self.kT = kT

    def _add(self, simulation):
        """Add the operation to a simulation.

        DPD uses RNGs. Warn the user if they did not set the seed.
        """
        if simulation is not None:
            simulation._warn_if_seed_unset()

        super()._add(simulation)


class DPDConservative(Pair):
    r"""DPD Conservative pair force.

    Args:
        nlist (`hoomd.md.nlist.NList`): Neighbor list
        r_cut (float): Default cutoff radius (in distance units).
        r_on (float): Default turn-on radius (in distance units).

    `DPDConservative` specifies the conservative part of the DPD pair potential
    should be applied between every non-excluded particle pair in the
    simulation. No thermostat (e.g. Drag Force and Random Force) is applied, as
    is in `DPD`.

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{DPD-C}}(r) = & A \cdot \left( r_{\mathrm{cut}} - r
          \right) - \frac{1}{2} \cdot \frac{A}{r_{\mathrm{cut}}} \cdot
          \left(r_{\mathrm{cut}}^2 - r^2 \right)
          & r < r_{\mathrm{cut}} \\
                              = & 0 & r \ge r_{\mathrm{cut}} \\
        \end{eqnarray*}


    `DPDConservative` does not implement and energy shift / smoothing modes due
    to the function of the force. Use `params` dictionary to set potential
    coefficients. The coefficients must be set per unique pair of particle
    types.

    .. py:attribute:: params

        The potential parameters. The dictionary has the following keys:

        * ``A`` (`float`, **required**) - :math:`A` (in force units)

        Type: `TypeParameter` [`tuple` [``particle_type``, ``particle_type``],
        `dict`]

    Example::

        nl = nlist.Cell()
        dpdc = pair.DPDConservative(nlist=nl, r_cut=3.0)
        dpdc.params[('A', 'A')] = dict(A=1.0)
        dpdc.params[('A', 'B')] = dict(A=2.0, r_cut = 1.0)
        dpdc.params[(['A', 'B'], ['C', 'D'])] = dict(A=3.0)
    """
    _cpp_class_name = "PotentialPairDPD"

    def __init__(self, nlist, r_cut=None, r_on=0., mode='none'):
        # initialize the base class
        super().__init__(nlist, r_cut, r_on, mode)
        params = TypeParameter('params', 'particle_types',
                               TypeParameterDict(A=float, len_keys=2))
        self._add_typeparam(params)


class DPDLJ(Pair):
    r"""Dissipative Particle Dynamics with a LJ conservative force.

    Args:
        nlist (`hoomd.md.nlist.NList`): Neighbor list
        kT (`hoomd.variant` or `float`): Temperature of
            thermostat (in energy units).
        r_cut (float): Default cutoff radius (in distance units).
        r_on (float): Default turn-on radius (in distance units).

    `DPDLJ` specifies that a DPD thermostat and a Lennard-Jones pair potential
    should be applied between every non-excluded particle pair in the
    simulation.

    `C. L. Phillips et. al. 2011 <http://dx.doi.org/10.1016/j.jcp.2011.05.021>`_
    describes the DPD implementation details in HOOMD-blue. Cite it if you
    utilize the DPD functionality in your work.

    .. math::
        :nowrap:

        \begin{eqnarray*}
        F = F_{\mathrm{C}}(r) + F_{\mathrm{R,ij}}(r_{ij}) +
            F_{\mathrm{D,ij}}(v_{ij}) \\
        \end{eqnarray*}

    .. math::
        :nowrap:

        \begin{eqnarray*}
        F_{\mathrm{C}}(r) = & \partial V_{\mathrm{LJ}} / \partial r \\
        F_{\mathrm{R, ij}}(r_{ij}) = & - \theta_{ij}\sqrt{3}
            \sqrt{\frac{2k_b\gamma T}{\Delta t}}\cdot w(r_{ij})  \\
        F_{\mathrm{D, ij}}(r_{ij}) = & - \gamma w^2(r_{ij})
            \left( \hat r_{ij} \circ v_{ij} \right)  \\
        \end{eqnarray*}

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{LJ}}(r) = & 4 \varepsilon \left[ \left(
            \frac{\sigma}{r} \right)^{12} -
             \left( \frac{\sigma}{r} \right)^{6} \right]
            & r < r_{\mathrm{cut}} \\
                            = & 0 & r \ge r_{\mathrm{cut}} \\
        \end{eqnarray*}

    .. math::
        :nowrap:

        \begin{eqnarray*}
        w(r_{ij}) = &\left( 1 - r/r_{\mathrm{cut}} \right)
            & r < r_{\mathrm{cut}} \\
                  = & 0 & r \ge r_{\mathrm{cut}} \\
        \end{eqnarray*}

    where :math:`\hat r_{ij}` is a normalized vector from particle i to
    particle j, :math:`v_{ij} = v_i - v_j`, and :math:`\theta_{ij}` is a
    uniformly distributed random number in the range [-1, 1].

    Use `params` dictionary to set potential coefficients. The coefficients must
    be set per unique pair of particle types.

    To use the DPD thermostat, an `hoomd.md.methods.NVE` integrator
    must be applied to the system and the user must specify a temperature.  Use
    of the dpd thermostat pair force with other integrators will result in
    unphysical behavior.

    .. py:attribute:: params

        The DPDLJ potential parameters. The dictionary has the following keys:

        * ``epsilon`` (`float`, **required**) - :math:`\varepsilon`
          (in energy units)
        * ``sigma`` (`float`, **required**) - :math:`\sigma`
          (in distance units)
        * ``gamma`` (`float`, **required**) - :math:`\gamma` (in units of
          force/velocity)

        Type: `TypeParameter` [`tuple` [``particle_type``, ``particle_type``],
        `dict`]

    Example::

        nl = nlist.Cell()
        dpdlj = pair.DPDLJ(nlist=nl, kT=1.0, r_cut=2.5)
        dpdlj.params[('A', 'A')] = dict(epsilon=1.0, sigma=1.0, gamma=4.5)
        dpdlj.params[(['A', 'B'], ['C', 'D'])] = dict(
            epsilon=3.0, sigma=1.0, gamma=1.2)
        dpdlj.r_cut[('B', 'B')] = 2.0**(1.0/6.0)
    """
    _cpp_class_name = "PotentialPairDPDLJThermoDPD"

    def __init__(self, nlist, kT, r_cut=None, r_on=0., mode='none'):
        if mode == 'xplor':
            raise ValueError("xplor smoothing is not supported with pair.DPDLJ")

        super().__init__(nlist, r_cut, r_on, mode)
        params = TypeParameter(
            'params', 'particle_types',
            TypeParameterDict(epsilon=float,
                              sigma=float,
                              gamma=float,
                              len_keys=2))
        self._add_typeparam(params)

        d = ParameterDict(kT=hoomd.variant.Variant,
                          mode=OnlyFrom(['none', 'shift']))
        self._param_dict.update(d)

        self.kT = kT
        self.mode = mode

    def _add(self, simulation):
        """Add the operation to a simulation.

        DPDLJ uses RNGs. Warn the user if they did not set the seed.
        """
        if simulation is not None:
            simulation._warn_if_seed_unset()

        super()._add(simulation)


class ForceShiftedLJ(Pair):
    r"""Force-shifted Lennard-Jones pair potential.

    Args:
        nlist (`hoomd.md.nlist.NList`): Neighbor list
        r_cut (float): Default cutoff radius (in distance units).
        r_on (float): Default turn-on radius (in distance units).
        mode (str): energy shifting/smoothing mode.

    `ForceShiftedLJ` specifies that a modified Lennard-Jones pair force should
    be applied between non-excluded particle pair in the simulation. The force
    differs from the one calculated by  `LJ` by the subtraction of the
    value of the force at :math:`r_{\mathrm{cut}}`, such that the force
    smoothly goes to zero at the cut-off. The potential is modified by a linear
    function. This potential can be used as a substitute for `LJ`,
    when the exact analytical form of the latter is not required but a smaller
    cut-off radius is desired for computational efficiency. See `Toxvaerd et.
    al. 2011 <http://dx.doi.org/10.1063/1.3558787>`_ for a discussion of this
    potential.

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V(r) = & 4 \varepsilon \left[ \left( \frac{\sigma}{r}
          \right)^{12} - \left( \frac{\sigma}{r} \right)^{6}
          \right] + \Delta V(r) & r < r_{\mathrm{cut}}\\
             = & 0 & r \ge r_{\mathrm{cut}} \\
        \end{eqnarray*}

    .. math::

        \Delta V(r) = -(r - r_{\mathrm{cut}}) \frac{\partial
          V_{\mathrm{LJ}}}{\partial r}(r_{\mathrm{cut}})

    See `Pair` for details on how forces are calculated and the
    available energy shifting and smoothing modes. Use `params` dictionary to
    set potential coefficients. The coefficients must be set per unique pair of
    particle types.

    .. py:attribute:: params

        The potential parameters. The dictionary has the following keys:

        * ``epsilon`` (`float`, **required**) - :math:`\varepsilon`
          (in energy units)
        * ``sigma`` (`float`, **required**) - :math:`\sigma`
          (in distance units)

        Type: `TypeParameter` [`tuple` [``particle_type``, ``particle_type``],
        `dict`]

    Example::

        nl = nlist.Cell()
        fslj = pair.ForceShiftedLJ(nlist=nl, r_cut=1.5)
        fslj.params[('A', 'A')] = dict(epsilon=1.0, sigma=1.0)
    """
    _cpp_class_name = "PotentialPairForceShiftedLJ"

    def __init__(self, nlist, r_cut=None, r_on=0., mode='none'):
        # initialize the base class
        super().__init__(nlist, r_cut, r_on, mode)

        params = TypeParameter(
            'params', 'particle_types',
            TypeParameterDict(sigma=float, epsilon=float, len_keys=2))
        self._add_typeparam(params)


class Moliere(Pair):
    r"""Moliere pair potential.

    Args:
        nlist (`hoomd.md.nlist.NList`): Neighbor list
        r_cut (float): Default cutoff radius (in distance units).
        r_on (float): Default turn-on radius (in distance units).
        mode (str): energy shifting/smoothing mode.

    `Moliere` specifies that a Moliere type pair potential should be applied
    between every non-excluded particle pair in the simulation.

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{Moliere}}(r)
          = & \frac{Z_i Z_j e^2}{4 \pi \epsilon_0 r_{ij}} \left[ 0.35 \exp
          \left( -0.3 \frac{r_{ij}}{a_F} \right) + \\
          0.55 \exp \left( -1.2 \frac{r_{ij}}{a_F} \right) + 0.10 \exp
          \left( -6.0 \frac{r_{ij}}{a_F} \right) \right]
          & r < r_{\mathrm{cut}} \\
          = & 0 & r > r_{\mathrm{cut}} \\
        \end{eqnarray*}

    Where each parameter is defined as:

    - :math:`Z_i` - *Z_i* - Atomic number of species i (unitless)
    - :math:`Z_j` - *Z_j* - Atomic number of species j (unitless)
    - :math:`e` - *elementary_charge* - The elementary charge (in charge units)
    - :math:`a_F = \frac{0.8853 a_0}{\left( \sqrt{Z_i} + \sqrt{Z_j}
      \right)^{2/3}}`, where :math:`a_0` is the Bohr radius (in distance units)

    See `Pair` for details on how forces are calculated and the available
    energy shifting and smoothing modes. Use `params` dictionary to set
    potential coefficients. The coefficients must be set per unique pair of
    particle types.

    .. py:attribute:: params

        The potential parameters. The dictionary has the following keys:

        * ``qi`` (`float`, **required**) -
          :math:`q_i = Z_i \frac{e}{\sqrt{4 \pi \epsilon_0}}`
          (in charge units)
        * ``qj`` (`float`, **required**) -
          :math:`q_j = Z_j \frac{e}{\sqrt{4 \pi \epsilon_0}}`
          (in charge units)
        * ``aF`` (`float`, **required**) -
          :math:`a_F = \frac{0.8853 a_0}{\left( \sqrt{Z_i} + \sqrt{Z_j}
          \right)^{2/3}}`

        Type: `TypeParameter` [`tuple` [``particle_type``, ``particle_type``],
        `dict`]

    Example::

        nl = nlist.Cell()
        moliere = pair.Moliere(r_cut = 3.0, nlist=nl)

        Zi = 54
        Zj = 7
        e = 1
        a0 = 1
        aF = 0.8853 * a0 / (np.sqrt(Zi) + np.sqrt(Zj))**(2/3)

        moliere.params[('A', 'B')] = dict(qi=Zi*e, qj=Zj*e, aF=aF)
    """
    _cpp_class_name = "PotentialPairMoliere"

    def __init__(self, nlist, r_cut=None, r_on=0., mode='none'):
        super().__init__(nlist, r_cut, r_on, mode)
        params = TypeParameter(
            'params', 'particle_types',
            TypeParameterDict(qi=float, qj=float, aF=float, len_keys=2))
        self._add_typeparam(params)


class ZBL(Pair):
    r"""ZBL pair potential.

    Args:
        nlist (`hoomd.md.nlist.NList`): Neighbor list
        r_cut (float): Default cutoff radius (in distance units).
        r_on (float): Default turn-on radius (in distance units).
        mode (str): energy shifting/smoothing mode.

    `ZBL` specifies that a Ziegler-Biersack-Littmark pair potential
    should be applied between every non-excluded particle pair in the
    simulation.

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{ZBL}}(r) =
          & \frac{Z_i Z_j e^2}{4 \pi \epsilon_0 r_{ij}} \left[ 0.1818
          \exp \left( -3.2 \frac{r_{ij}}{a_F} \right) \\
          + 0.5099 \exp \left( -0.9423 \frac{r_{ij}}{a_F} \right) \\
          + 0.2802 \exp \left( -0.4029 \frac{r_{ij}}{a_F} \right) \\
          + 0.02817 \exp \left( -0.2016 \frac{r_{ij}}{a_F} \right) \right],
          & r < r_{\mathrm{cut}} \\
          = & 0, & r > r_{\mathrm{cut}} \\
        \end{eqnarray*}

    Where each parameter is defined as:

    - :math:`Z_i` - *Z_i* - Atomic number of species i (unitless)
    - :math:`Z_j` - *Z_j* - Atomic number of species j (unitless)
    - :math:`e` - *elementary_charge* - The elementary charge (in charge units)
    - :math:`a_F = \frac{0.8853 a_0}{ Z_i^{0.23} + Z_j^{0.23} }`, where
      :math:`a_0` is the Bohr radius (in distance units)

    See `Pair` for details on how forces are calculated and the available
    energy shifting and smoothing modes. Use `params` dictionary to set
    potential coefficients.

    .. py:attribute:: params

        The ZBL potential parameters. The dictionary has the following keys:

        * ``q_i`` (`float`, **required**) - :math:`q_i=Z_i \frac{e}{\sqrt{4
          \pi \epsilon_0}}` (in charge units)
        * ``q_j`` (`float`, **required**) - :math:`q_j=Z_j \frac{e}{\sqrt{4
          \pi \epsilon_0}}` (in charge units)
        * ``a_F`` (`float`, **required**) -
          :math:`a_F = \frac{0.8853 a_0}{ Z_i^{0.23} + Z_j^{0.23} }`

        Type: `TypeParameter` [`tuple` [``particle_type``, ``particle_type``],
        `dict`]

    Example::

        nl = nlist.Cell()
        zbl = pair.ZBL(r_cut = 3.0, nlist=nl)

        Zi = 54
        Zj = 7
        e = 1
        a0 = 1
        aF = 0.8853 * a0 / (Zi**(0.23) + Zj**(0.23))

        zbl.params[('A', 'B')] = dict(qi=Zi*e, qj=Zj*e, aF=aF)
    """
    _cpp_class_name = "PotentialPairZBL"

    def __init__(self, nlist, r_cut=None, r_on=0., mode='none'):

        super().__init__(nlist, r_cut, r_on, mode)
        params = TypeParameter(
            'params', 'particle_types',
            TypeParameterDict(qi=float, qj=float, aF=float, len_keys=2))
        self._add_typeparam(params)


class Mie(Pair):
    r"""Mie pair potential.

    Args:
        nlist (`hoomd.md.nlist.NList`): Neighbor list
        r_cut (float): Default cutoff radius (in distance units).
        r_on (float): Default turn-on radius (in distance units).
        mode (str): energy shifting/smoothing mode.

    `Mie` specifies that a Mie pair potential should be applied between every
    non-excluded particle pair in the simulation.

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{mie}}(r)
          = & \left( \frac{n}{n-m} \right) {\left( \frac{n}{m}
          \right)}^{\frac{m}{n-m}} \varepsilon \left[ \left(
          \frac{\sigma}{r} \right)^{n} - \left( \frac{\sigma}{r}
          \right)^{m} \right] & r < r_{\mathrm{cut}} \\
          = & 0 & r \ge r_{\mathrm{cut}} \\
        \end{eqnarray*}

    `Pair` for details on how forces are calculated and the available energy
    shifting and smoothing modes. Use the `params` dictionary to set potential
    coefficients. The coefficients must be set per unique pair of particle
    types.

    .. py:attribute:: params

        The potential parameters. The dictionary has the following keys:

        * ``epsilon`` (`float`, **required**) - :math:`\varepsilon` (in units
          of energy)
        * ``sigma`` (`float`, **required**) - :math:`\sigma` (in distance
          units)
        * ``n`` (`float`, **required**) - :math:`n` (unitless)
        * ``m`` (`float`, **required**) - :math:`m` (unitless)

        Type: `TypeParameter` [`tuple` [``particle_type``, ``particle_type``],
        `dict`]

    Example::

        nl = nlist.Cell()
        mie = pair.Mie(nlist=nl, r_cut=3.0)
        mie.params[('A', 'A')] = dict(epsilon=1.0, sigma=1.0, n=12, m=6)
        mie.r_cut[('A', 'A')] = 2**(1.0/6.0)
        mie.r_on[('A', 'A')] = 2.0
        mie.params[(['A', 'B'], ['C', 'D'])] = dict(epsilon=1.5, sigma=2.0)
    """
    _cpp_class_name = "PotentialPairMie"

    def __init__(self, nlist, r_cut=None, r_on=0., mode='none'):

        super().__init__(nlist, r_cut, r_on, mode)
        params = TypeParameter(
            'params', 'particle_types',
            TypeParameterDict(epsilon=float,
                              sigma=float,
                              n=float,
                              m=float,
                              len_keys=2))

        self._add_typeparam(params)


class ReactionField(Pair):
    r"""Onsager reaction field pair potential.

    Args:
        nlist (`hoomd.md.nlist.NList`): Neighbor list
        r_cut (float): Default cutoff radius (in distance units).
        r_on (float): Default turn-on radius (in distance units).
        mode (str): energy shifting/smoothing mode

    `ReactionField` specifies that an Onsager reaction field pair potential
    should be applied between every non-excluded particle pair in the
    simulation.

    Reaction field electrostatics is an approximation to the screened
    electrostatic interaction, which assumes that the medium can be treated as
    an electrostatic continuum of dielectric constant :math:`\epsilon_{RF}`
    outside the cutoff sphere of radius :math:`r_{\mathrm{cut}}`. See: `Barker
    et. al. 1973 <http://dx.doi.org/10.1080/00268977300102101>`_.

    .. math::

       V_{\mathrm{RF}}(r) = \varepsilon \left[ \frac{1}{r} +
           \frac{(\epsilon_{RF}-1) r^2}{(2 \epsilon_{RF} + 1) r_c^3} \right]

    By default, the reaction field potential does not require charge or diameter
    to be set. Two parameters, :math:`\varepsilon` and :math:`\epsilon_{RF}`
    are needed. If :math:`\epsilon_{RF}` is specified as zero, it will
    represent infinity.

    If *use_charge* is set to True, the following formula is evaluated instead:

    .. math::

        V_{\mathrm{RF}}(r) = q_i q_j \varepsilon \left[ \frac{1}{r} +
          \frac{(\epsilon_{RF}-1) r^2}{(2 \epsilon_{RF} + 1) r_c^3} \right]

    where :math:`q_i` and :math:`q_j` are the charges of the particle pair.

    See `Pair` for details on how forces are calculated and the available
    energy shifting and smoothing modes.  Use the `params` dictionary to set
    potential coefficients. The coefficients must be set per unique pair of
    particle types.

    .. py:attribute:: params

        The potential parameters. The dictionary has the following keys:

        * ``epsilon`` (`float`, **required**) - :math:`\varepsilon` (in units
          of energy*distance)
        * ``eps_rf`` (`float`, **required**) - :math:`\epsilon_{RF}`
          (dimensionless)
        * ``use_charge`` (`boolean`, **optional**) - evaluate pair potential
          using particle charges (*default*: False)

        Type: `TypeParameter` [`tuple` [``particle_type``, ``particle_type``],
        `dict`]

    Example::

        nl = nlist.Cell()
        reaction_field = pair.reaction_field(nl, r_cut=3.0)
        reaction_field.params[('A', 'B')] = dict(epsilon=1.0, eps_rf=1.0)
        reaction_field.params[('B', 'B')] = dict(
            epsilon=1.0, eps_rf=0.0, use_charge=True)
    """
    _cpp_class_name = "PotentialPairReactionField"

    def __init__(self, nlist, r_cut=None, r_on=0., mode='none'):
        super().__init__(nlist, r_cut, r_on, mode)
        params = TypeParameter(
            'params', 'particle_types',
            TypeParameterDict(epsilon=float,
                              eps_rf=float,
                              use_charge=False,
                              len_keys=2))

        self._add_typeparam(params)


class DLVO(Pair):
    r"""DLVO colloidal interaction.

    Args:
        r_cut (float): Default cutoff radius (in distance units).
        nlist (`hoomd.md.nlist.NList`): Neighbor list
        name (str): Name of the force instance.
        d_max (float): Maximum diameter particles in the simulation will have
          (in distance units)

    `DLVO` specifies that a DLVO dispersion and electrostatic interaction should
    be applied between every non-excluded particle pair in the simulation.

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{DLVO}}(r)  = & - \frac{A}{6} \left[
            \frac{2a_1a_2}{r^2 - (a_1+a_2)^2} +
            \frac{2a_1a_2}{r^2 - (a_1-a_2)^2} \\
            + \log \left(
            \frac{r^2 - (a_1+a_2)^2}{r^2 - (a_1-a_2)^2} \right) \right]
            & \\
            & + \frac{a_1 a_2}{a_1+a_2} Z e^{-\kappa(r - (a_1+a_2))}
            & r < (r_{\mathrm{cut}} + \Delta) \\
            = & 0 & r \ge (r_{\mathrm{cut}} + \Delta)
        \end{eqnarray*}

    where :math:`a_i` is the radius of particle :math:`i`, :math:`\Delta = (d_i
    + d_j)/2` and :math:`d_i` is the diameter of particle :math:`i`.

    The first term corresponds to the attractive van der Waals interaction with
    :math:`A` being the Hamaker constant, the second term to the repulsive
    double-layer interaction between two spherical surfaces with Z proportional
    to the surface electric potential. See Israelachvili 2011, pp. 317.

    The DLVO potential does not need charge, but does need diameter. See
    `SLJ` for an explanation on how diameters are handled in the
    neighbor lists.

    Due to the way that DLVO modifies the cutoff condition, it will not function
    properly with the xplor shifting mode. See `Pair` for details on
    how forces are calculated and the available energy shifting and smoothing
    modes.

    .. py:attribute:: params

        The potential parameters. The dictionary has the following keys:

        * ``epsilon`` (`float`, **required**) - :math:`\varepsilon` (in units
          of energy)
        * ``kappa`` (`float`, **required**) - scaling parameter
          :math:`\kappa` (in units of 1/distance)
        * ``Z`` (`float`, **required**) - :math:`Z` (in units of 1/distance)
        * ``A`` (`float`, **required**) - :math:`A` (in units of energy)

        Type: `TypeParameter` [`tuple` [``particle_type``, ``particle_type``],
        `dict`]

    Example::

        nl = nlist.cell()
        dlvo = hoomd.md.pair.DLVO(nlist=nl)
        dlvo.params[('A', 'A')] = {"epsilon": 1.0, "kappa": 1.0}
        dlvo.params[('A', 'B')] = {
            "epsilon": 2.0, "kappa": 0.5, "r_cut": 3.0, "r_on": 2.0}
        dlvo.params[(['A', 'B'], ['C', 'D'])] = {"epsilon": 0.5, "kappa": 3.0}
    """
    _cpp_class_name = "PotentialPairDLVO"

    def __init__(self, nlist, r_cut=None, r_on=0., mode='none'):
        if mode == 'xplor':
            raise ValueError("xplor is not a valid mode for the DLVO potential")

        super().__init__(nlist, r_cut, r_on, mode)
        params = TypeParameter(
            'params', 'particle_types',
            TypeParameterDict(kappa=float, Z=float, A=float, len_keys=2))
        self._add_typeparam(params)

        # mode not allowed to be xplor, so re-do param dict entry without that
        # option
        param_dict = ParameterDict(mode=OnlyFrom(['none', 'shift']))
        self._param_dict.update(param_dict)
        self.mode = mode

        # this potential needs diameter shifting on
        self._nlist.diameter_shift = True


class Buckingham(Pair):
    r"""Buckingham pair potential.

    Args:
        nlist (`hoomd.md.nlist.NList`): Neighbor list
        r_cut (float): Default cutoff radius (in distance units).
        r_on (float): Default turn-on radius (in distance units).
        mode (str): energy shifting/smoothing mode

    `Buckingham` specifies that a Buckingham pair potential should be applied
    between every non-excluded particle pair in the simulation.

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{Buckingham}}(r) = & A \exp\left(-\frac{r}{\rho}\right)
          - \frac{C}{r^6} & r < r_{\mathrm{cut}} \\
          = & 0 & r \ge r_{\mathrm{cut}} \\
        \end{eqnarray*}

    See `Pair` for details on how forces are calculated and the available
    energy shifting and smoothing modes.  Use the `params` dictionary to set
    potential coefficients.

    .. py:attribute:: params

        The potential parameters. The dictionary has the following keys:

        * ``A`` (`float`, **required**) - :math:`A` (in energy units)
        * ``rho`` (`float`, **required**) - :math:`\rho` (in distance units)
        * ``C`` (`float`, **required**) - :math:`C` (in energy units)

        Type: `TypeParameter` [`tuple` [``particle_type``, ``particle_type``],
        `dict`]

    Example::

        nl = nlist.Cell()
        buck = pair.Buckingham(nl, r_cut=3.0)
        buck.params[('A', 'A')] = {'A': 2.0, 'rho'=0.5, 'C': 1.0}
        buck.params[('A', 'B')] = dict(A=1.0, rho=1.0, C=1.0)
        buck.params[('B', 'B')] = dict(A=2.0, rho=2.0, C=2.0)
    """

    _cpp_class_name = "PotentialPairBuckingham"

    def __init__(self, nlist, r_cut=None, r_on=0., mode='none'):
        super().__init__(nlist, r_cut, r_on, mode)
        params = TypeParameter(
            'params', 'particle_types',
            TypeParameterDict(A=float, rho=float, C=float, len_keys=2))
        self._add_typeparam(params)


class LJ1208(Pair):
    r"""Lennard-Jones 12-8 pair potential.

    Args:
        nlist (`hoomd.md.nlist.NList`): Neighbor list
        r_cut (float): Default cutoff radius (in distance units).
        r_on (float): Default turn-on radius (in distance units).
        mode (str): energy shifting/smoothing mode

    `LJ1208` specifies that a Lennard-Jones 12-8 pair potential should be
    applied between every non-excluded particle pair in the simulation.

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{LJ}}(r)
          = & 4 \varepsilon \left[ \left( \frac{\sigma}{r} \right)^{12} -
          \left( \frac{\sigma}{r} \right)^{8} \right]
          & r < r_{\mathrm{cut}} \\
          = & 0 & r \ge r_{\mathrm{cut}} \\
        \end{eqnarray*}

    See `Pair` for details on how forces are calculated and the available
    energy shifting and smoothing modes.  Use the `params` dictionary to set
    potential coefficients.

    .. py:attribute:: params

        The potential parameters. The dictionary has the following keys:

        * ``epsilon`` (`float`, **required**) - energy parameter
          :math:`\varepsilon` (in energy units)
        * ``sigma`` (`float`, **required**) - particle size :math:`\sigma`
          (in distance units)

        Type: `TypeParameter` [`tuple` [``particle_type``, ``particle_type``],
        `dict`]

    Example::

        nl = nlist.Cell()
        lj1208 = pair.LJ1208(nl, r_cut=3.0)
        lj1208.params[('A', 'A')] = {'sigma': 1.0, 'epsilon': 1.0}
        lj1208.params[('A', 'B')] = dict(epsilon=2.0, sigma=1.0)
    """
    _cpp_class_name = "PotentialPairLJ1208"

    def __init__(self, nlist, r_cut=None, r_on=0., mode='none'):
        super().__init__(nlist, r_cut, r_on, mode)
        params = TypeParameter(
            'params', 'particle_types',
            TypeParameterDict(epsilon=float, sigma=float, len_keys=2))
        self._add_typeparam(params)


class LJ0804(Pair):
    r"""Lennard-Jones 8-4 pair potential.

    Args:
        nlist (`hoomd.md.nlist.NList`): Neighbor list
        r_cut (float): Default cutoff radius (in distance units).
        r_on (float): Default turn-on radius (in distance units).
        mode (str): energy shifting/smoothing mode

    `LJ0804` specifies that a Lennard-Jones 8-4 pair potential should be
    applied between every non-excluded particle pair in the simulation.

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{LJ}}(r)
          = & 4 \varepsilon \left[ \left( \frac{\sigma}{r} \right)^{8} -
          \left( \frac{\sigma}{r} \right)^{4} \right]
          & r < r_{\mathrm{cut}} \\
          = & 0 & r \ge r_{\mathrm{cut}} \\
        \end{eqnarray*}

    See `Pair` for details on how forces are calculated and the
    available energy shifting and smoothing modes.  Use the `params` dictionary
    to set potential coefficients. The coefficients must be set per
    unique pair of particle types.

    .. py:attribute:: params

        The LJ potential parameters. The dictionary has the following keys:

        * ``epsilon`` (`float`, **required**) -
          energy parameter :math:`\varepsilon` (in energy units)
        * ``sigma`` (`float`, **required**) -
          particle size :math:`\sigma` (in distance units)

        Type: `TypeParameter` [`tuple` [``particle_type``, ``particle_type``],
        `dict`]

    Example::

        nl = nlist.Cell()
        lj0804 = pair.LJ0804(nl, r_cut=3.0)
        lj0804.params[('A', 'A')] = {'sigma': 1.0, 'epsilon': 1.0}
        lj0804.params[('A', 'B')] = dict(epsilon=2.0, sigma=1.0)
        lj0804.r_cut[('A', 'B')] = 3.0
    """
    _cpp_class_name = "PotentialPairLJ0804"

    def __init__(self, nlist, r_cut=None, r_on=0., mode='none'):
        super().__init__(nlist, r_cut, r_on, mode)
        params = TypeParameter(
            'params', 'particle_types',
            TypeParameterDict(epsilon=float, sigma=float, len_keys=2))
        self._add_typeparam(params)


class Fourier(Pair):
    r"""Fourier pair potential.

    Args:
        nlist (`hoomd.md.nlist.NList`): Neighbor list
        r_cut (float): Default cutoff radius (in distance units).
        r_on (float): Default turn-on radius (in distance units).
        mode (str): Energy shifting mode.

    `Fourier` specifies that a Fourier pair potential should be applied between
    every non-excluded particle pair in the simulation.

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{Fourier}}(r)
          = & \frac{1}{r^{12}} + \frac{1}{r^2}\sum_{n=1}^4
          [a_n cos(\frac{n \pi r}{r_{cut}}) +
          b_n sin(\frac{n \pi r}{r_{cut}})]
          & r < r_{\mathrm{cut}}  \\
          = & 0 & r \ge r_{\mathrm{cut}} \\
        \end{eqnarray*}

        where:
        \begin{eqnarray*}
        a_1 = \sum_{n=2}^4 (-1)^n a_n
        \end{eqnarray*}

        \begin{eqnarray*}
        b_1 = \sum_{n=2}^4 n (-1)^n b_n
        \end{eqnarray*}

        is calculated to enforce close to zero value at r_cut.

    See `Pair` for details on how forces are calculated and the available
    energy shifting and smoothing modes. Use `params` dictionary to set
    potential coefficients. The coefficients must be set per unique pair of
    particle types.

    .. py:attribute:: params

        The Fourier potential parameters. The dictionary has the following
        keys:

        * ``a`` (`float`, **required**) - array of 3 values corresponding to
          a2, a3 and a4 in the Fourier series, unitless)
        * ``b`` (`float`, **required**) - array of 3 values corresponding to
          b2, b3 and b4 in the Fourier series, unitless)

        Type: `TypeParameter` [`tuple` [``particle_type``, ``particle_type``],
        `dict`]

    Example::

        nl = nlist.Cell()
        fourier = pair.Fourier(r_cut=3.0, nlist=nl)
        fourier.params[('A', 'A')] = dict(a=[a2,a3,a4], b=[b2,b3,b4])
    """
    _cpp_class_name = "PotentialPairFourier"

    def __init__(self, nlist, r_cut=None, r_on=0., mode='none'):
        super().__init__(nlist, r_cut, r_on, mode)
        params = TypeParameter(
            'params', 'particle_types',
            TypeParameterDict(a=(float, float, float),
                              b=(float, float, float),
                              len_keys=2))
        self._add_typeparam(params)


class OPP(Pair):
    r"""Oscillating pair potential.

    Args:
        nlist (:py:mod:`hoomd.md.nlist.NList`): Neighbor list
        r_cut (float): Default cutoff radius (in distance units).
        r_on (float): Default turn-on radius (in distance units).
        mode (str): energy shifting/smoothing mode

    `OPP` specifies that an oscillating pair potential should be applied between
    every non-excluded particle pair in the simulation. The OPP potential can
    be used to model metallic interactions.

    .. math::
        :nowrap:

        \begin{equation*}
        V_{\mathrm{OPP}}(r) = C_1 r^{-\eta_1}
            + C_2 r^{-\eta_2} \cos{\left(k r - \phi\right)}
        \end{equation*}

    See `Pair` for details on how forces are calculated and the available energy
    shifting and smoothing modes.  Use `params` dictionary to set potential
    coefficients. The coefficients must be set per unique pair of particle
    types.

    The potential comes from Marek Mihalkovič and C. L. Henley 2012
    `paper link`_.

    .. _paper link: https://dx.doi.org/10.1103/PhysRevB.85.092102

    .. py:attribute:: params

        The OPP potential parameters. The dictionary has the following keys:

        * ``C1`` (`float`, **required**) -
          Energy scale of the first term :math:`C_1` (energy units)
        * ``C2`` (`float`, **required**) -
          Energy scale of the second term :math:`C_2` (energy units)
        * ``eta1`` (`float`, **required**) -
          The inverse power to take :math:`r` to in the first term,
          :math:`\eta_1` (unitless).
        * ``eta2`` (`float`, **required**) -
          The inverse power to take :math:`r` to in the second term
          :math:`\eta_2` (unitless).
        * ``k`` (`float`, **required**) -
          oscillation frequency :math:`k` (inverse distance units)
        * ``phi`` (`float`, **required**) -
          potential phase shift :math:`\phi` (unitless)

        Type: `TypeParameter` [`tuple` [``particle_type``, ``particle_type``],
        `dict`]

    Example::

        nl = nlist.Cell()
        opp = pair.OPP(nl, r_cut=3.0)
        opp.params[('A', 'A')] = {
            'C1': 1., 'C2': 1., 'eta1': 15,
            'eta2': 3, 'k': 1.0, 'phi': 3.14}
        opp.r_cut[('A', 'B')] = 3.0
    """
    _cpp_class_name = "PotentialPairOPP"

    def __init__(self, nlist, r_cut=None, r_on=0., mode='none'):
        super().__init__(nlist, r_cut, r_on, mode)
        params = TypeParameter(
            'params', 'particle_types',
            TypeParameterDict(C1=float,
                              C2=float,
                              eta1=float,
                              eta2=float,
                              k=float,
                              phi=float,
                              len_keys=2))
        self._add_typeparam(params)


class TWF(Pair):
    r"""Pair potential model for globular proteins.

    This potential was introduced by Ten-wolde and Daan Frenkel in 1997 for
    studying globular protein crystallization. The potential has the following
    form:

    .. math::
        :nowrap:

        \begin{equation}
        V_{\mathrm{TWF}}(r) = \frac{4 \epsilon}{\alpha^2} {\left[
        {\left(\frac{\sigma^2}{r^2} - 1 \right)}^6 -
        \alpha {\left(\frac{\sigma^2}{r^2} - 1 \right)}^3\right]}
        \end{equation}

    See `hoomd.md.pair.Pair` for details on how forces are calculated and the
    available energy shifting and smoothing modes.  Use `params` dictionary
    to set potential coefficients. The coefficients must be set per
    unique pair of particle types.

    The potential comes from Pieter Rein ten Wolde and Daan Frenkel 1997
    `paper link`_.

    .. _paper link: https://dx.doi.org/10.1126/science.277.5334.1975

    .. py:attribute:: params

        The LJ potential parameters. The dictionary has the following keys:

        * ``epsilon`` (`float`, **required**) -
          energy parameter :math:`\varepsilon` (units: [energy])
        * ``sigma`` (`float`, **required**) -
          particle size :math:`\sigma` (units: [length])
        * ``alpha`` (`float`, **required**) -
          controls well-width :math:`\alpha` (unitless)

        Type: `TypeParameter` [`tuple` [``particle_type``, ``particle_type``],
        `dict`]

    Example::

        nl = nlist.Cell()
        twf = hoomd.md.pair.TWF(nl, r_cut=3.0)
        twf.params[('A', 'A')] = {'sigma': 1.0, 'epsilon': 1.0, 'alpha': 50.0}
        twf.r_cut[('A', 'B')] = 3.0
    """
    _cpp_class_name = "PotentialPairTWF"

    def __init__(self, nlist, r_cut=None, r_on=0.0, mode='none'):
        super().__init__(nlist, r_cut, r_on, mode)
        params = TypeParameter(
            'params', 'particle_types',
            TypeParameterDict(epsilon=float,
                              sigma=float,
                              alpha=float,
                              len_keys=2))
        self._add_typeparam(params)
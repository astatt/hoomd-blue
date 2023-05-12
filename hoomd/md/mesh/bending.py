# Copyright (c) 2009-2022 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

r"""Mesh Bending potentials.

Mesh bending force classes apply a force and virial to every mesh vertex
particle based on the local curvature :math:`K` of the given mesh triangulation.

.. math::

    U_\mathrm{mesh bending} = \sum_{j \in \mathrm{mesh}}
    U_{j}(K(\mathbf{r}_j))

The curvature at each vertex particle :math:`j` is determined at its position
:math:`\mathbf{r}_j`.

See Also:
   See the documentation in `hoomd.mesh.Mesh` for more information on the
   initialization of the mesh object.


Note:
   The mesh bond forces are computed over the mesh data structure and not the
   separate bond data structure. Hence, the mesh bonds are defined exclusively
   by the mesh triangulation as HOOMD-blue automatically constructs the mesh
   bond pairs based on ``triangulation`` in the `hoomd.mesh.Mesh` object.
   The bonds should **not** be defined separately in the `hoomd.State` member
   ``bond_group``!

"""

from hoomd.md.mesh.potential import MeshPotential
from hoomd.data.typeparam import TypeParameter
from hoomd.data.parameterdicts import TypeParameterDict


class BendingRigidity(MeshPotential):
    r"""Bending potential.

    :py:class:`BendingRigidity` specifies a bending energy applied to
    all particles within the mesh.

    .. math::

        U(i) = \frac{1}{2} k \sum_{j \in \mathrm{Neigh}(i)}
        ( 1 - cos(\theta_{ij}))

    with :math:`\theta_{ij}` is the angle between the two normal
    directors of the bordering triangles of bong :math:`i` and :math:`j`.

    Args:
        mesh (:py:mod:`hoomd.mesh.Mesh`): Mesh data structure constraint.

    Attributes:
        parameter (TypeParameter[dict]):
            The parameter of the bending energy for the defined mesh.
            As the mesh can only have one type a type name does not have
            to be stated. The dictionary has the following keys:

            * ``k`` (`float`, **required**) - bending stiffness
              :math:`[\mathrm{energy}]`

    Examples::

        bending_potential = mesh.bond.BendingRigidity(mesh)
        bending_potential.parameter = dict(k=10.0)
    """
    _cpp_class_name = "BendingRigidityMeshForceCompute"

    def __init__(self, mesh):
        params = TypeParameter("params", "types",
                               TypeParameterDict(k=float, len_keys=1))
        self._add_typeparam(params)

        super().__init__(mesh)

class Helfrich(MeshPotential):
    r"""Helfrich bending potential.

    :py:class:`Helfrich` specifies a Helfrich bending energy applied to
    all particles within the mesh.

    .. math::

        U(i) = \frac{1}{2} k \frac{1}{\sigma_i}\left( \sum_{j \in
        \mathrm{Neigh}(i)} \frac{\sigma_{ij}}{l_{ij}} (\mathbf{r}_j
        - \mathbf{r}_k) \right)^2

    with the area of the dual cell of vertex i
    :math:`\sigma_i=(\sum_{j \in \mathrm{Neigh}(i)}\sigma_{ij})/4`, the
    length of the bond in the dual lattice  :math:`\sigma_{ij}=
    r_{ij}(\text{cot}\theta_1+\text{cot}\theta_2)/2` and the angles
    :math:`\theta_1` and :math:`\theta_2` opposite to the shared bond of
    vertex :math:`i` and :math:`j`.

    See Also:
    * `Gompper and Kroll 1996 <https://doi.org/10.1051/jp1:1996246>`__
    * `Helfrich 1973 <https://doi.org/10.1515/znc-1973-11-1209>`__

    Args:
        mesh (:py:mod:`hoomd.mesh.Mesh`): Mesh data structure constraint.

    Attributes:
        parameter (TypeParameter[dict]):
            The parameter of the Helfrich energy for the defined mesh.
            As the mesh can only have one type a type name does not have
            to be stated. The dictionary has the following keys:

            * ``k`` (`float`, **required**) - bending stiffness
              :math:`[\mathrm{energy}]`

    Examples::

        helfrich_potential = mesh.bond.Helfrich(mesh)
        helfrich_potential.params["mesh"] = dict(k=10.0)
    """
    _cpp_class_name = "HelfrichMeshForceCompute"

    def __init__(self, mesh):
        params = TypeParameter("params", "types",
                               TypeParameterDict(k=float, len_keys=1))
        self._add_typeparam(params)

        super().__init__(mesh)

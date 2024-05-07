# Copyright (c) 2009-2023 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

# Import the plugin module.
from hoomd import pair_plugin

# Import the hoomd Python package.
import hoomd

import itertools
import pytest
import numpy as np

#TODO: change this function to be squarewell, add parameters needed
def gcms_adj_force_and_energy(dx, w, sigma, a, q, r_cut, shift=False):

    dr = np.linalg.norm(dx)

    print(dr)

    if dr >= r_cut:
        return np.array([0.0, 0.0, 0.0], dtype=np.float64), 0.0

    exponent = np.exp(q*(dr - sig - w)/w)

    e = a*(-1*(w/(dr - sigma + w))**q + 1/(1 + exponent))
    f = a*(q*(w/(dr - sigma + w))**(q + 1) - q*exponent/(w + w*exponent)**2.0)
    if shift:
        e = e
    print(f, e)
    return f, e

#TODO: decide on distances and parameters to use for test

# Build up list of parameters.
distances = np.array([1.05, 1.6])
q = 16
ws = [0.1]
as = [10]
sigmas = [1.0]
# No need to test "xplor", as that is handled outside of the plugin impl.
modes = ["none"]

energies = [[9.981422110291863, 9.981422110291863], [-3.009063548896661e-13, -3.009063548896661e-13]]
forces = [[0.3739867969066353, -0.3739867969066353], [-6.877859540335224e-13, 6.877859540335224e-13]]

testdata = list(itertools.product(distances, ws, sigmas, as, q, modes, energies, forces))



@pytest.mark.parametrize("distance, w, sigma, a, q, mode, energies, forces", testdata)
def test_force_and_energy_eval(simulation_factory,
                               two_particle_snapshot_factory, distance, w, sigma,
                               a, q, mode, energies, forces):

    # Build the simulation from the factory fixtures defined in
    # hoomd/conftest.py.
    sim = simulation_factory(two_particle_snapshot_factory(d=distance))

    # Setup integrator and force.
    integrator = hoomd.md.Integrator(dt=0.001)
    nve = hoomd.md.methods.ConstantVolume(hoomd.filter.All())

    cell = hoomd.md.nlist.Cell(buffer=0.4)
    #TODO: change names to correct square well ones
    example_pair: hoomd.md.pair.Pair = pair_plugin.pair.GCMSAdjPair(
        cell, default_r_cut=3, mode=mode)
    example_pair.params[("A", "A")] = dict(w = w, sigma = sigma, a = a, q = q)
    integrator.forces = [example_pair]
    integrator.methods = [nve]

    sim.operations.integrator = integrator

    sim.run(0)
    snap = sim.state.get_snapshot()
    if snap.communicator.rank == 0:
        vec_dist = snap.particles.position[1] - snap.particles.position[0]

        shift = mode == "none"
        # Compute force and energy from Python
        f, e = gcms_adj_force_and_energy(vec_dist, w, sigma, a, q, 3, shift)
        # e /= 2.0

    # Test that the forces and energies match that predicted by the Python
    # implementation.
    # forces_hand = forces
    if snap.communicator.rank == 0:
        np.testing.assert_array_almost_equal(forces, [-f, f], decimal=6)

    # en_hand = energies

    # print(en_hand, [e, e])

    if snap.communicator.rank == 0:
        np.testing.assert_array_almost_equal(energies, [e, e], decimal=6)

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
def expanded_yukawa_force_and_energy(dx, epsilon, kappa, delta, r_cut, shift=False):

    dr = np.linalg.norm(dx)

    print(dr)

    if dr >= r_cut:
        return np.array([0.0, 0.0, 0.0], dtype=np.float64), 0.0

    rinv = 1/dr;
    delta_dist = dr - delta
    rinv_delt = 1/delta_dist
    kappa_delt = kappa*delta_dist
    exponent = np.exp(-kappa_delt)

    f = epsilon*exponent*(kappa_delt + 1)*rinv_delt*rinv_delt
    e = epsilon*exponent*rinv_delt
    if shift:
        e = e
    print(f, e)
    return f, e

#TODO: decide on distances and parameters to use for test 

# Build up list of parameters.
distances = np.array([1.01, 1.05, 2.0])
eps = [1.0, 2.0]
kaps = [1, 10]
delts = [0.5]

# No need to test "xplor", as that is handled outside of the plugin impl.
modes = ["none"]

energies = [[1.177442311, 1.177442311], [1.048999655, 1.048999655], [0.1487534401, 0.1487534401], [2.354884623, 2.354884623], [2.09799931, 2.09799931], [0.2975068802, 0.2975068802], [0.02390881006, 0.02390881006], [0.01486098705, 0.01486098705], [4.07*10**(-7.0), 4.07*10**(-7.0)], [0.01195440503, 0.01195440503], [0.007430493524, 0.007430493524], [2.03*10**(-7.0), 2.03*10**(-7.0)]]
forces = [[-3.486152726, 3.486152726], [-2.956271756, 2.956271756], [-0.2479224002, 0.2479224002], [-6.972305452, 6.972305452], [-5.912543511, 5.912543511], [-0.4958448003, 0.4958448003], [-0.2859681203, 0.2859681203], [-0.1756298469, 0.1756298469], [-4.35061078*10**(-6), 4.35061078*10**(-6)], [-0.1429840602, 0.1429840602], [-0.08781492347, 0.08781492347], [-2.17530539*10**(-6), 2.17530539*10**(-6)]]

testdata = [(distances[0], eps[0], kaps[0], delts[0], modes[0], energies[0], forces[0]), (distances[1], eps[0], kaps[0], delts[0], modes[0], energies[1], forces[1]), (distances[2], eps[0], kaps[0], delts[0], modes[0], energies[2], forces[2]), (distances[0], eps[1], kaps[0], delts[0], modes[0], energies[3], forces[3]), (distances[1], eps[1], kaps[0], delts[0], modes[0], energies[4], forces[4]), (distances[2], eps[1], kaps[0], delts[0], modes[0], energies[5], forces[5]), (distances[0], eps[1], kaps[1], delts[0], modes[0], energies[6], forces[6]), (distances[1], eps[1], kaps[1], delts[0], modes[0], energies[7], forces[7]), (distances[2], eps[1], kaps[1], delts[0], modes[0], energies[8], forces[8]), (distances[0], eps[0], kaps[1], delts[0], modes[0], energies[9], forces[9]), (distances[1], eps[0], kaps[1], delts[0], modes[0], energies[10], forces[10]), (distances[2], eps[0], kaps[1], delts[0], modes[0], energies[11], forces[11])]

@pytest.mark.parametrize("distance, epsilon, kappa, delta, mode, energies, forces", testdata)
def test_force_and_energy_eval(simulation_factory,
                               two_particle_snapshot_factory, distance, epsilon, kappa,
                               delta, mode, energies, forces):

    # Build the simulation from the factory fixtures defined in
    # hoomd/conftest.py.
    sim = simulation_factory(two_particle_snapshot_factory(d=distance))

    # Setup integrator and force.
    integrator = hoomd.md.Integrator(dt=0.001)
    nve = hoomd.md.methods.ConstantVolume(hoomd.filter.All())

    cell = hoomd.md.nlist.Cell(buffer=0.4)
    #TODO: change names to correct square well ones 
    example_pair: hoomd.md.pair.Pair = pair_plugin.pair.ExpandedYukawaPair(
        cell, default_r_cut=3, mode=mode)
    example_pair.params[("A", "A")] = dict(epsilon = epsilon, kappa = kappa, delta = delta)
    integrator.forces = [example_pair]
    integrator.methods = [nve]

    sim.operations.integrator = integrator

    sim.run(0)
    snap = sim.state.get_snapshot()
    if snap.communicator.rank == 0:
        vec_dist = snap.particles.position[1] - snap.particles.position[0]

        shift = mode == "none"
        # Compute force and energy from Python
        f, e = expanded_yukawa_force_and_energy(vec_dist, epsilon, kappa, delta, 3, shift)
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

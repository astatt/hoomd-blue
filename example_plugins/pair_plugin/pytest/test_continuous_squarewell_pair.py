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
def squarewell_force_and_energy(dx, n, m, lambda_val, r_cut, shift=False):

    dr = np.linalg.norm(dx)

    if dr >= r_cut:
        return np.array([0.0, 0.0, 0.0], dtype=np.float64), 0.0

    f = n/2.0*(1/dx)**(n + 1) - (m*2*dx - lambda_val - 1)*np.exp(-m*(dx - 1)*(dx - lambda_val))/(1 + np.exp(-m*(dx - 1)*(dx - lambda_val)))**2.0
    e = 1.0/2.0*((1/dx)**n + (1 - np.exp(-m*(dx - 1)*(dx - lambda_val)))/(1 + np.exp(-m*(dx - 1)*(r - lambda_val))))
    if shift:
        e = e

    return f, e

#TODO: decide on distances and parameters to use for test 

# Build up list of parameters.
distances = np.linspace(1.01, 1.05, 2)
ns = [400, 2500]
ms = [400, 20000]
lambda_vals = [1.5, 1.05]
# No need to test "xplor", as that is handled outside of the plugin impl.
modes = ["none"]

testdata = list(itertools.product(distances, ns, ms, lambda_vals, modes))


@pytest.mark.parametrize("distance, n, m, lamnbda_val, mode", testdata)
def test_force_and_energy_eval(simulation_factory,
                               two_particle_snapshot_factory, distance, n, m,
                               lambda_val, mode):

    # Build the simulation from the factory fixtures defined in
    # hoomd/conftest.py.
    sim = simulation_factory(two_particle_snapshot_factory(d=distance))

    # Setup integrator and force.
    integrator = hoomd.md.Integrator(dt=0.001)
    nve = hoomd.md.methods.ConstantVolume(hoomd.filter.All())

    cell = hoomd.md.nlist.Cell(buffer=0.4)
    #TODO: change names to correct square well ones 
    example_pair: hoomd.md.pair.Pair = pair_plugin.pair.ExamplePairContinuousSquareWell(
        cell, default_r_cut=3, mode=mode)
    example_pair.params[("A", "A")] = dict(n=n, m=m, lambda_val = lambda_val)
    integrator.forces = [example_pair]
    integrator.methods = [nve]

    sim.operations.integrator = integrator

    sim.run(0)
    snap = sim.state.get_snapshot()
    if snap.communicator.rank == 0:
        vec_dist = snap.particles.position[1] - snap.particles.position[0]

        # Compute force and energy from Python
        shift = mode == "shift"
        f, e = squarewell_force_and_energy(vec_dist, n, m, lambda_val, 3, shift)
        e /= 2.0

    # Test that the forces and energies match that predicted by the Python
    # implementation.
    forces = example_pair.forces
    if snap.communicator.rank == 0:
        np.testing.assert_array_almost_equal(forces, [-f, f], decimal=6)

    energies = example_pair.energies
    if snap.communicator.rank == 0:
        np.testing.assert_array_almost_equal(energies, [e, e], decimal=6)

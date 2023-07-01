# Copyright (c) 2009-2022 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Example Updater."""

# Import the C++ module.
from hoomd.updater_plugin import _updater_plugin

# Import the hoomd Python package.
import hoomd
from hoomd import operation
from hoomd.operation import Updater
from hoomd.box import Box
from hoomd.data.parameterdicts import ParameterDict
from hoomd.variant import Variant, Constant
from hoomd import _hoomd
from hoomd.filter import ParticleFilter, All
from hoomd.trigger import Periodic


class BoxResizeConstVolume(operation.Updater):
    """Resizes the box keeping the overall volume constant

    .. math::

        \\begin{align*}
        L_{x}' &= \\lambda L_{2x} + (1 - \\lambda) L_{1x} \\\\
        L_{y}' &= \\lambda L_{2y} + (1 - \\lambda) L_{1y} \\\\
        L_{z}' &= \\lambda L_{2z} + (1 - \\lambda) L_{1z} \\\\
        xy' &= \\lambda xy_{2} + (1 - \\lambda) xy_{2} \\\\
        xz' &= \\lambda xz_{2} + (1 - \\lambda) xz_{2} \\\\
        yz' &= \\lambda yz_{2} + (1 - \\lambda) yz_{2} \\\\
        \\end{align*}


    Args:
        trigger (hoomd.trigger.trigger_like): The trigger to activate this
            updater.
        box1 (hoomd.box.box_like): The box associated with the minimum of the
            passed variant.
        box2 (hoomd.box.box_like): The box associated with the maximum of the
            passed variant.
        variant (hoomd.variant.variant_like): A variant used to interpolate
            between the two boxes.
        filter (hoomd.filter.filter_like): The subset of particle positions
            to update.

    Attributes:
        box1 (hoomd.Box): The box associated with the minimum of the
            passed variant.
        box2 (hoomd.Box): The box associated with the maximum of the
            passed variant.
        variant (hoomd.variant.Variant): A variant used to interpolate between
            the two boxes.
        trigger (hoomd.trigger.Trigger): The trigger to activate this updater.
        filter (hoomd.filter.filter_like): The subset of particles to
            update.
    """

    def __init__(self, trigger, L1, L2, direction,variant, filter=All()):
        params = ParameterDict(L1=float,
                               L2=float,
                               direction=int,
                               variant=Variant,
                               filter=ParticleFilter)
        params['L1'] = L1
        params['L2'] = L2
        params['variant'] = variant
        params['direction'] = direction
        params['trigger'] = trigger
        params['filter'] = filter
        self._param_dict.update(params)
        super().__init__(trigger)

    def _attach_hook(self):
        group = self._simulation.state._get_group(self.filter)
        self._cpp_obj = _updater_plugin.BoxResizeConstVolumeUpdater(
            self._simulation.state._cpp_sys_def, self.trigger,
            self.L1, self.L2,self.direction, self.variant, group)

    def get_box(self, timestep):
        """Get the box for a given timestep.

        Args:
            timestep (int): The timestep to use for determining the resized
                box.

        Returns:
            Box: The box used at the given timestep.
            `None` before the first call to `Simulation.run`.
        """
        if self._attached:
            timestep = int(timestep)
            if timestep < 0:
                raise ValueError("Timestep must be a non-negative integer.")
            return Box._from_cpp(self._cpp_obj.get_current_box(timestep))
        else:
            return None

    @staticmethod
    def update(state, box, filter=All()):
        """Immediately scale the particle in the system state to the given box.

        Args:
            state (State): System state to scale.
            box (hoomd.box.box_like): New box.
            filter (hoomd.filter.filter_like): The subset of particles to
                update.
        """
        group = state._get_group(filter)
        updater = _updater_plugin.BoxResizeConstVolumeUpdater(state._cpp_sys_def, Periodic(1),
                                          state.box._cpp_obj, box._cpp_obj,
                                          Constant(1), group)
        updater.update(state._simulation.timestep)

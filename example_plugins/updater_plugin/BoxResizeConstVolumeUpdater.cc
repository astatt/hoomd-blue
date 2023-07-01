// Copyright (c) 2009-2022 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

/*! \file BoxResizeConstVolumeUpdater.cc
    \brief Defines the BoxResizeUpdater class
*/

#include "BoxResizeConstVolumeUpdater.h"

#include <iostream>
#include <math.h>
#include <stdexcept>

using namespace std;

namespace hoomd
    {
/*! \param sysdef System definition containing the particle data to set the box size on
    \param Lx length of the x dimension over time
    \param Ly length of the y dimension over time
    \param Lz length of the z dimension over time

    The default setting is to scale particle positions along with the box.
*/

BoxResizeConstVolumeUpdater::BoxResizeConstVolumeUpdater(std::shared_ptr<SystemDefinition> sysdef,
                                   std::shared_ptr<Trigger> trigger,
                                   Scalar L1,
                                   Scalar L2,
                                   unsigned int direction,
                                   std::shared_ptr<Variant> variant,
                                   std::shared_ptr<ParticleGroup> group)
    : Updater(sysdef, trigger), m_L1(L1), m_L2(L2), m_direction(direction), m_variant(variant), m_group(group)
    {
    assert(m_pdata);
    assert(m_variant);
    m_orig_box = m_pdata->getGlobalBox();
    m_exec_conf->msg->notice(5) << "Constructing BoxResizeConstVolumeUpdater" << endl;
    }

BoxResizeConstVolumeUpdater::~BoxResizeConstVolumeUpdater()
    {
    m_exec_conf->msg->notice(5) << "Destroying BoxResizeConstVolumeUpdater" << endl;
    }

/// Get box1
Scalar BoxResizeConstVolumeUpdater::getL1()
    {
    return m_L1;
    }

/// Set a new box1
void BoxResizeConstVolumeUpdater::setL1(Scalar L1)
    {
    m_L1 = L1;
    }

/// Get box2
Scalar BoxResizeConstVolumeUpdater::getL2()
    {
    return m_L2;
    }

void BoxResizeConstVolumeUpdater::setL2(Scalar L2)
    {
    m_L2 = L2;
    }

/// Get the current box based on the timestep
BoxDim BoxResizeConstVolumeUpdater::getCurrentBox(uint64_t timestep)
    {
    Scalar min = m_variant->min();
    Scalar max = m_variant->max();
    Scalar cur_value = (*m_variant)(timestep);
    Scalar scale = 0;
    if (cur_value == max)
        {
        scale = 1;
        }
    else if (cur_value > min)
        {
        scale = (cur_value - min) / (max - min);
        }
  
    unsigned int direction = m_direction;
      
    Scalar new_L = m_L2 * scale + m_L1 * (1.0 - scale);
    Scalar current_lambda = new_L/m_L1;
   
    Scalar new_scaling_off_directions = 1.0/sqrt(current_lambda);
    Scalar3 b = m_orig_box.getL();
    BoxDim new_box;

    if (direction==0)
    {
    new_box = BoxDim(make_scalar3(
            new_L, 
            new_scaling_off_directions*b.y,
            new_scaling_off_directions*b.z));

    }else if(direction==1)
    {
    new_box = BoxDim(make_scalar3(
            new_scaling_off_directions*b.x,
            new_L,
            new_scaling_off_directions*b.z));

    }
    else
    {
    new_box = BoxDim(make_scalar3(
            new_scaling_off_directions*b.x,
            new_scaling_off_directions*b.y,
            new_L));

    }
    return new_box;
    }

/** Perform the needed calculations to scale the box size
    \param timestep Current time step of the simulation
*/
void BoxResizeConstVolumeUpdater::update(uint64_t timestep)
    {
    Updater::update(timestep);
    m_exec_conf->msg->notice(10) << "Box resize update const volume" << endl;

    // first, compute the new box
    BoxDim new_box = getCurrentBox(timestep);

    // check if the current box size is the same
    BoxDim cur_box = m_pdata->getGlobalBox();

    // only change the box if there is a change in the box dimensions
    if (new_box != cur_box)
        {
        // set the new box
        m_pdata->setGlobalBox(new_box);

        // scale the particle positions (if we have been asked to)
        // move the particles to be inside the new box
        ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(),
                                   access_location::host,
                                   access_mode::readwrite);

        for (unsigned int group_idx = 0; group_idx < m_group->getNumMembers(); group_idx++)
            {
            unsigned int j = m_group->getMemberIndex(group_idx);
            // obtain scaled coordinates in the old global box
            Scalar3 fractional_pos = cur_box.makeFraction(
                make_scalar3(h_pos.data[j].x, h_pos.data[j].y, h_pos.data[j].z));

            // intentionally scale both rigid body and free particles, this
            // may waste a few cycles but it enables the debug inBox checks
            // to be left as is (otherwise, setRV cannot fixup rigid body
            // positions without failing the check)
            Scalar3 scaled_pos = new_box.makeCoordinates(fractional_pos);
            h_pos.data[j].x = scaled_pos.x;
            h_pos.data[j].y = scaled_pos.y;
            h_pos.data[j].z = scaled_pos.z;
            }

        // ensure that the particles are still in their
        // local boxes by wrapping them if they are not
        ArrayHandle<int3> h_image(m_pdata->getImages(),
                                  access_location::host,
                                  access_mode::readwrite);

        const BoxDim& local_box = m_pdata->getBox();

        for (unsigned int i = 0; i < m_pdata->getN(); i++)
            {
            // need to update the image if we move particles from one side
            // of the box to the other
            local_box.wrap(h_pos.data[i], h_image.data[i]);
            }
        }
    }

namespace detail
    {
void export_BoxResizeConstVolumeUpdater(pybind11::module& m)
    {
    pybind11::class_<BoxResizeConstVolumeUpdater, Updater, std::shared_ptr<BoxResizeConstVolumeUpdater>>(
        m,
        "BoxResizeConstVolumeUpdater")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>,
                            std::shared_ptr<Trigger>,
                            Scalar,
                            Scalar,
                            unsigned int,
                            std::shared_ptr<Variant>,
                            std::shared_ptr<ParticleGroup>>())
        .def_property("L1", &BoxResizeConstVolumeUpdater::getL1, &BoxResizeConstVolumeUpdater::setL1)
        .def_property("L2", &BoxResizeConstVolumeUpdater::getL2, &BoxResizeConstVolumeUpdater::setL2)
        .def_property("variant", &BoxResizeConstVolumeUpdater::getVariant, &BoxResizeConstVolumeUpdater::setVariant)
        .def_property_readonly("filter",
                               [](const std::shared_ptr<BoxResizeConstVolumeUpdater> method)
                               { return method->getGroup()->getFilter(); })
        .def("get_current_box", &BoxResizeConstVolumeUpdater::getCurrentBox);
    }

    } // end namespace detail

    } // end namespace hoomd

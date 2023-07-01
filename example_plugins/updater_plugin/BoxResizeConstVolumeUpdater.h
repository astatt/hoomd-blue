// Copyright (c) 2009-2022 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

/*! \file BoxResizeUpdater.h
    \brief Declares an updater that resizes the simulation box of the system
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#include "hoomd/BoxDim.h"
#include "hoomd/ParticleGroup.h"
#include "hoomd/Updater.h"
#include "hoomd/Variant.h"

#include <memory>
#include <pybind11/pybind11.h>
#include <stdexcept>
#include <string>

#ifndef __BOXRESIZECONSTVOLUMEUPDATER_H__
#define __BOXRESIZECONSTVOLUMEUPDATER_H__

namespace hoomd
    {
/// Updates the simulation box over time
/** This simple updater gets the box lengths from specified variants and sets
 * those box sizes over time. As an option, particles can be rescaled with the
 * box lengths or left where they are. Note: rescaling particles does not work
 * properly in MPI simulations.
 * \ingroup updaters
 */
class PYBIND11_EXPORT BoxResizeConstVolumeUpdater : public Updater
    {
    public:
    /// Constructor
    BoxResizeConstVolumeUpdater(std::shared_ptr<SystemDefinition> sysdef,
                     std::shared_ptr<Trigger> trigger,
                     Scalar  L1,
                     Scalar  L2,
                     unsigned int direction,
                     std::shared_ptr<Variant> variant,
                     std::shared_ptr<ParticleGroup> m_group);

    /// Destructor
    virtual ~BoxResizeConstVolumeUpdater();

    /// Get the current m_box2
    Scalar getL1();

    /// Set a new m_box_1
    void setL1(Scalar L1);

    /// Get the current m_box2
    Scalar getL2();

    /// Set a new m_box_2
    void setL2(Scalar L2);

    /// Gets particle scaling filter
    std::shared_ptr<ParticleGroup> getGroup()
        {
        return m_group;
        }

    /// Set the variant for interpolation
    void setVariant(std::shared_ptr<Variant> variant)
        {
        m_variant = variant;
        }

    /// Get the variant for interpolation
    std::shared_ptr<Variant> getVariant()
        {
        return m_variant;
        }

    /// Get the current box for the given timestep
    BoxDim getCurrentBox(uint64_t timestep);

    /// Update box interpolation based on provided timestep
    virtual void update(uint64_t timestep);

    private:
    Scalar m_L1;         ///< C++ box assoc with min
    Scalar m_L2;         ///< C++ box assoc with max
    BoxDim m_orig_box;
    unsigned int m_direction;         ///< C++ box assoc with max
    std::shared_ptr<Variant> m_variant;     //!< Variant that interpolates between boxes
    std::shared_ptr<ParticleGroup> m_group; //!< Selected particles to scale when resizing the box.
    };

namespace detail
    {
/// Export the BoxResizeUpdater to python
void export_BoxResizeConstVolumeUpdater(pybind11::module& m);
    } // end namespace detail
    } // end namespace hoomd
#endif

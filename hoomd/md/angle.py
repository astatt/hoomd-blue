# -- start license --
# Highly Optimized Object-oriented Many-particle Dynamics -- Blue Edition
# (HOOMD-blue) Open Source Software License Copyright 2009-2016 The Regents of
# the University of Michigan All rights reserved.

# HOOMD-blue may contain modifications ("Contributions") provided, and to which
# copyright is held, by various Contributors who have granted The Regents of the
# University of Michigan the right to modify and/or distribute such Contributions.

# You may redistribute, use, and create derivate works of HOOMD-blue, in source
# and binary forms, provided you abide by the following conditions:

# * Redistributions of source code must retain the above copyright notice, this
# list of conditions, and the following disclaimer both in the code and
# prominently in any materials provided with the distribution.

# * Redistributions in binary form must reproduce the above copyright notice, this
# list of conditions, and the following disclaimer in the documentation and/or
# other materials provided with the distribution.

# * All publications and presentations based on HOOMD-blue, including any reports
# or published results obtained, in whole or in part, with HOOMD-blue, will
# acknowledge its use according to the terms posted at the time of submission on:
# http://codeblue.umich.edu/hoomd-blue/citations.html

# * Any electronic documents citing HOOMD-Blue will link to the HOOMD-Blue website:
# http://codeblue.umich.edu/hoomd-blue/

# * Apart from the above required attributions, neither the name of the copyright
# holder nor the names of HOOMD-blue's contributors may be used to endorse or
# promote products derived from this software without specific prior written
# permission.

# Disclaimer

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS ``AS IS'' AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND/OR ANY
# WARRANTIES THAT THIS SOFTWARE IS FREE OF INFRINGEMENT ARE DISCLAIMED.

# IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
# OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# -- end license --

# Maintainer: joaander / All Developers are free to add commands for new features

R""" Angle potentials.

Angles add forces between specified triplets of particles and are typically used to
model chemical angles between two bonds.

By themselves, angles that have been specified in an initial configuration do nothing. Only when you
specify an angle force (i.e. angle.harmonic), are forces actually calculated between the
listed particles.
"""

from hoomd import _hoomd
from hoomd.md import _md
from hoomd.md import force
import hoomd

import math;
import sys;

class coeff:
    R""" Define angle coefficients.

    The coefficients for all angle force are specified using this class. Coefficients are
    specified per angle type.

    There are two ways to set the coefficients for a particular angle potential.
    The first way is to save the angle potential in a variable and call :py:meth:`set()` directly.
    See below for an example of this.

    The second method is to build the coeff class first and then assign it to the
    angle potential. There are some advantages to this method in that you could specify a
    complicated set of angle potential coefficients in a separate python file and import
    it into your job script.

    Example::

        my_coeffs = hoomd.md.angle.coeff();
        my_angle_force.angle_coeff.set('polymer', k=330.0, r=0.84)
        my_angle_force.angle_coeff.set('backbone', k=330.0, r=0.84)

    """

    ## \internal
    # \brief Initializes the class
    # \details
    # The main task to be performed during initialization is just to init some variables
    # \param self Python required class instance variable
    def __init__(self):
        self.values = {};
        self.default_coeff = {}

    ## \var values
    # \internal
    # \brief Contains the vector of set values in a dictionary

    ## \var default_coeff
    # \internal
    # \brief default_coeff['coeff'] lists the default value for \a coeff, if it is set

    ## \internal
    # \brief Sets a default value for a given coefficient
    # \details
    # \param name Name of the coefficient to for which to set the default
    # \param value Default value to set
    #
    # Some coefficients have reasonable default values and the user should not be burdened with typing them in
    # all the time. set_default_coeff() sets
    def set_default_coeff(self, name, value):
        self.default_coeff[name] = value;

    def set(self, type, **coeffs):
        R""" Sets parameters for angle types.

        Args:
            type (str): Type of angle (or a list of type names)
            coeffs: Named coefficients (see below for examples)

        Calling :py:meth:`set()` results in one or more parameters being set for a angle type. Types are identified
        by name, and parameters are also added by name. Which parameters you need to specify depends on the angle
        potential you are setting these coefficients for, see the corresponding documentation.

        All possible angle types as defined in the simulation box must be specified before executing run().
        You will receive an error if you fail to do so. It is not an error, however, to specify coefficients for
        angle types that do not exist in the simulation. This can be useful in defining a potential field for many
        different types of angles even when some simulations only include a subset.

        Examples::

            my_angle_force.angle_coeff.set('polymer', k=330.0, r0=0.84)
            my_angle_force.angle_coeff.set('backbone', k=1000.0, r0=1.0)
            my_angle_force.angle_coeff.set(['angleA','angleB'], k=100, r0=0.0)

        Note:
            Single parameters can be updated. If both ``k`` and ``r0`` have already been set for a particle type,
            then executing ``coeff.set('polymer', r0=1.0)`` will update the value of ``r0`` and leave the other
            parameters as they were previously set.

        """
        hoomd.util.print_status_line();

        # listify the input
        if isinstance(type, str):
            type = [type];

        for typei in type:
            self.set_single(typei, coeffs);

    ## \internal
    # \brief Sets a single parameter
    def set_single(self, type, coeffs):
        type = str(type);

        # create the type identifier if it hasn't been created yet
        if (not type in self.values):
            self.values[type] = {};

        # update each of the values provided
        if len(coeffs) == 0:
            hoomd.context.msg.error("No coefficents specified\n");
        for name, val in coeffs.items():
            self.values[type][name] = val;

        # set the default values
        for name, val in self.default_coeff.items():
            # don't override a coeff if it is already set
            if not name in self.values[type]:
                self.values[type][name] = val;

    ## \internal
    # \brief Verifies that all values are set
    # \details
    # \param self Python required self variable
    # \param required_coeffs list of required variables
    #
    # This can only be run after the system has been initialized
    def verify(self, required_coeffs):
        # first, check that the system has been initialized
        if not hoomd.init.is_initialized():
            hoomd.context.msg.error("Cannot verify angle coefficients before initialization\n");
            raise RuntimeError('Error verifying force coefficients');

        # get a list of types from the particle data
        ntypes = hoomd.context.current.system_definition.getAngleData().getNTypes();
        type_list = [];
        for i in range(0,ntypes):
            type_list.append(hoomd.context.current.system_definition.getAngleData().getNameByType(i));

        valid = True;
        # loop over all possible types and verify that all required variables are set
        for i in range(0,ntypes):
            type = type_list[i];

            if type not in self.values.keys():
                hoomd.context.msg.error("Angle type " +str(type) + " not found in angle coeff\n");
                valid = False;
                continue;

            # verify that all required values are set by counting the matches
            count = 0;
            for coeff_name in self.values[type].keys():
                if not coeff_name in required_coeffs:
                    hoomd.context.msg.notice(2, "Notice: Possible typo? Force coeff " + str(coeff_name) + " is specified for type " + str(type) + \
                          ", but is not used by the angle force\n");
                else:
                    count += 1;

            if count != len(required_coeffs):
                hoomd.context.msg.error("Angle type " + str(type) + " is missing required coefficients\n");
                valid = False;

        return valid;

    ## \internal
    # \brief Gets the value of a single angle potential coefficient
    # \detail
    # \param type Name of angle type
    # \param coeff_name Coefficient to get
    def get(self, type, coeff_name):
        if type not in self.values.keys():
            hoomd.context.msg.error("Bug detected in force.coeff. Please report\n");
            raise RuntimeError("Error setting angle coeff");

        return self.values[type][coeff_name];

    ## \internal
    # \brief Return metadata
    def get_metadata(self):
        return self.values

class harmonic(force._force):
    R""" Harmonic angle potential.

    The command angle.harmonic specifies a harmonic potential energy between every triplet of particles
    with an angle specified between them.

    .. math::

        V(\theta) = \frac{1}{2} k \left( \theta - \theta_0 \right)^2

    where :math:`\theta` is the angle between the triplet of particles.

    Coefficients:

    - :math:`\theta_0` - rest angle  ``t0`` (in radians)
    - :math:`k` - potential constant ``k`` (in units of energy/radians^2)

    Coefficients :math:`k` and :math:`\theta_0` must be set for each type of angle in the simulation using the
    method :py:meth:`angle_coeff.set() <coeff.set()>`.

    Examples::

        harmonic = angle.harmonic()
        harmonic.angle_coeff.set('polymer', k=3.0, t0=0.7851)
        harmonic.angle_coeff.set('backbone', k=100.0, t0=1.0)

    """
    def __init__(self):
        hoomd.util.print_status_line();
        # check that some angles are defined
        if hoomd.context.current.system_definition.getAngleData().getNGlobal() == 0:
            hoomd.context.msg.error("No angles are defined.\n");
            raise RuntimeError("Error creating angle forces");

        # initialize the base class
        force._force.__init__(self);

        # setup the coefficient vector
        self.angle_coeff = coeff();

        # create the c++ mirror class
        if not hoomd.context.exec_conf.isCUDAEnabled():
            self.cpp_force = _md.HarmonicAngleForceCompute(hoomd.context.current.system_definition);
        else:
            self.cpp_force = _md.HarmonicAngleForceComputeGPU(hoomd.context.current.system_definition);

        hoomd.context.current.system.addCompute(self.cpp_force, self.force_name);

        self.required_coeffs = ['k', 't0'];

    ## \internal
    # \brief Update coefficients in C++
    def update_coeffs(self):
        coeff_list = self.required_coeffs;
        # check that the force coefficients are valid
        if not self.angle_coeff.verify(coeff_list):
           hoomd.context.msg.error("Not all force coefficients are set\n");
           raise RuntimeError("Error updating force coefficients");

        # set all the params
        ntypes = hoomd.context.current.system_definition.getAngleData().getNTypes();
        type_list = [];
        for i in range(0,ntypes):
            type_list.append(hoomd.context.current.system_definition.getAngleData().getNameByType(i));

        for i in range(0,ntypes):
            # build a dict of the coeffs to pass to proces_coeff
            coeff_dict = {};
            for name in coeff_list:
                coeff_dict[name] = self.angle_coeff.get(type_list[i], name);

            self.cpp_force.setParams(i, coeff_dict['k'], coeff_dict['t0']);

    ## \internal
    # \brief Get metadata
    def get_metadata(self):
        data = force._force.get_metadata(self)

        # make sure coefficients are up-to-date
        self.update_coeffs()

        data['angle_coeff'] = self.angle_coeff
        return data

class cgcmm(force._force):
    R""" CGCMM angle potential.

    The command angle.cgcmm defines a regular harmonic potential energy between every defined triplet
    of particles in the simulation, but in addition in adds the repulsive part of a CGCMM pair potential
    between the first and the third particle.

    `B. Levine et. al. 2011 <http://dx.doi.org/10.1021/ct2005193>`_ describes the CGCMM implementation details in
    HOOMD-blue. Cite it if you utilize the CGCMM potential in your work.

    The total potential is thus:

    .. math::

        V(\theta) = \frac{1}{2} k \left( \theta - \theta_0 \right)^2

    where :math:`\theta` is the current angle between the three particles
    and either:

    .. math::

        V_{\mathrm{LJ}}(r_{13}) -V_{\mathrm{LJ}}(r_c) \mathrm{~with~~~} V_{\mathrm{LJ}}(r) = 4 \varepsilon \left[
        \left( \frac{\sigma}{r} \right)^{12} - \left( \frac{\sigma}{r} \right)^{6} \right]
        \mathrm{~~~~for~} r <= r_c \mathrm{~~~} r_c = \sigma \cdot 2^{\frac{1}{6}}


    .. math::

        V_{\mathrm{LJ}}(r_{13}) -V_{\mathrm{LJ}}(r_c) \mathrm{~with~~~}
        V_{\mathrm{LJ}}(r) = \frac{27}{4} \varepsilon \left[ \left( \frac{\sigma}{r} \right)^{9} -
        \left( \frac{\sigma}{r} \right)^{6} \right]
        \mathrm{~~~~for~} r <= r_c \mathrm{~~~} r_c = \sigma \cdot \left(\frac{3}{2}\right)^{\frac{1}{3}}

    .. math::

        V_{\mathrm{LJ}}(r_{13}) -V_{\mathrm{LJ}}(r_c) \mathrm{~with~~~}
        V_{\mathrm{LJ}}(r) = \frac{3\sqrt{3}}{2} \varepsilon \left[ \left( \frac{\sigma}{r} \right)^{12} -
        \left( \frac{\sigma}{r} \right)^{4} \right]
        \mathrm{~~~~for~} r <= r_c \mathrm{~~~} r_c = \sigma \cdot 3^{\frac{1}{8}}

    with :math:`r_{13}` being the distance between the two outer particles of the angle.

    Coefficients:

    - :math:`\theta_0` - rest angle ``t0`` (in radians)
    - :math:`k` - potential constant ``k`` (in units of energy/radians^2)
    - :math:`\varepsilon` - strength of potential ``epsilon`` (in energy units)
    - :math:`\sigma` - distance of interaction ``sigma`` (in distance units)

    Coefficients :math:`k, \theta_0, \varepsilon``, and :math:`\sigma` and Lennard-Jones exponents pair must be set for
    each type of angle in the simulation using :py:meth:`set_coeff()`.
    """
    def __init__(self):
        hoomd.util.print_status_line();
        # check that some angles are defined
        if hoomd.context.current.system_definition.getAngleData().getNGlobal() == 0:
            hoomd.context.msg.error("No angles are defined.\n");
            raise RuntimeError("Error creating CGCMM angle forces");

        # initialize the base class
        force._force.__init__(self);

        # create the c++ mirror class
        if not hoomd.context.exec_conf.isCUDAEnabled():
            self.cpp_force = _md.CGCMMAngleForceCompute(hoomd.context.current.system_definition);
        else:
            self.cpp_force = _md.CGCMMAngleForceComputeGPU(hoomd.context.current.system_definition);

        hoomd.context.current.system.addCompute(self.cpp_force, self.force_name);

        # variable for tracking which angle type coefficients have been set
        self.angle_types_set = [];

    def set_coeff(self, angle_type, k, t0, exponents, epsilon, sigma):
        R""" Sets the CG-CMM angle coefficients for a particular angle type.

        Args:
            angle_type (str): Angle type to set coefficients for
            k (float): Coefficient :math:`k` (in units of energy/radians^2)
            t0 (float): Coefficient :math:`\theta_0` (in radians)
            exponents (str): is the type of CG-angle exponents we want to use for the repulsion.
            epsilon (float): is the 1-3 repulsion strength (in energy units)
            sigma (float): is the CG particle radius (in distance units)

        Examples::

            cgcmm.set_coeff('polymer', k=3.0, t0=0.7851, exponents=126, epsilon=1.0, sigma=0.53)
            cgcmm.set_coeff('backbone', k=100.0, t0=1.0, exponents=96, epsilon=23.0, sigma=0.1)
            cgcmm.set_coeff('residue', k=100.0, t0=1.0, exponents='lj12_4', epsilon=33.0, sigma=0.02)
            cgcmm.set_coeff('cg96', k=100.0, t0=1.0, exponents='LJ9-6', epsilon=9.0, sigma=0.3)

        """
        hoomd.util.print_status_line();
        cg_type=0

        # set the parameters for the appropriate type
        if (exponents == 124) or  (exponents == 'lj12_4') or  (exponents == 'LJ12-4') :
            cg_type=2;

            self.cpp_force.setParams(hoomd.context.current.system_definition.getAngleData().getTypeByName(angle_type),
                                     k,
                                     t0,
                                     cg_type,
                                     epsilon,
                                     sigma);

        elif (exponents == 96) or  (exponents == 'lj9_6') or  (exponents == 'LJ9-6') :
            cg_type=1;

            self.cpp_force.setParams(hoomd.context.current.system_definition.getAngleData().getTypeByName(angle_type),
                                     k,
                                     t0,
                                     cg_type,
                                     epsilon,
                                     sigma);

        elif (exponents == 126) or  (exponents == 'lj12_6') or  (exponents == 'LJ12-6') :
            cg_type=3;

            self.cpp_force.setParams(hoomd.context.current.system_definition.getAngleData().getTypeByName(angle_type),
                                     k,
                                     t0,
                                     cg_type,
                                     epsilon,
                                     sigma);
        else:
            raise RuntimeError("Unknown exponent type.  Must be 'none' or one of MN, ljM_N, LJM-N with M/N in 12/4, 9/6, or 12/6");

        # track which particle types we have set
        if not angle_type in self.angle_types_set:
            self.angle_types_set.append(angle_type);

    def update_coeffs(self):
        # get a list of all angle types in the simulation
        ntypes = hoomd.context.current.system_definition.getAngleData().getNTypes();
        type_list = [];
        for i in range(0,ntypes):
            type_list.append(hoomd.context.current.system_definition.getAngleData().getNameByType(i));

        # check to see if all particle types have been set
        for cur_type in type_list:
            if not cur_type in self.angle_types_set:
                hoomd.context.msg.error(str(cur_type) + " coefficients missing in angle.cgcmm\n");
                raise RuntimeError("Error updating coefficients");



def _table_eval(theta, V, T, width):
      dth = (math.pi) / float(width-1);
      i = int(round((theta)/dth))
      return (V[i], T[i])

class table(force._force):
    R""" Tabulated angle potential.

    Args:

        width (int): Number of points to use to interpolate V and F (see documentation above)
        name (str): Name of the force instance

    :py:class:`table` specifies that a tabulated  angle potential should be added to every bonded triple of particles
    in the simulation.

    The torque :math:`T` is (in units of force * distance) and the potential :math:`V(\theta)` is (in energy units):

    .. math::

        T(\theta)     = & T_{\mathrm{user}}(\theta) \\
        V(\theta)     = & V_{\mathrm{user}}(\theta)

    where :math:`\theta` is the angle from A-B to B-C in the triple.

    :math:`T_{\mathrm{user}}(\theta)` and :math:`V_{\mathrm{user}}(\theta)` are evaluated on *width* grid points
    between :math:`0` and :math:`\pi`. Values are interpolated linearly between grid points.
    For correctness, you must specify: :math:`T = -\frac{\partial V}{\partial \theta}`

    Parameters:

    - :math:`T_{\mathrm{user}}(\theta)` and :math:`V_{\mathrm{user}}(\theta)` - evaluated by `func` (see example)
    - coefficients passed to `func` - `coeff` (see example)

    The table *width* is set once when :py:class:`table` is specified. There are two ways to specify the other
    parameters.

    .. rubric:: Set table from a given function

    When you have a functional form for T and F, you can enter that
    directly into python. :py:class:`table` will evaluate the given function over *width* points between :math:`0` and :math:`\pi`
    and use the resulting values in the table::

        def harmonic(theta, kappa, theta_0):
            V = 0.5 * kappa * (theta-theta_0)**2;
            T = -kappa*(theta-theta_0);
            return (V, T)

        btable = angle.table(width=1000)
        btable.angle_coeff.set('angle1', func=harmonic, coeff=dict(kappa=330, theta_0=0))
        btable.angle_coeff.set('angle2', func=harmonic,coeff=dict(kappa=30, theta_0=0.1))

    .. rubric:: Set a table from a file

    When you have no function for for *T* or *F*, or you otherwise have the data listed in a file, :py:class:`table` can use the given
    values directly. You must first specify the number of rows in your tables when initializing :py:class:`table`. Then use
    :py:meth:`set_from_file()` to read the file::

        btable = angle.table(width=1000)
        btable.set_from_file('polymer', 'angle.dat')

    """
    def __init__(self, width, name=None):
        hoomd.util.print_status_line();

        # initialize the base class
        force._force.__init__(self, name);


        # create the c++ mirror class
        if not hoomd.context.exec_conf.isCUDAEnabled():
            self.cpp_force = _md.TableAngleForceCompute(hoomd.context.current.system_definition, int(width), self.name);
        else:
            self.cpp_force = _md.TableAngleForceComputeGPU(hoomd.context.current.system_definition, int(width), self.name);

        hoomd.context.current.system.addCompute(self.cpp_force, self.force_name);

        # setup the coefficent matrix
        self.angle_coeff = coeff();

        # stash the width for later use
        self.width = width;

    def update_angle_table(self, atype, func, coeff):
        # allocate arrays to store V and F
        Vtable = _hoomd.std_vector_scalar();
        Ttable = _hoomd.std_vector_scalar();

        # calculate dth
        dth = math.pi / float(self.width-1);

        # evaluate each point of the function
        for i in range(0, self.width):
            theta = dth * i;
            (V,T) = func(theta, **coeff);

            # fill out the tables
            Vtable.append(V);
            Ttable.append(T);

        # pass the tables on to the underlying cpp compute
        self.cpp_force.setTable(atype, Vtable, Ttable);


    def update_coeffs(self):
        # check that the angle coefficents are valid
        if not self.angle_coeff.verify(["func", "coeff"]):
            hoomd.context.msg.error("Not all angle coefficients are set for angle.table\n");
            raise RuntimeError("Error updating angle coefficients");

        # set all the params
        ntypes = hoomd.context.current.system_definition.getAngleData().getNTypes();
        type_list = [];
        for i in range(0,ntypes):
            type_list.append(hoomd.context.current.system_definition.getAngleData().getNameByType(i));


        # loop through all of the unique type angles and evaluate the table
        for i in range(0,ntypes):
            func = self.angle_coeff.get(type_list[i], "func");
            coeff = self.angle_coeff.get(type_list[i], "coeff");

            self.update_angle_table(i, func, coeff);

    def set_from_file(self, anglename, filename):
        R""" Set a angle pair interaction from a file.

        Args:
            anglename (str): Name of angle
            filename (str): Name of the file to read

        The provided file specifies V and F at equally spaced theta values::

            #t  V    T
            0.0 2.0 -3.0
            1.5707 3.0 -4.0
            3.1414 2.0 -3.0

        Warning:
            The theta values are not used by the code.  It is assumed that a table that has N rows will start at 0, end at :math:`\pi`
            and that :math:`\delta \theta = \pi/(N-1)`. The table is read
            directly into the grid points used to evaluate :math:`T_{\mathrm{user}}(\theta)` and :math:`V_{\mathrm{user}}(\theta)`.

        """
        hoomd.util.print_status_line();

        # open the file
        f = open(filename);

        theta_table = [];
        V_table = [];
        T_table = [];

        # read in lines from the file
        for line in f.readlines():
            line = line.strip();

            # skip comment lines
            if line[0] == '#':
                continue;

            # split out the columns
            cols = line.split();
            values = [float(f) for f in cols];

            # validate the input
            if len(values) != 3:
                hoomd.context.msg.error("angle.table: file must have exactly 3 columns\n");
                raise RuntimeError("Error reading table file");

            # append to the tables
            theta_table.append(values[0]);
            V_table.append(values[1]);
            T_table.append(values[2]);

        # validate input
        if self.width != len(theta_table):
            hoomd.context.msg.error("angle.table: file must have exactly " + str(self.width) + " rows\n");
            raise RuntimeError("Error reading table file");


        # check for even spacing
        dth = math.pi / float(self.width-1);
        for i in range(0,self.width):
            theta =  dth * i;
            if math.fabs(theta - theta_table[i]) > 1e-3:
                hoomd.context.msg.error("angle.table: theta must be monotonically increasing and evenly spaced\n");
                raise RuntimeError("Error reading table file");

        hoomd.util.quiet_status();
        self.angle_coeff.set(anglename, func=_table_eval, coeff=dict(V=V_table, T=T_table, width=self.width))
        hoomd.util.unquiet_status();

    ## \internal
    # \brief Get metadata
    def get_metadata(self):
        data = force._force.get_metadata(self)

        # make sure coefficients are up-to-date
        self.update_coeffs()

        data['angle_coeff'] = self.angle_coeff
        return data

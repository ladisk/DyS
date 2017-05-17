"""
Created on 9. maj. 2013

@author: luka.skrinjar
"""
import itertools
import os
from pprint import pprint
import shutil


import numpy as np
import vtk


from MBD_system.MBD_system_items import BodyItem
from simulation_control_widget.vtk_widget.vtk_coordinate_system.vtk_coordinate_system import VTKCoordinateSystem
from MBD_system.q2dR_i import q2dR_i
from MBD_system.q2R_i import q2R_i
from MBD_system.q2dtheta_i import q2dtheta_i
from MBD_system.q2q_body import q2q_body
from MBD_system.Ai_ui_P import Ai_ui_P_vector
from MBD_system.q2dq_body import q2dq_body


class Body(BodyItem):
    """
    classdocs
    """
    _id = itertools.count(0)

    def __init__(self, name="", file_path=None, parent=None):
        """
        Constructor of body class
        :param name:                    body name (string)
        :param MBD_folder_abs_path:     absolute path to MBD system project folder
        :param density:                 density of the material of the body
        :param volume:                  volume of the body (as float) in m^3
        :param mass:                    mass of the body (as float) in kg
        :param J_zz:                    mass moment of inertia of a body against z-axis (as float) in kg*m^2
        :param u_CAD:                   a vector to mass center of a body in body CAD CS (as array) in m
        :param r_CAD:                   a vector to body CAD CS in GCS of a system (as array) in m
        :param theta:                   orientation angles (numpy array) in degrees
        :param dR:                      a vector of velocities (numpy array) in m/s
        :param dtheta:                  a vector of angular velocities (numpy array) in deg/s
        :param color:                   a color vector (RGB)
        :param transparent:             transparency (float) a value from 0 to 1
        :param visible:                 true or false
        :param display_style:           a display style (as string) options: filled, wireframe, points-check last option
        :param properties_file:         a path to mass and geometry properties data in .dat file (todo)
        :param geometry_data_file:      a path to geometry .stl or .obj file (todo)
        """
        super(Body, self).__init__(name, parent=parent)

        #   parent
        self._parent = parent

        #   id
        self.body_id = None

        #   type of body
        self.body_type = None

        #    name of a body
        self._name = name

        #   option to save results
        self.save_results = False
        self._results_folder_path = None

        #    geometry and physical properties
        self.mass = 0.
        self.J_zz = 0.
        self.q = None
        self.q_i_size = None
        self.q_i_dim = None
        self.M = None
        self.M_size = None

        #   material properties
        self.density = 0.
        self.volume = 0.

        #   coordinate system properties to transform from CAD to LCS of a body
        #   CAD coordinate system in GCS CAD_LCS_GCS -> r_CAD (or r_CADCS_GCS)
        self.r_CAD = np.zeros(3)
        #   LCS (CM - center of mass coordinate system) in CAD coordinate system
        self.u_CAD = np.zeros(3)

        #    dynamic properties
        #    transform with respect to selected CS
        #    options: CAD, LCS
        self.transformCS = "CAD"
        self.R = np.zeros(3)
        self.CM = self.R

        #    (initial) coordinates and angles (in degrees)
        self.theta = np.zeros(3)
        #    (initial) translational and rotational velocities
        self.dR = np.zeros(3)
        self.dtheta = np.zeros(3)

        #   energy
        self._kinetic_energy = 0.
        self._potential_energy = 0.

        #   visualization properties
        self._visible = True
        self.geometry = None
        self.geometry_list = []
        self._idVisible = False
        self.color = np.ones(3, dtype="float32")
        self.transparent = 1.
        self.display_style = "filled"
        self.display_styles = ["filled",
                               "wireframe",
                               "points"]

        #   visualization properties
        self.scale = 1E-3
        self.vtk_reader = None
        self.vtk_mapper = None
        self.vtk_actor = None
        self.vtk_center_of_mass = None

        #    AABB tree object assigned to body object as attribute if body is specified to be in contact with other body
        self.AABBtree = None
        self.AABBtree_created = False

        #    contact properties
        self.max_penetration_depth = 1E-7
        self.uP_i_max = None

        #   list of forces
        self.forces = []

        #   list for markers
        self.markers = []

        #   list of contact geometry data (nodes, edges) objects that are created for contact detection
        self.contact_geometries = []

        #   flexible body properties
        self.mesh = None

        #   default z dimension for planar dimension
        self.z_dim = 0

        #   geometry properties
        self._geometry_type = None
        self._geometry_types = ["line", 
                                ".stl",
                                "lines"] #    TODO: lines - for complex parametric models
        self.geometry = None
        self.geometry_filename = None
        self.geometry_file_extension = ".stl"

        #    file properties
        self._dict = {}
        self.file_path = file_path
        if hasattr(self._parent, "_typeInfo"):
            if self._parent._typeInfo == "MBDsystem":
                #    ground body object has parent MBD system
                self.properties_file_extension = self._parent._project_filetype
            
            elif hasattr(self._parent._parent, "_typeInfo"):
                if self._parent._typeInfo == "group":
                    #    ground body object has parent MBD system
                    self.properties_file_extension = self._parent._parent._data_filetype

        else:
            self.properties_file_extension = ".dat"

        if self._name is not None and self.properties_file_extension is not None:
            self.properties_file_with_extension = self._name + self.properties_file_extension
        else:
            self.properties_file_with_extension = None

        #   create coordinate systems
        #   local coordinate system at center of gravity of a body
        self.LCS = None
        self.LCS_actor = None

        #   coordinate system of geometry
        self.geometry_CS = None
        self.geometry_CS_actor = None

    @classmethod
    def _count(cls):
        return Body._id.next()

    def _excel_header_q(self):
        """

        :return:
        """
        id = str(self.body_id)
        return ["Rx_" + id, "Ry_" + id, "theta_" + id]

    def _excel_header_dq(self):
        """

        :return:
        """
        id = str(self.body_id)
        return ["dRx_" + id, "dRy_" + id, "dtheta_" + id]

    def _save_results(self, folder_path):
        """

        :return:
        """
        #   check if options to save
        if self.save_results:
            if os.path.exists(folder_path):
                shutil.rmtree(folder_path)

            os.makedirs(folder_path)

    def _set_vtk_data(self):
        """
        Define in subclass
        :return:
        """

    def set_vtk_data(self):
        """
        Defined in subclass
        :return:
        """
        #   set vtk data defined in subclass
        self._set_vtk_data()

        #   common settings
        self.vtk_actor.GetProperty().SetColor(self.color)
        self.vtk_actor.GetProperty().SetOpacity(self.transparent)

        self.vtk_actor.SetPosition(self.R)
        self.vtk_actor.SetOrientation(np.rad2deg(self.theta))

        #   marker of LCS
        self.set_vtk_LCS()

        #   marker of geometry CS
        self.set_vtk_geometry_CS()

    def _update_vtk_data(self, t, q):
        """
        Menthod defined in subclass.
        :return:
        """

    def update_vtk_data(self, t, q):
        """
        Defined in subclass
        :param t:
        :param q:
        :return:
        """
        self._update_vtk_data(t, q)

        #   LCS marker
        self.update_vtk_LCS()

    def set_vtk_LCS(self):
        """
        Marker of LCS
        :return:
        """
        self.LCS_actor = VTKCoordinateSystem(parent=self)
        transform = vtk.vtkTransform()
        transform.Translate(self.R)
        theta = np.rad2deg(self.theta)
        transform.RotateZ(theta[2])
        transform.RotateY(theta[1])
        transform.RotateX(theta[0])
        self.LCS_actor.SetUserTransform(transform)

    def update_vtk_LCS(self, R=None, theta=None):
        """

        :return:
        """
        if self.LCS_actor is not None:
            transform = vtk.vtkTransform()
            if R is None:
                R = self.R

            if theta is None:
                theta = self.theta

            if len(R) == 2:
                R = np.append(R, 0.)

            transform.Translate(R)
            theta = np.rad2deg(theta)
            transform.RotateZ(theta[2])
            transform.RotateY(theta[1])
            transform.RotateX(theta[0])
            self.LCS_actor.SetUserTransform(transform)

    def set_vtk_geometry_CS(self):
        """
        Marker of geometry CS
        :return:
        """
        self.geometry_CS_actor = VTKCoordinateSystem(parent=self)
        transform = vtk.vtkTransform()
        transform.Translate(self.r_CAD)
        theta = np.rad2deg(self.theta)
        transform.RotateZ(theta[2])
        transform.RotateY(theta[1])
        transform.RotateX(theta[0])
        self.geometry_CS_actor.SetUserTransform(transform)

    def update_vtk_geometry_CS(self, uR=None, theta=None):
        """

        :return:
        """
        if self.geometry_CS_actor is not None:
            transform = vtk.vtkTransform()
            if uR is None:
                rG = self.r_CAD

            if theta is None:
                theta = self.theta

            if len(rG) == 2:
                rG = np.append(rG, 0.)

            transform.Translate(rG)
            theta = np.rad2deg(theta)
            transform.RotateZ(theta[2])
            transform.RotateY(theta[1])
            transform.RotateX(theta[0])
            self.geometry_CS_actor.SetUserTransform(transform)

    def highlightSelected(self):
        """

        :return:
        """
        self.vtk_actor.GetProperty().SetColor(1.0, 0.0, 0.0)
        self.vtk_actor.GetProperty().SetDiffuse(1.0)
        self.vtk_actor.GetProperty().SetSpecular(0.0)

    def unHighlightSelected(self):
        """

        :return:
        """
        self.vtk_actor.GetProperty().SetColor(self.color)
        self.vtk_actor.GetProperty().SetOpacity(self.transparent)
        self.vtk_actor.GetProperty().SetSpecular(0.0)

    def evaluate_M(self):
        """
        Defined in subclass
        :return:
        """
        return self.M

    def evaluate_CM(self, q):
        """

        :return:
        """
        self.R[0:2] = q2R_i(q, self.body_id)

        return self.R[0:2]

    def print_M(self):
        """

        :return:
        """
        print self.evaluate_M()

    def print_q(self):
        """

        :return:
        """
        print self.evaluate_q()

    def print_dq(self):
        """

        :return:
        """
        print self.evaluate_dq()

    def vtkCreateBoxWidget(self, interactor):
        """

        :param interactor:
        :return:
        """
        boxWidget = vtk.vtkBoxWidget()
        boxWidget.SetInteractor(interactor)
        boxWidget.SetProp3D(self.vtk_actor)
        boxWidget.SetPlaceFactor(1.) # Make the box 1.25x larger than the actor
        boxWidget.PlaceWidget()
        boxWidget.On()

    def __getstate__(self):
        """

        :return:
        """
        if "VBO_created" in self.__dict__:
            self.__dict__["VBO_created"] = False

        return self.__dict__

    def __setstate__(self, state):
        """

        :param state:
        :return:
        """
        self.__dict__ = state

    def add_attributes_from_dict(self, dictionary):
        """

        :return:
        """
        input_angle_units = "rad"
        if hasattr(self._parent, "_parent"):
            if hasattr(self._parent._parent, "_angle_units"):
                input_angle_units = self._parent._parent._input_angle_units

        for key, val in dictionary.iteritems():
            if key == "name":
                key = "_name"

            if key == "theta" and input_angle_units == "deg":
                val = np.deg2rad(val)

            #   cross section properties
            _cross_section_str = "cross_section."
            if self.mesh is not None and _cross_section_str in key:
                _key = key[key.index(_cross_section_str)+len(_cross_section_str)::]
                setattr(self.mesh.cross_section, _key, val)

            #   finite element properties
            _mes_str = "mesh."
            if self.mesh is not None and _mes_str in key:
                _key = key[key.index(_mes_str)+len(_mes_str)::]
                setattr(self.mesh, _key, val)

            setattr(self, key, val)

    def _evaluate_uP_i_max(self):
        """
        Function finds max coordinate of a point of body geometry in x and y axis
        :return:
        """
        self.uP_i_max = np.zeros(3)
        if self.geometry is not None:
            if self._geometry_type == ".stl":
                self.uP_i_max = np.amax(abs(self.geometry.geom_data.vertices), axis=0)
            elif self._geometry_type == "line":
                self.uP_i_max = np.amax(abs(self.geometry.vertices), axis=0)
            else:
                print "uPmax not evaluated"

    def get_uP_i_max(self):
        """

        :return:
        """
        self._evaluate_uP_i_max()
        return self.uP_i_max

    def _show(self):
        """
        
        """
        if self._visible:
            self._visible = False
        else:
            self._visible = True

    def _showID(self):
        """

        :return:
        """
        self.printID()
        if self._idVisible:
            self._idVisible = False
        else:
            self._idVisible = True

    def printID(self):
        """

        :return:
        """
        print "ID =", self.body_id

    def _show_AABB(self):
        """

        :return:
        """
        if self.AABBtree._visible:
            self.AABBtree._visible = False
        else:
            self.AABBtree._visible = True

    def get_q(self):
        """
        Used by user via GUI only!
        :return:
        """
        theta = self.theta[2]
        #    check angle units for display
        if self._parent._parent._angle_units == "deg":
            theta = np.rad2deg(theta)

        print np.append(self.R[0:2], theta)

    def get_dq(self):
        """

        :return:
        """
        print np.append(self.dR[0:2], self.dtheta[2])

    def get_qdq(self):
        """

        :return:
        """
        print np.array([np.append(self.R[0:2], self.theta[2]), np.append(self.dR[0:2], self.dtheta[2])]).flatten()

    def evaluate_M_size(self):
        """
        Function returns size of mass matrix of body
        :return:
        """
        return self.M_size

    def evaluate_q_i_size(self):
        """

        :return:
        """
        return self.q_i_dim

    def evaluate_q0(self):
        """

        :return:
        """
        self.q = self.evaluate_q()
        return self.q

    def evaluate_q(self):
        """
        Used only in program and also during numerical integration.
        :return:
        """
        return np.append(self.R[0:2], self.theta[2])

    def evaluate_dq(self):
        """

        :return:
        """
        return np.append(self.dR[0:2], self.dtheta[2])

    def evaluate_dq0(self):
        """

        :return:
        """
        self.dq = self.evaluate_dq()
        return self.dq

    def evaluate_qdq(self):
        """

        :return:
        """
        print np.array([np.append(self.R[0:2], self.theta[2]), np.append(self.dR[0:2], self.dtheta[2])]).flatten()

    def update_coordinates_and_angles_2D(self, q_i):
        """
        Function updates coordinates Rx, Ry and angle theta and can be displayed in VTK widget
        Args:
            q_i - vector of coordinates of i-th body
        Returns:
            none
        """
        self.R[0:2] = q_i[0:2]
        self.theta[2] = q_i[2]

    def update_velocities_2D(self, dq):
        """
        Function updates velocities (translational and rotational) dRx, dRy and angle dtheta and can be displayed in VTK widget
        Args:
            dq - vector of velocities
        Returns:
            none
        """
        self.dR[0:2] = dq[0:2]
        self.dtheta[2] = dq[2]

    def evaluate_Q_g(self, g):
        """
        Function evaluates gravity force vector for rigid body
        :return:
        """
        Q_g = np.dot(self.M, g)
        return Q_g

    def evaluate_mechanical_energy(self, q=None, gravity=None):
        """
        Function evaluates mechanical energy of the body
        :return:
        """
        self._mechanical_energy = self.evaluate_kinetic_energy(q) + self.evaluate_potential_energy(q, gravity)
        return self._mechanical_energy

    def evaluate_kinetic_energy(self, q):
        """

        :return:
        """
        #   vector of translational velocity
        if q is None:
            dR = self.dR[0:2]
        else:
            if len(q) == self.q_i_size:
                dR = q[0:2]
            else:
                dR = q2dR_i(q, self.body_id)

        #   float of angular velocity
        if q is None:
            dtheta = self.dtheta[2]
        else:
            if len(q) == self.q_i_size:
                dtheta = q[2]
            else:
                dtheta = q2dtheta_i(q, self.body_id)

        #   kinetic energy
        self._kinetic_energy = 0.5 * self.mass * (np.linalg.norm(dR)**2) + .5 * self.J_zz * (dtheta**2)
        return self._kinetic_energy

    def evaluate_potential_energy(self, q, gravity):
        """

        :return:
        """
        #   coordinates vector
        if q is None:
            R = self.R[0:2]
        else:
            R = q[0:2]

        #   gravity
        if gravity is None:
            g = 0
        else:
            g = np.linalg.norm(gravity)

        #   potential energy
        self._potential_energy = (self.mass * g * R[1])
        return self._potential_energy

    def evaluate_q2q_body(self, q):
        """
        Check if input vector of q is equal to mesh vector
        :param q:
        :return:
        """
        #
        if len(q) == self.q_i_size:
            q_b = q
        else:
            q_b = q2q_body(q, self.body_id)

        return q_b

    def evaluate_q2dq_body(self, q):
        """

        :param q:
        :return:
        """
        if len(q) == self.q_i_size:
            dq_b = q
        else:
            dq_b = q2dq_body(q, self.body_id)

        return dq_b

    def evaluate_rCAD(self):
        """
        Function evaluates vector of geometry CS in GCS, based on body LCS position and rotation data (R, theta=
        :return:
        """
        self.r_CAD = self.R - Ai_ui_P_vector(self.u_CAD, self.theta[2])


if __name__ == "__main__":
    a = Body()
    print "update only 1 component"
    a.R[2] = 99
    pprint(vars(a))
    print a.R[0:2]
    print a.R[2]
    print hasattr(a, "AABBtree")

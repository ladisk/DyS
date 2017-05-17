# coding=utf-8
"""

created by: lskrinjar
date of creation: 25/07/2016
time of creation: 22:34
"""
import itertools
import os
from pprint import pprint


import numpy as np
import vtk


import read_body_data_file.read_body_data_file as read_body_data_file
from MBD_system.Ai_ui_P import Ai_ui_P_vector
from MBD_system.extract_from_dictionary_by_string_in_key import extract_from_dictionary_by_string_in_key
from geometry.geometry import Geometry
from geometry.geometry_2D import Geometry2D
from MBD_system.q2R_i import q2R_i
from MBD_system.q2theta_i import q2theta_i
from MBD_system.body.body import Body
from MBD_system.body.geometry.line.line import Line


class RigidBody(Body):
    """
    classdocs
    """

    def __init__(self, name, file_path="", density=0, volume=0, mass=0, J_zz=0, u_CAD=np.zeros(3), r_CAD=np.zeros(3), theta=np.zeros(3), dR=np.zeros(3), dtheta=np.zeros(3), color=np.ones(3,dtype="float32"), _dict={}, connected_to_ground=False, parent=None):
        """
        Constructor of body class
        :param name:                    body name (string)
        :param filename:                absolute path file of body properties
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
        :param properties_file:         a path to mass and geometry properties data in .dat file (todo)
        :param geometry_data_file:      a path to geometry .stl or .obj file (todo)
        """
        super(RigidBody, self).__init__(name=name, file_path=file_path, parent=parent)

        #   type of body
        self.body_type = "rigid body"

        #   body id
        self.body_id = self._count()

        #   body coordinates
        self.q_i_size = 3

        #    geometry and physical properties
        self.mass = mass
        self.J_zz = J_zz
        #   size of mass matrix
        self.M_size = self.q_i_dim = 3

        #   material properties
        self.density = density
        self.volume = volume

        #   visualization properties
        #   size
        self.size = 1.

        #   coordinate system properties
        self.u_CAD = u_CAD
        self.r_CAD = r_CAD

        #    dynamic properties
        #    transform with respect to selected CS
        #    options: CAD, LCS
        self.transformCS = "CAD"
        self.R = self.u_CAD + self.r_CAD

        #    (initial) coordinates and angles (in degrees)
        self.theta = theta
        #    (initial) translational and rotational velocities
        self.dR = dR
        self.dtheta = dtheta

        #   visualization properties
        self.color = color

        #    connected to ground
        self._connected_to_ground = connected_to_ground

        #   set directory to read body data from file
        self.file_path = file_path
        # os.chdir(MBD_folder_abs_path)

        #    read body properties file
        if os.path.isfile(self.file_path):
            self._dict = read_body_data_file.read_body_data_file(self.file_path)
            self.add_attributes_from_dict(self._dict)

            #    check if both files exist
            if not os.path.isfile(self.file_path):
                raise IOError, "Properties file not found!"

            if self._geometry_type == self.geometry_file_extension and self.geometry_filename is not None:
                if not os.path.isfile(self.geometry_filename):
                    print "Geometry file %s not found!"%self.geometry_filename

        #    additional translation due to rotation with respect to CAD CS
        # self.u_CAD[0:2] = Ai_ui_P_vector(self.u_CAD[0:2], 0)#self.theta[2]
        _R = self.u_CAD - Ai_ui_P_vector(self.u_CAD, self.theta[2])#np.zeros(3)#

        #   reevaluate R, based on body data from .dat file
        if all(self.u_CAD == np.zeros(3)) and all(self.r_CAD == np.zeros(3)):
            pass
        else:
            self.R = self.u_CAD + self.r_CAD

        #    create geometry object
        #    read geometry file and save vertices and normals
        if self._parent is not None:
            os.chdir(self._parent._parent.MBD_folder_abs_path)

        if self.geometry_filename is not None:
            if os.path.isfile(self.geometry_filename):
                #   get extension
                self._geometry_type = os.path.splitext(self.geometry_filename)[1]

                if self._geometry_type == ".stl":
                    self.geometry = Geometry(self.geometry_filename, parent=self)

                elif self._geometry_type == ".txt":
                    self.geometry = Geometry2D(self.geometry_filename, parent=self)

                else:
                    raise ValueError, "Object attribute _geometry_type not correct!"

        elif self._geometry_type == "line":
            self.geometry = Line(parent=self)

        else:
            if self.geometry is None:
                print "Body geometry file %s not found! Attribute self.geometry for body %s not created."%(self.geometry_filename, self._name)

        #   add additional attributes to geometry object
        _dict_geometry = extract_from_dictionary_by_string_in_key(_dict, "geometry.")
        if self.geometry is not None and _dict_geometry:
            self.geometry.add_attributes_from_dict(_dict_geometry)

            if (self.r_CAD == np.zeros(3)).all() and (self.u_CAD == np.zeros(3)).all():
                self.r_CAD = self.R

    def _set_vtk_data(self):
        """

        :return:
        """
        if self.geometry_filename is not None:
            if os.path.isfile(self.geometry_filename):
                self.scale = 1E-3

                self.vtk_reader = vtk.vtkSTLReader()
                self.vtk_reader.SetFileName(self.geometry_filename)

                # create a transform that rotates the cone
                transform = vtk.vtkTransform()
                transform.Translate(-self.u_CAD / self.scale)

                #   transform filter
                transformFilter = vtk.vtkTransformPolyDataFilter()
                transformFilter.SetTransform(transform)
                transformFilter.SetInputConnection(self.vtk_reader.GetOutputPort())
                transformFilter.Update()

                #   vtk mapper
                self.vtk_mapper = vtk.vtkPolyDataMapper()
                self.vtk_mapper.SetInputConnection(transformFilter.GetOutputPort())

                #   vtk actor
                self.vtk_actor = vtk.vtkActor()
                self.vtk_actor.SetMapper(self.vtk_mapper)
                self.vtk_actor.SetScale(self.scale)

                self.evaluate_rCAD()
                # print "GetOrigin () =", self.vtk_actor.GetOrigin()
                # print "GetPosition () =", self.vtk_actor.GetPosition()

                # self.vtk_actor.SetOrigin(self.u_CAD)


                # print "GetOrigin () =", self.vtk_actor.GetOrigin()
                # self.vtk_actor.SetPosition(self.rG)
                # print "GetOrigin () =", self.vtk_actor.GetOrigin()
                # self.vtk_actor.SetPosition(self.R)

                # self.vtk_actor.SetPosition(self.R)
                # self.vtk_actor.SetOrientation(np.rad2deg(self.theta))

        elif self._geometry_type == "line":
            self.geometry.set_vtk_data()
            self.vtk_actor = self.geometry.vtk_actor

        elif self._geometry_type == "cube":
            cube = vtk.vtkCubeSource()
            cube.SetXLength(self.size * self.scale)
            cube.SetYLength(self.size * self.scale)
            cube.SetZLength(self.size * self.scale)
            cube.Update()

            #   mapper
            self.vtk_mapper = vtk.vtkPolyDataMapper()
            self.vtk_mapper.SetInputConnection(cube.GetOutputPort())

            #   actor
            self.vtk_actor = vtk.vtkActor()
            self.vtk_actor.SetMapper(self.vtk_mapper)

        else:
            print "Body geometry not known!"

    # def set_vtk_data(self):
    #     """
    #
    #     :return:
    #     """
    #     self._set_vtk_data()
    #
    #     #   marker of LCS
    #     self.set_vtk_LCS()
    #
    #     #   marker of geometry CS
    #     self.set_vtk_geometry_CS()

    def _update_vtk_data(self, t, q):
        """

        :return:
        """
        #   body coordinates R
        self.R[0:2] = q2R_i(q, self.body_id)
        #   body angles theta
        self.theta[-1] = q2theta_i(q, self.body_id)
        if self.vtk_actor is not None:
            self.evaluate_rCAD()

            self.vtk_actor.SetPosition(self.R)
            self.vtk_actor.SetOrientation(np.rad2deg(self.theta))

        for geom in self.geometry_list:
            geom.vtk_actor.SetPosition(self.R)
            geom.vtk_actor.SetOrientation(np.rad2deg(self.theta))

        # #   LCS marker
        # self.update_vtk_LCS()

        #   geometry marker
        self.update_vtk_geometry_CS()

        #   update AABB
        if hasattr(self.AABBtree, "vtk_actor"):
            self.AABBtree.update_vtk_data(q)

    def get_q(self):
        """
        Function returns a vector of absolute coordinated for rigid body
        :return:
        """
        return np.append(self.R[0:2], self.theta[-1])

    def evaluate_M(self):
        """

        :return:
        """
        self.M = np.array([[self.mass, 0, 0],
                      [0, self.mass, 0],
                      [0, 0, self.J_zz]])

        return self.M

    def evaluate_M_size(self):
        """
        Function evaluates and returns size of mass matrix of body
        :return:
        """
        return self.M_size

    def update(self):
        """
        Function updates simulation parameters
        :return:
        """
        q = self.evaluate_q()
        self.AABBtree.update_frame_geometry(q=q)

    def testing(self):
        """

        :return:
        """
        print "geometry CS =", self.r_CAD
        print "LCS =", self.R

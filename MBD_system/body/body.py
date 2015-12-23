'''
Created on 9. maj. 2013

@author: reti-luka
'''
import time
import sys
import os
from pprint import pprint
import itertools
import logging

import numpy as np
from PyQt4 import QtCore, QtGui
from OpenGL import *
from OpenGL.GL import *
from OpenGL.GLU import *


from geometry.geometry import Geometry
import read_body_data_file.read_body_data_file as read_body_data_file
from simulation_control_widget.opengl_widget.coordinate_system.coordinate_system import CoordinateSystem
from simulation_control_widget.opengl_widget.marker.marker import Marker


from MBD_system.Ai_ui_P import Ai_ui_P_vector
from MBD_system.MBD_system_items import BodyItem


class Body(BodyItem):
    """
    classdocs
    """
    __id = itertools.count(0)

    def __init__(self, body_name, MBD_folder_abs_path=None, density=0, volume=0, mass=0, J_zz=0, CM_CAD_LCS=np.zeros(3), CAD_LCS_GCS=np.zeros(3), theta=np.zeros(3), dR=np.zeros(3), dtheta=np.zeros(3), color_GL=np.ones(3), transparent_GL=1, visible=True, display_style="filled", connected_to_ground = False, parent=None):
        super(Body, self).__init__(body_name, parent)
        """
        Constructor of body class
        Args:
            _name - body name as string
            properties_file - mass and geometry properties data in .txt file
            geometry_data_file - geometry .stl or .obj file
            
            R - position numpy array
            theta - orientation angles in degrees numpy array
            dR - velocitiy numpy array
            dtheta - angular velocity numpy array
            
            color_GL - RGB numpy array
            transparent_GL - a transparency factor of a body between 0 and 1
            display_style - filled or wireframe as string
            
            density - specific density of a body material
        Returns:
            body object with its properties
            OpenGL object - VBO
        """
        self._parent = parent

        #    set working directory
        if body_name.lower() == "ground":
            self.body_id = -1
        else:
            #    body id
            self.body_id = self.__id.next()
            # self.body_id = (len(self._parent._children) - 1)


            if MBD_folder_abs_path != None:
                os.chdir(MBD_folder_abs_path)
            else:
                print "MBS data folder not found. Define MBS data folder."
#             raise IOError, "MBS data folder not found. Define MBS data folder."
        

#         self.body_id = (len(bodies_list) - 1) + 1
        
        #    name as string
        self._name = body_name                                  
        
        #    check if body is ground
        if not body_name == "ground":
            self.mass = mass
            self.J_zz = J_zz
        else:
            #    override default geometry and physical properties
            self.mass = 0
            self.J_zz = 0
            
        self.density = density
        self.volume = volume
        self.CM_CAD_LCS = CM_CAD_LCS
        self.CAD_LCS_GCS = CAD_LCS_GCS

        #    dynamic properties
        self.R = self.CM_CAD_LCS + self.CAD_LCS_GCS
        #    (initial) coordinates and angles (in degrees)
        self.theta = theta
        #    (initial) translational and rotational velocities
        self.dR = dR
        self.dtheta = dtheta

        #    connected to ground
        self._connected_to_ground = connected_to_ground

        #   list of forces
        self.forces = []
        #   list for markers
        self.markers = []
        #   list of contact geometry data (nodes, edges) objects that are created for contact detection
        self.contact_geometries = []

        #   opengl visualization properties
        self._visible = visible
        self.geom = None
        self.VBO_created = False

        #    AABB tree object assigned to body object as attribute if body is specified to be in contact with other body
        self.AABBtree = None
        self.AABBtree_created = False
        
        #    file properties
        if hasattr(self._parent, "_typeInfo"):  # == 
            if self._parent._typeInfo == "MBDsystem":  #    ground body object has parent MBD system
                self.properties_file_extension = self._parent._project_filetype
            
            elif hasattr(self._parent._parent, "_typeInfo"):
                if self._parent._typeInfo == "group":  #    ground body object has parent MBD system
                    self.properties_file_extension = self._parent._parent._data_filetype
        else:
            print 'Data file not found for body %s' % self._name

        self.geometry_file_extension = ".stl"

        self.properties_file_with_extension = self._name + self.properties_file_extension

        #   opengl properties
        #   create coordinate systems
        #   local coordinate system at center of gravity of a body
        self.LCS = CoordinateSystem(parent=self)
        self.LCS._visible = False

        if body_name.lower() == "ground":
            pass
        else:
            self.color_GL = color_GL
            self.transparent_GL = transparent_GL
            self.display_style = display_style

            #    contact properties
            self.max_penetration_depth = 1E-7
            self.uP_i_max = None

            #    if body is not ground read body properties data
            if body_name == "ground":
                self._visible = False
                self.actual_body_name = "ground"
            else:
                #    read body properties file
                if os.path.isfile(self.properties_file_with_extension):
                    # self.actual_body_name, self.density, self.volume, self.mass, self.J_zz, self.CM_CAD_LCS, self.CAD_LCS_GCS, self.theta, self.dR, self.dtheta, self.body_geometry_filename, self.color_GL, self.transparent_GL, self.display_style = read_body_data_file.read_data(self.properties_file_with_extension)
                    _dict = read_body_data_file.read_data(self.properties_file_with_extension)
                    self.add_attributes_from_dict(_dict)

                    #    check if both files exist
                    if not os.path.isfile(self.properties_file_with_extension):
                        raise IOError, "Properties file not found!"

                    if not os.path.isfile(self.geometry_filename):
                        raise IOError, "Geometry file %s not found!"%self.geometry_filename

            self.CM_CAD_LCS_ = self.CM_CAD_LCS
            self.CM_CAD_LCS[0:2] = Ai_ui_P_vector(self.CM_CAD_LCS[0:2], self.theta[2])  # np.deg2rad(self.theta[2])

            self.R = self.CM_CAD_LCS + self.CAD_LCS_GCS

            #   cad coordinate system of geometry
            # print "self.CM_CAD_LCS =", self.CM_CAD_LCS
            self.CAD_CS = Marker(-self.CM_CAD_LCS, parent=self)#-self.CM_CAD_LCS-self.CAD_LCS_GCS
            self.CAD_CS._visible = False

            #    create geom object
            #    read geometry file and save vertices and normals
            if os.path.isfile(self.geometry_filename):
                self.geom = Geometry(self.geometry_filename, parent=self)
            else:
                print "Body geometry file not found! Attribute self.geom not created."

        #    construct and display body geometry from stl data file
        # if body_name != "ground":
        #     #    check if opengl is running without errors and creates VBOs
        #     try:
        #         #if glGetError() == GL_NO_ERROR:
        #         self.create_VBO()
        #
        #     except:
        #         self.create_VBO()
        #         ValueError
        #         logging.getLogger("DyS_logger").error("OpenGL error - is geometry created and displayed?")
        #         raise

    def add_attributes_from_dict(self, dict):
        """

        :return:
        """
        for key in dict:
            val = dict[key]
            if key == "theta":
                val = np.deg2rad(val)

            setattr(self, key, val)

    def _show(self):
        """
        
        """
        if self._visible:
            self._visible = False
        else:
            self._visible = True

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

        :return:
        """
        print np.append(self.R[0:2], self.theta[2])

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

    def _update_VBO(self):
        """
        
        """
        self.geom._update_VBO_color(self.color_GL, self.transparent_GL)

    def create_VBO(self):
        """
        Method is called after the opengl is initialised for each body to construct opengl VBO.
        """
        # print "=================================="
        # print "create_VBO() =", self._name, "id =", self.body_id
        #    construct a body shape OpenGL object - VBO
        if self.geom is not None:
            self.geom.create_VBO()
            #    is buffer - status if VBO is created and can be displayed
            self.VBO_created = True
            #   create body LCS as vbo
            self.LCS._create_VBO()
            # #    create body CAD LCS vbo
            self.CAD_CS._create_VBO()

        # print "self.markers =", self.markers
        # for marker in self.markers:
            # pprint(vars(marker))
            # marker._create_VBO()
            # print "node =", marker.node
        # marker = self.markers[1]
        # if self.body_id == 0:
        #     print "node =", marker.node
        #     marker._create_VBO()

        for contact_geometry in self.contact_geometries:
            contact_geometry.create_VBO()
        #
        if not self.AABBtree_created and self.AABBtree is not None:
            self.AABBtree.create_VBO_tree()
            self.AABBtree_created = True
            self.AABBtree._visible = False
            self.AABBtree.create_VBO()

    def _paintGL_VBO_AABBtree(self, shader_program):
        """
        
        """
        if self.AABBtree is not None and self.AABBtree._visible and self.AABBtree._VBO_created:
            shader_program.bind()
            color_location = shader_program.uniformLocation("vertex_color")
            shader_program.setUniformValue(color_location, QtGui.QVector3D(self.AABBtree.color_GL[0], self.AABBtree.color_GL[1], self.AABBtree.color_GL[2]))  # cg.cgtypes.vec3(body.AABBtree.color_GL) ctypes.c_floatQVector3D(1, 0, 0)
            self.AABBtree.paintGL_VBO_tree()

        glUseProgram(0)

    def paintGL_VBO(self, step=None, shader=None):
        """
        Paint body VBO
        """
        glTranslatef(self.R[0], self.R[1], self.R[2])
        glRotatef(np.rad2deg(self.theta[2]), 0, 0, 1)
        glRotatef(np.rad2deg(self.theta[0]), 1, 0, 0)
        glRotatef(np.rad2deg(self.theta[1]), 0, 1, 0)

        self._paintGL_VBO_AABBtree(shader)
        
        for contact_geometry in self.contact_geometries:
            contact_geometry.paintGL_VBO()

        if self._visible:

            for force in self.forces:
                force._paint_GL(step)

            if self.LCS._visible and self.LCS._VBO_created:
                self.LCS._paintGL_VBO()

            if hasattr(self, "CAD_CS"):
                if self.CAD_CS._visible and self.CAD_CS._VBO_created:
                    self.CAD_CS._paintGL_VBO()

            #   display markers of a body
            for marker in self.markers:
                if marker._VBO_created and marker._visible:
                    try:
                        marker._paintGL_VBO()
                    except:
                        pass
                else:
                    marker._create_VBO()

            if self.VBO_created:
#                 print "buffer check(body) =", glIsBuffer(self.geom.vbo_id)
                #   bind to buffer
                glBindBuffer(GL_ARRAY_BUFFER, self.geom.vbo_id)

                #    pointer - offset in bits
                v_pointer = ctypes.c_void_p(0)
                n_pointer = ctypes.c_void_p(12)
                c_pointer = ctypes.c_void_p(24)

                #    stride in bits (1 float = 4 bits)
                stride_in_bits = 40

                glVertexPointer(3, GL_FLOAT, stride_in_bits, v_pointer)
                glColorPointer(4, GL_FLOAT, stride_in_bits, c_pointer)
                glNormalPointer(GL_FLOAT, stride_in_bits, n_pointer)

                #    draw (fill or wireframe)
                if self.display_style == "wireframe":
                    glDisable(GL_LIGHTING)
                    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)  # GL_FRONT, GL_FRONT_AND_BACK
                    glDrawArrays(GL_TRIANGLES, 0, self.geom.N_vertices)
                    glEnable(GL_LIGHTING)
                elif self.display_style == "filled":
                    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
                    glDrawArrays(GL_TRIANGLES, 0, self.geom.N_vertices)
                elif self.display_style == "points":
                    glDisable(GL_LIGHTING)
                    glPointSize(0.2)
                    glPolygonMode(GL_FRONT_AND_BACK, GL_POINT)
                    glDrawArrays(GL_TRIANGLES, 0, self.geom.N_vertices)
                    glPointSize(1)
                    glEnable(GL_LIGHTING)
                else:
                    None

                glBindBuffer(GL_ARRAY_BUFFER, 0)

        # for contact_geometry in self.contact_geometries:
        #     contact_geometry.paintGL_VBO()

                # for marker in self.markers:
                #     if marker._VBO_created and marker._visible:
                #             marker._paintGL_VBO()
                #     else:
                #         marker._create_VBO()

    def update_coordinates_and_angles_2D(self, q_i):
        """
        Function updates coordinates Rx, Ry and angle theta and can be displayed in opengl widget
        Args:
            q_i - vector of coordinates of i-th body
        Returns:
            none
        """
        self.R[0:2] = q_i[0:2]
        self.theta[2] = q_i[2]

    def update_velocities_2D(self, dq):
        """
        Function updates velocities (translational and rotational) dRx, dRy and angle dtheta and can be displayed in opengl widget
        Args:
            dq - vector of velocities
        Returns:
            none
        """
        self.dR[0:2] = dq[0:2]
        self.dtheta[2] = dq[2]

    def delete_VBO(self):
        self.geom.__del__()
#        if self.AABBtree != None:
#            for _AABBtree in self.AABBtree.children:
#                _AABBtree.__del__()

    def mechanical_energy(self, q=None, dq=None, gravity=None):
        """
        Function evaluates mechanical energy of the body
        :return:
        """
        if q is None:
            R = self.R[0:2]
        else:
            R = q[0:2]

        if dq is None:
            dR = self.dR[0:2]
            dtheta = self.dR[2]
        else:
            dR = dq[0:2]
            dtheta = dq[2]

        if gravity is None:
            g = 0
        else:
            g = gravity

        #   mechanical energy (kinematic and potential)
        _energy = 0.5 * (self.mass * (np.linalg.norm(dR)**2) + self.J_zz * (dtheta**2)) + (self.mass * gravity * R[1])

        return _energy

if __name__ == "__main__":
    a = Body()
    print "update only 1 component"
    a.R[2] = 99
    pprint(vars(a))
    print a.R[0:2]
    print a.R[2]
    print hasattr(a, "AABBtree")

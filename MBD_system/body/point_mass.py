# coding=utf-8
"""

created by: lskrinjar
date of creation: 09/03/2017
time of creation: 22:53
"""
import itertools
import os
from pprint import pprint


import numpy as np
import vtk


from MBD_system.body.body import Body
import read_body_data_file.read_body_data_file as read_body_data_file
from MBD_system.q2R_i import q2R_i


class PointMass(Body):
    """
    classdocs
    """

    def __init__(self, name, file_path="", mass=0., R=np.zeros(3, dtype=float), dR=np.zeros(3, dtype=float), color=np.ones(3, dtype="float32"), _dict={}, connected_to_ground=False, parent=None):
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
        super(PointMass, self).__init__(name=name, file_path=file_path, parent=parent)

        #   type of body
        self.body_type = "point mass"

        #   body id
        self.body_id = self._count()

        #   body coordinates
        self.q_i_dim = 2

        #    geometry and physical properties
        self.mass = mass
        #   size of mass matrix
        self.M_size = self.q_i_dim = 2

        #    dynamic properties
        #    transform with respect to selected CS
        #    options: CAD, LCS
        self.transformCS = "CAD"
        self.R = R
        #    (initial) translational and rotational velocities
        self.dR = dR

        #   sphere radius
        self.radius = 0.5

        #   scale
        self.scale = 1E-3

        #    read body properties file
        if os.path.isfile(self.file_path):
            self._dict = read_body_data_file.read_body_data_file(self.file_path)
            self.add_attributes_from_dict(self._dict)

    def _excel_header_q(self):
        """

        :return:
        """
        id = str(self.body_id)
        return ["Rx_" + id, "Ry_" + id]

    def _excel_header_dq(self):
        """

        :return:
        """
        id = str(self.body_id)
        return ["dRx_" + id, "dRy_" + id]

    def _set_vtk_data(self):
        """

        :return:
        """
        if self.geometry_filename is None:
            self.mass_center = vtk.vtkSphereSource()
            self.mass_center.SetCenter(0., 0., 0.)
            self.mass_center.SetRadius(self.radius * self.scale)

            #   vtk mapper
            self.vtk_mapper = vtk.vtkPolyDataMapper()
            self.vtk_mapper.SetInputConnection(self.mass_center.GetOutputPort())

        else:
            print "TODO@",__name__

        self.vtk_actor = vtk.vtkActor()
        self.vtk_actor.SetMapper(self.vtk_mapper)
        self.vtk_actor.SetPosition(self.R)

    def _update_vtk_data(self, t, q):
        """
        Update of properties and attributes used by this class
        :return:
        """
        #   body coordinates R
        self.R[0:2] = q2R_i(q, self.body_id)

        self.vtk_actor.SetPosition(self.R)

    def evaluate_M(self):
        """

        :return:
        """
        self.M = np.array([[self.mass, 0.],
                           [0., self.mass]])

        return self.M

    def evaluate_Q_g(self, g):
        """
        Function evaluates gravity force vector for rigid body
        :return:
        """
        Q_g = np.dot(self.M, g[0:2])
        return Q_g

    def update_coordinates_and_angles_2D(self, q_i):
        """
        Function updates coordinates Rx, Ry
        Args:
            q_i - vector of coordinates of i-th body
        Returns:
            none
        """
        self.R[0:2] = q_i[0:2]

    def update_velocities_2D(self, dq):
        """
        Function updates velocities (translational only) dRx, dRy
        Args:
            dq - vector of velocities
        Returns:
            none
        """
        self.dR[0:2] = dq[0:2]

    def evaluate_q(self):
        """

        :return:
        """
        return self.R[0:2]

    def evaluate_dq(self):
        """

        :return:
        """
        return self.dR[0:2]


if __name__ == "__main__":
    mass = PointMass()
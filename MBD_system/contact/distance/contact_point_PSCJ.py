# coding=utf-8
"""

created by: lskrinjar
date of creation: 23/01/2017
time of creation: 08:23
"""
import os
import inspect
import time
from pprint import pprint
import numpy as np
from matplotlib import pyplot as plt


from MBD_system.contact.distance.distance_line_node import DistanceLineNode
from MBD_system.transform_cs import gcs2cm_lcs, cm_lcs2gcs, uP_lcs2gcs
from MBD_system.q2R_i import q2R_i
from MBD_system.q2theta_i import q2theta_i
from MBD_system.q2dR_i import q2dR_i
from MBD_system.q2dtheta_i import q2dtheta_i
from MBD_system.n2t import n2t
from MBD_system.dr_contact_point_uP import dr_contact_point_uP
from MBD_system.force.force import Force
from MBD_system.contact.distance.distance_PSCJ import DistancePSCJ


class ContactPointPSCJ(DistancePSCJ):
    """
    classdocs
    """

    def __init__(self, dist_obj, scale=1E-3, parent=None):
        """
        :param r_iP:        center point on body i - slot
        :param r_iR:        center point on body i - slot
        :param r_jP:        center point of pin on body j - pin
        :param tangent:     tangent in GCS
        """
        super(ContactPointPSCJ, self).__init__(dist_obj.r_iP, dist_obj.r_iR, dist_obj.r_jP, parent=parent)

        #   assign all distance object attributes to this object attributes
        if hasattr(dist_obj, "__dict__"):
            self.__dict__ = dist_obj.__dict__

        self.__attrs = ["_dq_n", "_dq_t", "_dq0_n", "_dq0_t"]
        for attr in self.__attrs:
            if not hasattr(self, attr):
                setattr(self, attr, None)

        #   contact point velocity
        self._dq_n = None
        self._dq_t = None
        #   initial
        self._dq0_n = None
        self._dq0_t = None

        #   status
        self.active = True

        #   visualization properties
        self.scale = scale

        #   coordinates
        self.u_P_LCS_list = [np.zeros(2, dtype="float32"), np.zeros(2, dtype="float32")]
        self.r_P_GCS_list = [np.zeros(2, dtype="float32"), np.zeros(2, dtype="float32")]

        #   contact points
        self.u_CP_list = [np.zeros(2, dtype="float32"), np.zeros(2, dtype="float32")]
        self.r_CP_list = [np.zeros(2, dtype="float32"), np.zeros(2, dtype="float32")]

        #   normals and tangents in GCS
        self._n_GCS_list = [np.zeros(2, dtype="float32"), np.zeros(2, dtype="float32")]
        self._t_GCS_list = [np.zeros(2, dtype="float32"), np.zeros(2, dtype="float32")]

        #   contact geometry
        self._n_GCS = None
        self._t_GCS = None

        #   evaluate contact geometry in LCS
        self.contact_geometry_LCS()

        #   contact force
        self.Fn = 0.
        self.Ft = 0.
        self.F = np.zeros(2)
        self._Fn_list = self._create_forces(direction="normal")
        self._Ft_list = self._create_forces(direction="tangent")

    def _create_forces(self, direction=""):
        """
        Method creates forces objects for each body in contact
        :return:
        """
        F_list = []
        if hasattr(self._parent, "body_id_list"):
            #   create force object
            if direction == "normal":
                d = "n"
            elif direction == "tangent":
                d = "t"
            else:
                d = direction

            # for body_id, visibility in zip(self._parent.body_id_list, self._parent._visible_F_list):
            for body_id, visibility in zip(self.body_id_list, self._parent._visible_F_list): #  corrected from self._parent.body_id_list, self._parent._visible_F_list to self.body_id_list, self._parent._visible_F_list

                F = Force(body_id, force_name=self._parent._name + "_F" + d + "_on_body_" + str(body_id), visible=visibility, scale=self.scale, parent=self._parent)
                F_list.append(F)

                #   create force visualization data
                ren = self._parent._parent._parent._parent.dys.simulation_control_widget.vtkWidget.renderer
                F.set_vtk_data(renderer=ren)

                #   append force object to list of forces of MBD system object
                self._parent._parent._parent.forces.append(F)

                #   append force object to list of forces of body object
                # self._parent._parent._parent.bodies[body_id].forces.append(F)

        return F_list

    def _delete_forces(self, F_list):
        """
        Delete force objects
        :return:
        """
        for F in F_list:
            del F

    def _remove_forces(self, F_list):
        """
        Remove force objects from lists
        :return:
        """
        print "F_list =", F_list
        print "self._parent.body_id_list =", self._parent.body_id_list
        for F, body_id in zip(F_list, self._parent.body_id_list):
            print "F =", F
            if F in self._parent._parent._parent.forces:
                print "force removed!"
                self._parent._parent._parent.forces.remove(F)

            if F in self._parent._parent._parent.bodies[body_id].forces:
                self._parent._parent._parent.bodies[body_id].forces.remove(F)

    def _clear_VTK_data(self, F_list):
        """

        :param F_list:
        :return:
        """
        for F in F_list:
            if F.vtk_actor is not None and F.renderer is not None:
                F.remove_vtk_data()

    def deactivate_forces(self):
        """

        :return:
        """
        #   join both lists to one
        F_list = self._Fn_list + self._Ft_list

        for i, F in enumerate(F_list):
            # print "i =", i, F._name
            F.deactivate()

            # F._clear_vtk_data()

        self._Fn_list = []
        self._Ft_list = []

    def set_dq0(self, dq0_n=None, dq0_t=None):
        """
        Function saves the initial normal and tangential contact velocity at impact
        :param dq0_n    a normal contact velocity (scalar, float value)
        :param dq0_t    a tangential contact velocity (scalar, float value)
        """
        print "set_dq0()"
        if dq0_n is None:
            if self._dq0_n is None:
                self._dq0_n = self._dq_n
        else:
            self._dq0_n = dq0_n

        if dq0_t is None:
            if self._dq0_t is None:
                self._dq0_t = self._dq_t
        else:
            self._dq0_t = dq0_t

        print "self._dq0_n =", self._dq0_n
        print "self._dq0_t =", self._dq0_t

        if self._dq0_n < 0.:
            print ValueError, "Normal component of initial relative contact velocity is negative!"
            print "pin in section =", self.pin_in_section
            print "n =", self._n_GCS
            print "t =", self._t_GCS

    def set_dq(self, dq_n=None, dq_t=None):
        """

        :param dq_n:
        :param dq_t:
        :return:
        """
        self._dq_n = dq_n
        self._dq_t = dq_t

    def evaluate_F(self, Fn, Ft):
        """

        Returns:
        """
        self.Fn = Fn
        self.Ft = Ft
        self.F = self.Fn * self.normal + self.Ft * self.tangent

    def contact_geometry_LCS(self):
        """

        :param q:
        :return:
        """
        n = self.n_list[self.index]
        t = self.t_list[self.index]
        for i, (body_id, sign) in enumerate(zip(self.body_id_list, self.sign_list)):
            self._n_GCS_list[i] = sign * n
            self._t_GCS_list[i] = sign * t

    def evaluate_index(self):
        """

        :return:
        """
        if self.pin_in_section == "iPiR":
            index = self._distance[0:2].index(min(self._distance[0:2]))

        if self.pin_in_section == "iP":
            index = 2

        if self.pin_in_section == "iR":
            index = 3

        return index

    def update_contact_point(self, r_iP, r_iR, r_jP):
        """

        :return:
        """
        # print "update_contact_point()@",__name__
        #   coordinates in GCS
        self.r_iP = r_iP
        self.r_iR = r_iR
        self.r_jP = r_jP

        self.r_iPR_list = [self.r_iP, self.r_iR]

        #   distance vector
        self._distance_vector = self._evaluate_distance_vector()

        #   evaluate pin position in a slot
        self.pin_in_section = self.evaluate_pin_in_section()

        #   distance and distance sign
        self._evaluate_distance()

        #   evaluate index
        self.index = self.evaluate_index()

    def contact_points_LCS(self, q):
        """

        :param q:
        :return:
        """
        # print "contact_points_LCS()"
        for i, (body_id, rP, Fn, Ft) in enumerate(zip(self.body_id_list, self.r_CP_list, self._Fn_list, self._Ft_list)):
            R_i = q2R_i(q, body_id)
            theta_i = q2theta_i(q, body_id)
            uP = gcs2cm_lcs(rP, R=R_i, theta=theta_i)

            self.u_P_LCS_list[i] = uP

            #   update position vectors of force object at contact point
            Fn._update_u_P(uP)
            Ft._update_r_P(rP)

    def contact_points_GCS(self, frame_nodes_GCS=np.zeros([4, 2])):
        """

        :return:
        """
        # print "contact_points_GCS()@",__name__
        # print "self.pin_in_section =", self.pin_in_section, "self.index =", self.index
        if self.pin_in_section == "iPiR":
            if self.index == 0:
                r_frame = frame_nodes_GCS[0]

            if self.index == 1:
                r_frame = frame_nodes_GCS[2]

            # print "self._contact_point_obj.index =", self._contact_point_obj.index
            # print "r_frame =", r_frame
            # print "self._contact_point_obj.r_jPiL_list[self._contact_point_obj.index] =", self._contact_point_obj.r_jPiL_list[self._contact_point_obj.index]
            r_Pi = r_frame + np.dot(self.r_jPiL_list[self.index], self._t_GCS_list[self.index]) * self._t_GCS_list[self.index]
            r_Pj = self.r_jP - self.R0_j * self._n_GCS_list[1]

        if self.pin_in_section in ["iP", "iR"]:
            rPi = self.r_iPR_list[["iP", "iR"].index(self.pin_in_section)]
            rPj = self.r_jP

            e = rPj - rPi
            # print "self._n_GCS_list =", self._n_GCS_list
            #   contact point on cylindrical section of slot
            r_Pi = rPi + self.R0_i * self._n_GCS_list[0]

            #   contact point on pin
            r_Pj = rPj + self.R0_j * self._n_GCS_list[0]

        #   list of actual contact points in GCS on each body
        self.r_CP_list = [r_Pi, r_Pj]
        # print "self.r_CP_list =", self.r_CP_list
        # self.contact_points_LCS(q)

    def set_r_P_GCS_list(self, r_P_GCS_list):
        """

        :param r_P_GCS_list:
        :return:
        """
        self.r_P_GCS_list = r_P_GCS_list

    def contact_velocity(self, q):
        """
        positive value - bodies are approaching
        negative values - bodies are moving apart
        :return:
        """
        dr_P = [np.zeros(2, dtype="float32"), np.zeros(2, dtype="float32")]

        for i, (body_id, u_P) in enumerate(zip(self.body_id_list, self.u_P_LCS_list)):
            #   body velocity, R, theta
            dR = q2dR_i(q, body_id)
            theta = q2theta_i(q, body_id)
            #    dtheta - omega
            dtheta = q2dtheta_i(q, body_id)
            #    point velocity
            dr_P_i = dr_contact_point_uP(dR, theta, dtheta, u_P)
            #    update list
            dr_P[i] = dr_P_i

        #    relative contact velocity vector
        _dq = dr_P[1] - dr_P[0]

        #   relative contact velocity
        #   normal direction
        self._n_GCS = self.n_list[self.index]
        self._dq_n = np.dot(_dq, self._n_GCS)

        #   tangent direction
        self._t_GCS = self.t_list[self.index]
        self._dq_t = np.dot(_dq, self._t_GCS)

        return self._dq_n, self._dq_t


if __name__ == "__main__":
    r_iP = np.array([0., 0.])
    r_iR = np.array([1E-3, 0.])
    r_jP = np.array([0.5E-3, -.2E-3])

    L = 0.95E-3
    R0j = 0.75E-3  # 0.1
    hi = 1.6E-3  # 0.3
    R0i = hi / 2.
    c = R0i - R0j

    d = DistancePSCJ()

    cp = ContactPointPSCJ()





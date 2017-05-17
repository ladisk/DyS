"""

created by: lskrinjar
date of creation: 21/02/2017
time of creation: 21:38
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
from MBD_system.contact.distance.distance_RCJ import DistanceRCJ


class ContactPointRCJ(DistanceRCJ):
    """
    classdocs
    """

    def __init__(self, dist_obj, scale=1E-3, parent=None):
        """
        :param r_iP:        free a point in GCS on body i
        :param r_jP:        line point in GCS on body j
        :param normal:      normal in GCS
        :param r_jR:        line point in GCS on body j
        :param tangent:     tangent in GCS
        """
        super(ContactPointRCJ, self).__init__(dist_obj.r_iP, dist_obj.r_jP, parent=parent)

        #   assign all distance object attributes to this object attributes
        if hasattr(dist_obj, "__dict__"):
            self.__dict__ = dist_obj.__dict__

        #   contact point velocity
        self._dq_n = None
        self._dq_t = None
        #   initial
        self._dq0_n = None
        self._dq0_t = None

        #   status of contact point
        self.active = True

        #   visualization properties
        self.scale = scale

        #   distance with sign
        if dist_obj._distance > parent._radial_clearance:
            self._distance_sign = -dist_obj._distance

        #   coordinates in GCS
        self.r_P_GCS_list = [np.zeros(2, dtype="float32"), np.zeros(2, dtype="float32")]
        self.u_P_LCS_list = [np.zeros(2, dtype="float32"), np.zeros(2, dtype="float32")]

        #   normals and tangents in GCS
        self._n_GCS_list = [np.zeros(2, dtype="float32"), np.zeros(2, dtype="float32")]
        self._t_GCS_list = [np.zeros(2, dtype="float32"), np.zeros(2, dtype="float32")]

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
            # for body_id, visibility in zip(self._parent.body_id_list, self._parent._visible_F_list):
            for body_id, visibility in zip(self.body_id_list, self._parent._visible_F_list): #  corrected from self._parent.body_id_list, self._parent._visible_F_list to self.body_id_list, self._parent._visible_F_list
                #   create force object
                if direction == "normal":
                    dir = "n"
                elif direction == "tangent":
                    dir = "t"
                else:
                    dir = direction

                F = Force(body_id, force_name=self._parent._name + "_F" + dir + "_on_body_" + str(body_id), visible=visibility, scale=self.scale, parent=self._parent)
                F_list.append(F)

                #   create force visualization datad
                F.set_vtk_data()
                #   add force actor to renderer
                if F.vtk_actor is not None:
                    self._parent._parent._parent._parent.dys.simulation_control_widget.vtkWidget.renderer.AddActor(F.vtk_actor)
                else:
                    print Warning, "Force visualizations not created!"

                #   append force object to list of forces of MBD system object
                self._parent._parent._parent.forces.append(F)

                #   append force object to list of forces of body object
                # self._parent._parent._parent.bodies[body_id].forces.append(F)

        return F_list

    def deactivate_forces(self):
        """

        :return:
        """
        # print "deactivate_forces()@",__name__
        #   join both lists to one
        F_list = self._Fn_list + self._Ft_list

        for i, F in enumerate(F_list):
            # print "i =", i, F._name
            F.deactivate()

            # F._clear_vtk_data()

        self._Fn_list = []
        self._Ft_list = []

    def _remove_forces(self, F_list):
        """

        :param F_list:
        :return:
        """
        if hasattr(self._parent, "body_id_list"):
            for F, body_id in zip(F_list, self._parent.body_id_list):
                #   remove from renderer
                if F.vtk_actor is not None:
                    F.force.Delete()
                    F.vtk_mapper.Delete()

                    self._parent._parent._parent._parent.dys.simulation_control_widget.vtkWidget.renderer.RemoveViewProp(F.vtk_actor)
                    self._parent._parent._parent._parent.dys.simulation_control_widget.vtkWidget.renderer.RemoveActor(F.vtk_actor)

                else:
                    print Warning, "Force visualizations not removed!"

                #   append force object to list of forces of MBD system object
                self._parent._parent._parent.forces.remove(F)

                #   append force object to list of forces of body object
                self._parent._parent._parent.bodies[body_id].forces.remove(F)

                del F

    def set_dq0(self, dq0_n=None, dq0_t=None):
        """
        Function saves the initial normal and tangential contact velocity at impact
        :param dq0_n    a normal contact velocity (scalar, float value)
        :param dq0_t    a tangential contact velocity (scalar, float value)
        """
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


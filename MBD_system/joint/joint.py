"""
Created on 18. mar. 2014

@author: lskrinjar
"""
import time
import itertools
from copy import copy
from pprint import pprint

import numpy as np


from MBD_system.Ai_ui_P import Ai_ui_P_vector
from MBD_system.Ai_theta_ui_P import Ai_theta_ui_P_vector
from MBD_system.Ai_hi import Ai_hi
from MBD_system.Ai_theta_hi import Ai_theta_hi
from MBD_system.hi_Ai_theta_ui_P import hi_Ai_theta_ui_P_constant
from MBD_system.r_ij_P_Ai_theta_hi import r_ij_P_Ai_theta_hi_constant
from MBD_system.r_ij_P import r_ij_P
from MBD_system.u_P_cad2cm_lcs import u_P_cad2cm_lcs
from MBD_system.q2dtheta_i import q2dtheta_i
from MBD_system.q2R_i import q2R_i
from MBD_system.q2dR_i import q2dR_i
from MBD_system.MBD_system_items import JointItem
from MBD_system.q2theta_i import q2theta_i
from MBD_system.force.force import Force
from MBD_system.joint.joint_C_q_matrix import Joint_C_q_matrix
from MBD_system.joint.joint_Q_d_vector import Joint_Q_d_vector
from simulation_control_widget.opengl_widget.marker.marker import Marker


class Joint(JointItem):
    """
    classdocs
    """
    __id = itertools.count(0)

    def __init__(self, joint_type, body_id_i, body_id_j, u_iP_CAD = np.array([0, 0]), u_jP_CAD=np.array([0, 0]), u_iQ_CAD=np.array([0, 0]), parent=None):
        """
        Creates a joint object
        in:
            body_i = body_id
            body_j = body_id
            type (string) - fixed, revolute, prismatic
            u_iP_CAD - in CAD LCS of a body
            u_jP_CAD - in CAD LCS of a body
        """
        super(Joint, self).__init__(joint_type, parent)

        #    number
        self.joint_id = self.__id.next()  # len(joints_list) + 1
        self._parent = parent
        # print "self._parent =", self._parent._children
        # print "self._parent._parent =", self._parent._parent
        # pprint(vars(self._parent._parent))
        # joints_list = self._parent._children
        # print "joints_list =", joints_list
        joints_list = self._parent._parent.joints
        self.joint_type = joint_type

        self.joint_types = ["fixed",
                            "revolute",
                            "prismatic"]

        #   visualization (opengl) properties
        self.z_dim = 0
        #   get q vector of current MBD system at joint object initialization
        self.q0 = self._parent._parent.get_q()

        if self.joint_type == "fixed":
            self.fixed_joint_id = (sum(1 for joint in joints_list if joint.joint_type == "fixed") - 1) + 1
            self._name = self.joint_type + "_" + str(self.fixed_joint_id)

        elif self.joint_type == "revolute":
            self.revolute_joint_id = (sum(1 for joint in joints_list if joint.joint_type == "revolute") - 1) + 1
            self._name = self.joint_type + "_" + str(self.revolute_joint_id)

        elif self.joint_type == "prismatic":
            self.prismatic_joint_id = (sum(1 for joint in joints_list if joint.joint_type == "prismatic") - 1) + 1
            self._name = self.joint_type + "_" + str(self.prismatic_joint_id)

        else:
            raise ValueError, "Joint type not correct!"

        # swap body_id that if body is connected to ground that ground is always the last item in list
        if body_id_i == "ground":
            self.body_id_i = body_id_j
            self.body_id_j = body_id_i
            self.u_iP_CAD = u_jP_CAD
            self.u_jP_CAD = u_iP_CAD

        else:
            self.body_id_i = body_id_i
            self.body_id_j = body_id_j
            self.u_iP_CAD = u_iP_CAD
            self.u_jP_CAD = u_jP_CAD

        #   body id list
        self.body_id_list = [self.body_id_i, self.body_id_j]

        #   list of point vectors to joint constraint in CAD lCS of a body
        self.u_P_CAD_list = [self.u_iP_CAD, self.u_jP_CAD]

        #   predefined empty list to store point vectors of joint constraint in LCS (center of gravity)
        #   of each body
        self.u_P_LCS_list = []

        self.u_iQ_CAD = u_iQ_CAD

        #   predefined variable of vector C(q, t) constraint equation
        self.C0 = None
        self.theta0 = None

        # if self.body_id_i == "ground" or self.body_id_j == "ground":
        #     _body_id = copy(self.body_id_list)
        #     _body_id.remove("ground")
        #     _body_id = _body_id[0]
        #
        #     indx = self.body_id_list.index(_body_id)
        #     C1 = q2R_i(self.q0, _body_id) + Ai_ui_P_vector(self.u_P_CAD_list[indx], q2theta_i(self.q0, _body_id))
        #
        #     C2 = q2theta_i(self.q0, self.body_id_list[indx])
        #
        #     self.C0 = np.append(C1, C2)
        # else:
        #     self.C0 = np.zeros(3)

        #   predefined empty list
        self.C_q_list = []
        self.Q_d_list = []

        #   list of joint points
        #   1   uPi - body i
        #   2   uQi - body i
        #   3   uPj - body j
        self.u_QP_LCS_list = []

        #   pair of contact force list
        self.contact_bodies_added_to_list = False
        self.list_of_contact_force_objects_constructed = False
        self.force_list = []

        #   calculate u_P vector of every body in body LCS
        #   build constraint forces object
        #   build position markers of force
        #   create marker list

        for body_id, _u_P in zip(self.body_id_list, self.u_P_CAD_list):
            if body_id == "ground" or body_id == -1:
                # print "ground"
                u_P_LCS = u_P_cad2cm_lcs(body_id, self._parent._parent.ground, _u_P)
                u_P_LCS = _u_P
            else:

                #   pointer to body object
                _body = self._parent._parent.bodies[body_id]
                #   calculate point vector in body LCS (center of gravity)
                u_P_LCS = u_P_cad2cm_lcs(body_id, _body, _u_P)
                if body_id == 1:
                    print "_u_P(CAD) =", _u_P
                    print "u_P_LCS =", u_P_LCS
            self.u_P_LCS_list.append(u_P_LCS)

        [self.u_iP_LCS, self.u_jP_LCS] = self.u_P_LCS_list

        self.list_of_contact_force_objects_constructed = True
        self.additional_params_calulated = True

        #   create forces
        self._create_Fn_forces()

    def _create_Fn_forces(self):
        """
        Function creates list of contact forces as objects
        :return:
        """
        self.Fn = 0
        self._Fn_list = []
        if self._parent is not None:
            for body_id in self.body_id_list:
                if body_id != "ground":
                    _Fn = Force(body_id, force_name=self._name + "_Fn_on_body_" + str(body_id))
                    #   add pair of contact forces to forces list of MBD system
                    self._parent._parent.forces.append(_Fn)
                    #   add pair of contact forces to forces list of contact
                    self._Fn_list.append(_Fn)
                    self._parent._parent.bodies[body_id].forces.append(_Fn)
            self._Fn_list[0]._visible = False
            self.list_of_contact_force_objects_constructed = True

    def transform_cad2cm_lcs(self, points):
        """
        Function tranforms body points from body CAD CS to body center of mass LCS
        :return:
        """

    def _create_markers(self):
        """
        Function creates markers
        :return:
        """
        markers = []
        for body_id, u_P in zip(self.body_id_list, self.u_P_LCS_list):
            _node = np.array(np.append(u_P, self.z_dim), dtype="float32")

            #   create marker object
            marker = Marker(_node)
            #   append marker object to body markers
            if body_id == "ground":
                self._parent._parent.ground.markers.append(marker)
            else:
                self._parent._parent.bodies[body_id].markers.append(marker)
            #   append marker object to joint markers
            markers.append(marker)

        return markers

    def evaluate_rijP(self, q):
        """

        :param q:
        :return:
        """
        return r_ij_P(q2R_i(q, self.body_id_i), q2theta_i(q, self.body_id_i), self.u_iP_LCS, q2R_i(q, self.body_id_j), q2theta_i(q, self.body_id_j), self.u_jP_LCS)

    def evaluate_C_t(self, q):
        """

        :param q:
        :return:
        """
        C_t = np.zeros(3 - self.joint_DOF)
        return C_t

    def evaluate_d(self, q):
        """

        :return:
        """
        rijP = self.evaluate_rijP(q)
        d = np.linalg.norm(rijP)
        return d

if __name__ == "__main__":
    #     a = Joint(joint_type = "revolute", body_id_i = 1, body_id_j = 2, u_iP_CAD = np.array([1, 2]), u_jP_CAD = np.array([33, 44]))
    a = Joint(joint_type="prismatic", body_id_i=1, body_id_j=2, u_iP_CAD=np.array([1, 2]), u_jP_CAD=np.array([33, 44]),
              u_iQ=np.array([-2, -2]))
    print a
    print a.create_joint_C_q_matrix()
# b = Joint(joint_type = "fixed", body_id_i = 1, body_id_j = 2, u_iP_CAD = np.array([1, 2]), u_jP_CAD = np.array([33, 44]))
#     pprint(vars(b))
#     b.create_joint_C_q_matrix()

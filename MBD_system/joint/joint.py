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

    def __init__(self, joint_type, body_id_i, body_id_j, u_iP_CAD = np.array([0, 0]), u_jP_CAD=np.array([0, 0]),
                 u_iQ_CAD=np.array([0, 0]), parent=None):
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

        #   create marker list
        self.markers = []

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

        self.body_id_list = [self.body_id_i, self.body_id_j]

        #   list of point vectors to joint constraint in CAD lCS of a body
        self.u_P_CAD_list = [self.u_iP_CAD, self.u_jP_CAD]

        #   predefined empty list to store point vectors of joint constraint in LCS (center of gravity)
        #   of each body
        self.u_P_list = []
        #    only for prismatic joint - a normal vector to translational direction
        if self.joint_type == "prismatic":
            self.u_iQ = u_iQ_CAD
        else:
            self.u_iQ = []

        self.z_dim = 0

        #   get q vector of current MBD system at joint object initialization
        self.q0 = self._parent._parent.get_q()

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

        # predefined empty list
        self.C_q_list = []
        self.Q_d_list = []


        #   pair of contact force list
        self.contact_bodies_added_to_list = False
        self.list_of_contact_force_objects_constructed = False
        self.force_list = []

        #   calculate u_P vector of every body in body LCS
        #   build constraint forces object
        #   build position markers of force
        for body_id, _u_P in zip(self.body_id_list, self.u_P_CAD_list):
            if body_id == "ground" or body_id == -1:
                u_P_LCS = u_P_cad2cm_lcs(body_id, self._parent._parent.ground, _u_P=_u_P)
            else:
                #   create pointer to body
                _body = self._parent._parent.bodies[body_id]

                #   calculate point vector in body LCS (center of gravity)
                u_P_LCS = u_P_cad2cm_lcs(body_id, _body, _u_P=_u_P)

                #   create force object
                _force = Force(body_id, force_name=self._name + "_on_body_" + str(body_id))
                #   add pair of contact forces to forces list of MBD system
                self._parent._parent.forces.append(_force)

                #   add pair of contact forces to forces list of contact
                self.force_list.append(_force)
                self._parent._parent.bodies[body_id].forces.append(_force)

                #   create markers for u_P points on every body in kinematic constraint
                _node = np.array(np.append(u_P_LCS, self.z_dim), dtype='float32')
                _marker = Marker(_node, visible=True)
                self._parent._parent.bodies[body_id].markers.append(_marker)
                self.markers.append(_marker)

            self.u_P_list.append(u_P_LCS)

        [self.u_iP, self.u_jP] = self.u_P_list
        # print "self.u_iP, self.u_jP =", self.u_iP, self.u_jP
        self.list_of_contact_force_objects_constructed = True
        self.additional_params_calulated = True

    def evaluate_C_q(self, q):
        """
        Function creates C_q matrix for the specified joint type.
        Args:
        bodies - list of all bodies in MBD system
        q - vector of absolute coordinates (positions and velocities)
        Returns:
        C_q matrix for each body in joint
        Raises:
        """

        if self.joint_type == "fixed":
            return self._create_fixed_joint_C_q(q, bodies)
        elif self.joint_type == "revolute":
            return self._create_revolute_joint_C_q(q)
        elif self.joint_type == "prismatic":
            return self._create_prismatic_joint_C_q(q)
        else:
            raise IOError, "joint_type not correct. Please specify correct joint type."

    def _create_fixed_joint_C_q(self, q_, bodies):  # self, bodies
        """
        Create fixed joint matrix for each body connected with the joint.
        """
        self.joint_DOF = 0
        #    calculate additional parameters only once
        #    assigns properties of body object to variable (object) from list of bodies and body index
        #    calculates vector of point of joint on each body from cad lcs to cm lcs of a body
        # if not self.additional_params_calulated:
        #     self.u_iP_cm = u_P_cad2cm_lcs(_body_id=self.body_id_i, _u_P=self.u_iP_CAD, _bodies=bodies)
        #     self.u_jP_cm = u_P_cad2cm_lcs(_body_id=self.body_id_j, _u_P=self.u_jP_CAD, _bodies=bodies)
        #
        #     self.u_P_list = [self.u_iP_cm, self.u_jP_cm]
        #
        #     self.additional_params_calulated = True


        self.C_q_list = []
        #    construct C_q matrix for fixed joint
        for body_id in self.body_id_list:
            if body_id == "ground":
                joint_C_q_matrix_body = Joint_C_q_matrix(joint_type_=self.joint_type)
            else:
                joint_C_q_matrix_body = Joint_C_q_matrix(joint_type_=self.joint_type)

                joint_C_q_matrix_body.matrix[[0, 1], -1] = Ai_theta_ui_P_vector(
                    u_P=self.u_P_list[joint_C_q_matrix_body.id],
                    theta=q_[3 * self.body_id_list[joint_C_q_matrix_body.id] + 2], _id=joint_C_q_matrix_body.id)

            self.C_q_list.append(joint_C_q_matrix_body)

        [i_joint_C_q_matrix, j_joint_C_q_matrix] = self.C_q_list

        return i_joint_C_q_matrix, j_joint_C_q_matrix

    def _create_revolute_joint_C_q(self, q_):
        """
        Create revolute joint C_q matrix for each body connected with the joint.
        """
        self.joint_DOF = 1
        #    predefine empty list to store joint C_q matrix for each body in joint
        self.C_q_list = []
        #    construct C_q matrix for revolute joint
        for body_id, _u_P, _id in zip(self.body_id_list, self.u_P_list, [0, 1]):

            if body_id == "ground":
                joint_C_q_matrix_body = Joint_C_q_matrix(self.joint_type)
            else:
                joint_C_q_matrix_body = Joint_C_q_matrix(self.joint_type)
                #   get body theta from q
                _theta = q2theta_i(q_, body_id)
                #   create body Cq matrix
                joint_C_q_matrix_body.matrix[:, -1] = Ai_theta_ui_P_vector(_u_P, _theta, _id)

            # append joint C_q matrix to list
            self.C_q_list.append(joint_C_q_matrix_body)

        [i_joint_C_q_matrix, j_joint_C_q_matrix] = self.C_q_list

        return i_joint_C_q_matrix, j_joint_C_q_matrix

    def _create_prismatic_joint_C_q(self, q_):
        """
        Create prismatic joint matrix for each body connected with the joint.
        """
        self.joint_DOF = 1
        #    calculate additional parameters only once
        #    assigns properties of body object to variable (object) from list of bodies and body index
        #    calculates vector of point of joint on each body from cad lcs to cm lcs of a body
        if not self.additional_params_calulated:
            self.u_iP_cm = u_P_cad2cm_lcs(_body_id=self.body_id_i, _u_P=self.u_iP_CAD)
            self.u_jP_cm = u_P_cad2cm_lcs(_body_id=self.body_id_j, _u_P=self.u_jP_CAD)

            self.u_iQ_cm = u_P_cad2cm_lcs(_body_id=self.body_id_i, _u_P=self.u_iQ)

            self.h_i_cm = self.u_iP_cm - self.u_iQ_cm

            self.u_P_list = [self.u_iP_cm, self.u_jP_cm]

            self.additional_params_calulated = True

        # calculate h_i
        self.h_i = Ai_hi(h=self.h_i_cm, theta=q2theta_i(q_, self.body_id_list[0]))
        hi_x, hi_y = self.h_i

        #    predefine empty list to store joint C_q matrix for each body in joint
        self.C_q_list = []
        #    construct C_q matrix for prismatic joint
        for body_id in self.body_id_list:

            if body_id == "ground":
                joint_C_q_matrix_body = Joint_C_q_matrix(joint_type_=self.joint_type)
            else:
                joint_C_q_matrix_body = Joint_C_q_matrix(joint_type_=self.joint_type)

                joint_C_q_matrix_body.matrix[0, 0] = hi_x
                joint_C_q_matrix_body.matrix[0, 1] = hi_y
                joint_C_q_matrix_body.matrix[0, 2] = hi_Ai_theta_ui_P_constant(self.h_i,
                                                                               q2theta_i(q_, self.body_id_list[0]),
                                                                               self.u_P_list[
                                                                                   self.body_id_list.index(body_id)])
                joint_C_q_matrix_body.matrix[1, 2] = 1


                #    body i matrix
                if joint_C_q_matrix_body.id == 0:
                    self.rijP_ = r_ij_P(R_i=q2R_i(q_, self.body_id_i), theta_i=q2theta_i(q_, self.body_id_i),
                                        u_iP_CAD=self.u_iP_cm, R_j=q2R_i(q_, self.body_id_j),
                                        theta_j=q2theta_i(q_, self.body_id_j), u_jP_CAD=self.u_jP_cm)
                    joint_C_q_matrix_body.matrix[0, 2] = joint_C_q_matrix_body.matrix[
                                                             0, 2] + r_ij_P_Ai_theta_hi_constant(rijP=self.rijP_,
                                                                                                 theta_i=q2theta_i(q_,
                                                                                                                   self.body_id_i),
                                                                                                 hi=self.h_i_cm)

                # body j matrix
                elif joint_C_q_matrix_body.id == 1:
                    joint_C_q_matrix_body.matrix = -joint_C_q_matrix_body.matrix



            # append joint C_q matrix to list
            self.C_q_list.append(joint_C_q_matrix_body)

        [i_joint_C_q_matrix, j_joint_C_q_matrix] = self.C_q_list

        return i_joint_C_q_matrix, j_joint_C_q_matrix

    def evaluate_Q_d(self, q):
        """
        Function creates Q_d vector for joint if joint type is fixed or revolute.
        Args:
            bodies - list of all bodies
            q - array of ppositions and velosities
            N - number of bodies
        Returns:
            Q_d - vector Q_d for joint
        """

        if self.joint_type == "fixed":
            return self._create_fixed_joint_Q_d(bodies_=bodies, q_=q, N_b_=N_b)
        elif self.joint_type == "revolute":
            return self._create_revolute_joint_Q_d(q)
        elif self.joint_type == "prismatic":
            return self._create_prismatic_joint_Q_d(bodies_=bodies, q_=q, N_b_=N_b)
        else:
            raise IOError, "joint_type not correct. Please specify right joint type."

    def _create_fixed_joint_Q_d(self, bodies_, q_, N_b_):
        """
        Q_d equation according to eq. 3.133 in Computational Dynamics 3rd Ed. (A.A.Shabana)
        expanded with z component of moment
        """
        #    if body i and j are not ground
        if self.body_id_i != "ground" and self.body_id_j != "ground":
            self.Q_d_list = []

            for body_id in self.body_id_list:
                joint_Q_d_vector_body = Joint_Q_d_vector(joint_type_=self.joint_type, body_connected_to_ground=False)

                Ai_ui_P_vector_ = Ai_ui_P_vector(u_P=self.u_P_list[joint_Q_d_vector_body.id],
                                                 theta=q2theta_i(q_, body_id))

                dtheta_body = q2dtheta_i(q_, body_id, N_b_)

                joint_Q_d_vector_body.joint_Q_d_vector[0:2] = Ai_ui_P_vector_ * (dtheta_body ** 2)

                #    add calculated joint matrix of a body to list
                self.Q_d_list.append(joint_Q_d_vector_body)

            Q_d_vector_ = +(
            self.Q_d_list[0].joint_Q_d_vector - self.Q_d_list[1].joint_Q_d_vector)


        elif self.body_id_i == "ground":
            joint_Q_d_vector_body = Joint_Q_d_vector(joint_type_=self.joint_type, body_connected_to_ground=True)

            joint_Q_d_vector_body.id = 1

            body_id_ = self.body_id_list[joint_Q_d_vector_body.id]
            #    calculate Q_d vector of only non-ground body
            Ai_ui_P_vector_ = Ai_ui_P_vector(u_P=self.u_P_list[joint_Q_d_vector_body.id], theta=q2theta_i(q_, body_id_))

            dtheta_body = q2dtheta_i(q_, body_id_, N_b_)

            joint_Q_d_vector_body.joint_Q_d_vector[0:2] = Ai_ui_P_vector_ * (dtheta_body ** 2)

            Q_d_vector_ = joint_Q_d_vector_body.joint_Q_d_vector


        elif self.body_id_j == "ground":
            joint_Q_d_vector_body = Joint_Q_d_vector(joint_type_=self.joint_type, body_connected_to_ground=True)

            body_id_ = self.body_id_list[joint_Q_d_vector_body.id]

            #    calculate Q_d vector of only non-ground body
            Ai_ui_P_vector_ = Ai_ui_P_vector(u_P=self.u_P_list[joint_Q_d_vector_body.id], theta=q2theta_i(q_, body_id_))

            dtheta_body = q2dtheta_i(q_, body_id_)

            joint_Q_d_vector_body.joint_Q_d_vector[0:2] = Ai_ui_P_vector_ * dtheta_body ** 2

            Q_d_vector_ = joint_Q_d_vector_body.joint_Q_d_vector
        else:
            raise ValueError, "Q_d vector for fixed joint not constructed. Check calculation process."

        return Q_d_vector_

    def _create_revolute_joint_Q_d(self, q):
        """
        Q_d equation according to eq. 3.133 in Computational Dynamics 3rd Ed. (A.A.Shabana)
        """
        #    if body i and body j are not ground
        if self.body_id_i != "ground" and self.body_id_j != "ground":

            self.Q_d_list = []
            for body_id in self.body_id_list:
                joint_Q_d_vector_body = Joint_Q_d_vector(joint_type_=self.joint_type)

                Ai_ui_P_vector_ = Ai_ui_P_vector(u_P=self.u_P_list[joint_Q_d_vector_body.id],
                                                 theta=q2theta_i(q, body_id))

                dtheta_body = q2dtheta_i(q, body_id)

                joint_Q_d_vector_body.joint_Q_d_vector = Ai_ui_P_vector_ * (dtheta_body ** 2)

                #    add calculated joint matrix of a body to list
                self.Q_d_list.append(joint_Q_d_vector_body)

            Q_d_vector_ = self.Q_d_list[0].joint_Q_d_vector - self.Q_d_list[1].joint_Q_d_vector


        # construct Q_d vector if body i is ground
        elif self.body_id_i == "ground":
            joint_Q_d_vector_body = Joint_Q_d_vector(joint_type_=self.joint_type, body_connected_to_ground=True)

            joint_Q_d_vector_body.id = 1

            _body_id = self.body_id_list[joint_Q_d_vector_body.id]

            # _theta = q2theta_i(q_, _body_id)
            #    calculate Q_d vector of only non-ground body
            Ai_ui_P_vector_ = Ai_ui_P_vector(u_P=self.u_P_list[joint_Q_d_vector_body.id], _theta=q2theta_i(q, _body_id))

            dtheta_body = q2dtheta_i(q, _body_id)

            joint_Q_d_vector_body.joint_Q_d_vector = Ai_ui_P_vector_ * (dtheta_body ** 2)

            Q_d_vector_ = joint_Q_d_vector_body.joint_Q_d_vector

        # construct Q_d vector if body j is ground
        elif self.body_id_j == "ground":

            joint_Q_d_vector_body = Joint_Q_d_vector(joint_type_=self.joint_type, body_connected_to_ground = True)

            _body_id = self.body_id_list[joint_Q_d_vector_body.id]

            #   get evaluated vector to point of kinematic constraint u_P
            _u_P = self.u_P_list[joint_Q_d_vector_body.id]
            #   theta of body from vector q
            _theta = q2theta_i(q, _body_id)
            #    calculate Q_d vector of only non-ground body
            Ai_ui_P_vector_ = Ai_ui_P_vector(_u_P, _theta)
            #   omega of body from vector q
            dtheta_body = q2dtheta_i(q, _body_id)

            joint_Q_d_vector_body.joint_Q_d_vector = Ai_ui_P_vector_ * dtheta_body ** 2

            Q_d_vector_ = joint_Q_d_vector_body.joint_Q_d_vector

        else:
            raise ValueError, "Q_d vector for revolute joint not constructed. Check calculation process."

        return Q_d_vector_

    def _create_prismatic_joint_Q_d(self, bodies_, q_, N_b_):
        """
        Q_d equation according to eq. 3.133 (example 3.14) in Computational Dynamics 3rd Ed. (A.A.Shabana)
        """
        Q_d_vector_ = np.zeros(2)

        theta_i = q2theta_i(q_, self.body_id_i)
        theta_j = q2theta_i(q_, self.body_id_j)

        R_i = q2R_i(q_, self.body_id_i)

        dR_i = q2dR_i(q_, self.body_id_i, N_b_)
        dtheta_i = q2dtheta_i(q_, self.body_id_i, N_b_)

        dR_j = q2dR_i(q_, self.body_id_j, N_b_)
        dtheta_j = q2dtheta_i(q_, self.body_id_j, N_b_)

        Q_d_21 = 2 * dtheta_i * np.dot(Ai_theta_hi(self.h_i_cm, theta_i), (dR_i - dR_j))
        Q_d_22 = (dtheta_i ** 2) * (np.dot(self.u_iP_cm, self.h_i_cm) - np.dot(self.rijP_, self.h_i))
        Q_d_23 = 2 * dtheta_i * dtheta_j * np.dot(Ai_theta_hi(self.h_i_cm, theta_i),
                                                  Ai_ui_P_vector(self.u_jP_cm, theta_j))
        Q_d_24 = (dtheta_j ** 2) * np.dot(self.h_i, R_i + Ai_ui_P_vector(self.u_jP_cm, theta_j))

        Q_d_vector_[0] = -Q_d_21 - Q_d_22 + Q_d_23 - Q_d_24

        return Q_d_vector_

    def evaluate_C(self, q_):
        """
        Function creates C vector for the specified joint type.
        Args:
        bodies - list of all bodies in MBD system
        q - vector of absolute coordinates (positions and velocities)
        Returns:
        C_q matrix for each body in joint
        Raises:
        """
        if self.joint_type == "fixed":
            return self._create_C_fixed_joint(q_)
        elif self.joint_type == "revolute":
            return self._create_C_revolute_joint(q_)
        elif self.joint_type == "prismatic":
            return self._create_C_prismatic_joint(q_)
        else:
            raise IOError, "joint_type not correct. Please specify correct joint type."

    def _create_C_fixed_joint(self, q_):
        """
        
        """
        body_id_i = self.body_id_list[0]
        body_id_j = self.body_id_list[1]

        #   if body not connected to ground
        if self.body_id_i != "ground" and self.body_id_j != "ground":
            C1 = self._create_C_revolute_joint(q_)
            C2 = q2theta_i(q_, body_id_i) - q2theta_i(q_, body_id_j)
            C = np.append(C1, C2)
        # body connected to ground
        else:
            __body_id = copy(self.body_id_list)
            __body_id.remove("ground")
            __body_id = __body_id[0]

            indx = self.body_id_list.index(__body_id)
            C1 = q2R_i(q_, self.body_id_list[indx])
            C2 = q2theta_i(q_, self.body_id_list[indx])

            C = np.append(C1, C2) - self.C
        return C

    def _create_C_revolute_joint(self, q_):
        """
        
        """
        body_id_i = self.body_id_list[0]
        body_id_j = self.body_id_list[1]

        #   revolute joint
        if self.body_id_i != "ground" and self.body_id_j != "ground":
            C = q2R_i(q_, body_id_i) + Ai_ui_P_vector(self.u_P_CAD_list[0], q2theta_i(q_, body_id_i)) - q2R_i(q_,
                                                                                                              body_id_j) - Ai_ui_P_vector(
                self.u_P_CAD_list[1], q2theta_i(q_, body_id_j))
        # revolute joint to ground
        else:
            #   check which body is and which body is not ground
            __body_id = copy(self.body_id_list)
            __body_id.remove("ground")
            __body_id = __body_id[0]

            indx = self.body_id_list.index(__body_id)

            C = q2R_i(q_, __body_id) + Ai_ui_P_vector(self.u_P_CAD_list[indx], q2theta_i(q_, __body_id)) - self.C[0:2]

        return C

    def _create_C_prismatic_joint(self):
        """
        
        """

    def evaluate_rijP(self, q):
        """

        :param q:
        :return:
        """
        return r_ij_P(q2R_i(q, self.body_id_i), q2theta_i(q, self.body_id_i), self.u_iP, q2R_i(q, self.body_id_j), q2theta_i(q, self.body_id_j), self.u_jP)

if __name__ == "__main__":
    #     a = Joint(joint_type = "revolute", body_id_i = 1, body_id_j = 2, u_iP_CAD = np.array([1, 2]), u_jP_CAD = np.array([33, 44]))
    a = Joint(joint_type="prismatic", body_id_i=1, body_id_j=2, u_iP_CAD=np.array([1, 2]), u_jP_CAD=np.array([33, 44]),
              u_iQ=np.array([-2, -2]))
    print a
    print a.create_joint_C_q_matrix()
# b = Joint(joint_type = "fixed", body_id_i = 1, body_id_j = 2, u_iP_CAD = np.array([1, 2]), u_jP_CAD = np.array([33, 44]))
#     pprint(vars(b))
#     b.create_joint_C_q_matrix()

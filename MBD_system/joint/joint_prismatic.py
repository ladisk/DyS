"""

created by: lskrinjar
date of creation: 01/02/2016
time of creation: 10:11
"""

import numpy as np

from MBD_system.A import A_matrix
from MBD_system.A_theta import A_theta_matrix
from MBD_system.joint.joint import Joint
from MBD_system.joint.joint_C_q_matrix import Joint_C_q_matrix
from MBD_system.q2dR_i import q2dR_i
from MBD_system.q2dtheta_i import q2dtheta_i
from MBD_system.q2theta_i import q2theta_i
from MBD_system.transform_cs import u_P_cad2cm_lcs
from simulation_control_widget.opengl_widget.marker.marker import Marker


class JointPrismatic(Joint):
    """
    classdocs
    """
    def __init__(self, joint_type, body_id_i, body_id_j, u_iP_CAD=np.array([0, 0]), u_jP_CAD=np.array([0, 0]), u_iQ_CAD=np.array([0, 0]), parent=None):
        """

        :param joint_type:
        :param body_id_i:
        :param body_id_j:
        :param u_iP_CAD:
        :param u_jP_CAD:
        :param u_iQ:
        :param parent:
        :return:
        """
        self.joint_type = "prismatic"
        super(JointPrismatic, self).__init__(joint_type, body_id_i, body_id_j, u_iP_CAD=u_iP_CAD, u_jP_CAD=u_jP_CAD, u_iQ_CAD=u_iQ_CAD, parent=parent)

        #   joint DOF
        self.joint_DOF = 1

        #   evaluate initial angle between connected bodies
        self.theta0 = self.evaluate_C0(self.q0)

        #   add point u_iQ to an extended list
        self.u_iQ_LCS = u_P_cad2cm_lcs(self.body_id_i, self._parent._parent.bodies[self.body_id_i], self.u_iQ_CAD)
        #   this extended list has elements in form of:
        #   [uiP, ujP, uiQ]
        self.u_QP_LCS_list = self.u_P_LCS_list
        self.u_QP_LCS_list.append(self.u_iQ_LCS)

        #   vector h_i
        self.h_i_LCS = self.u_iP_LCS - self.u_iQ_LCS

        #   markers
        self.markers = self._create_markers()

    def _create_markers(self):
        """
        Function creates markers
        :return:
        """
        markers = []

        for i, u_P in enumerate(self.u_QP_LCS_list):
            _node = np.array(np.append(u_P, self.z_dim), dtype="float32")
            if i == 0 or i == 2:
                body_id = self.body_id_i
            if i == 1:
                body_id = self.body_id_j
            #   create marker object
            marker = Marker(_node)
            #   append marker object to body markers
            if isinstance(body_id, int):
                self._parent._parent.bodies[body_id].markers.append(marker)
            else:
                self._parent._parent.ground.markers.append(marker)
            #   append marker object to joint markers
            markers.append(marker)

        return markers

    def evaluate_C_q(self, q):
        """
        Function evaluates C_q matrix of prismatic joint
        :return:
        """
        #    calculate additional parameters only once
        #    assigns properties of body object to variable (object) from list of bodies and body index
        #    calculates vector of point of joint on each body from cad lcs to cm lcs of a body
        # if not self.additional_params_calulated:
        #     self.u_iP_cm = u_P_cad2cm_lcs(_body_id=self.body_id_i, _u_P=self.u_iP_CAD)
        #     self.u_jP_cm = u_P_cad2cm_lcs(_body_id=self.body_id_j, _u_P=self.u_jP_CAD)
        #
        #     self.u_iQ_cm = u_P_cad2cm_lcs(_body_id=self.body_id_i, _u_P=self.u_iQ)
        #
        #     self.h_i_cm = self.u_iP_cm - self.u_iQ_cm
        #
        #     self.u_P_list = [self.u_iP_cm, self.u_jP_cm]
        #
        #     self.additional_params_calulated = True

        #   evaluate h_i
        self.h_i = self.evaluate_hi(q)
        hi_x, hi_y = self.h_i

        self.rijP = self.evaluate_rijP(q)

        #    predefine empty list to store joint C_q matrix for each body in joint
        self.C_q_list = []
        #    construct C_q matrix for prismatic joint
        for body_id in self.body_id_list:
            if body_id == "ground":
                C_q_body = Joint_C_q_matrix(joint_type_=self.joint_type)
            else:
                C_q_body = Joint_C_q_matrix(joint_type_=self.joint_type)

                C_q_body.matrix[0, 0] = hi_x
                C_q_body.matrix[0, 1] = hi_y
                C_q_body.matrix[1, 2] = 1

                #    body i matrix
                if body_id == self.body_id_i:
                    self.rijP = self.evaluate_rijP(q)
                    C_q_body.matrix[0, 2] = np.dot(self.rijP, np.dot(A_theta_matrix(q2theta_i(q, body_id)), self.h_i_LCS)) + np.dot(self.h_i, np.dot(A_theta_matrix(q2theta_i(q, body_id)), self.u_iP_LCS))
                    #r_ij_P_Ai_theta_hi_constant(self.rijP, q2theta_i(q, self.body_id_i), self.h_i_LCS)

                # body j matrix
                if body_id == self.body_id_j:
                    C_q_body.matrix = -C_q_body.matrix
                    C_q_body.matrix[0, 2] = - np.dot(self.h_i, np.dot(A_theta_matrix(q2theta_i(q, body_id)), self.u_jP_LCS))
                    #hi_Ai_theta_ui_P_constant(self.h_i, q2theta_i(q, self.body_id_list[0]), self.u_P_LCS_list[self.body_id_list.index(body_id)])

            #   append joint C_q matrix to list
            self.C_q_list.append(C_q_body)
        [C_qi, C_qj] = self.C_q_list
        return C_qi.matrix, C_qj.matrix

    def evaluate_Q_d(self, q):
        """
        Function evaluates Q_d vector of a joint
        Q_d equation according to eq. 3.133 (example 3.14) in Computational Dynamics 3rd Ed. (A.A.Shabana)
        Args:
            bodies - list of all bodies
            q - array of ppositions and velosities
            N - number of bodies
        Returns:
            Q_d - vector Q_d for joint
        """
        Q_d = np.zeros(2)

        #   get theta i, j
        theta_i = q2theta_i(q, self.body_id_i)
        theta_j = q2theta_i(q, self.body_id_j)

        #   get translational absolute velocity of body - R i, j
        dR_i = q2dR_i(q, self.body_id_i)
        dR_j = q2dR_i(q, self.body_id_j)

        #   get rotational absolute velocity of body i, j
        dtheta_i = q2dtheta_i(q, self.body_id_i)
        dtheta_j = q2dtheta_i(q, self.body_id_j)

        #   evaluate h_i
        hi = self.evaluate_hi(q)
        #   evaluate distance between
        rijP = self.evaluate_rijP(q)

        Q_d_21 = 2 * theta_i * np.dot(np.dot(A_theta_matrix(theta_i), hi), (dR_i - dR_j))
        Q_d_22 = (dtheta_i ** 2) * (np.dot(self.u_iP_LCS, self.u_iQ_LCS) - np.dot(rijP, hi))
        Q_d_23 = 2 * dtheta_i * dtheta_j * np.dot(np.dot(A_theta_matrix(theta_i), hi), np.dot(A_theta_matrix(theta_j), self.u_jP_LCS))
        Q_d_24 = (dtheta_j ** 2) * np.dot(hi, self.u_jP_LCS)
        Q_d[0] = -Q_d_21 - Q_d_22 + Q_d_23 - Q_d_24

        return Q_d

    def evaluate_C(self, q, t=None):
        """
        Function creates C vector for the joint
        Args:
        bodies - list of all bodies in MBD system
        q - vector of absolute coordinates (positions and velocities)
        Returns:
        C_q matrix for each body in joint
        Raises:
        """
        C = np.zeros(2)

        #   dot product of orthogonal vectors hi and rijP is 0
        rijP = self.evaluate_rijP(q)
        hi = self.evaluate_hi(q)
        C[0] = np.dot(hi, rijP)

        #   angle between connected bodies is constant
        C[1] = q2theta_i(q, self.body_id_i) - q2theta_i(q, self.body_id_j) - self.theta0
        return C

    def evaluate_hi(self, q):
        """

        :param q:
        :return:
        """
        h_i = np.dot(A_matrix(q2theta_i(q, self.body_id_i)), self.h_i_LCS)
        return h_i

    def evaluate_C0(self, q):
        """

        :param q:
        :return:
        """
        return q2theta_i(q, self.body_id_i) - q2theta_i(q, self.body_id_j)



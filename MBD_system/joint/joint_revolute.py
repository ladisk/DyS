"""

created by: lskrinjar
date of creation: 01/02/2016
time of creation: 10:11
"""

import numpy as np

from MBD_system.Ai_theta_ui_P import Ai_theta_ui_P_vector
from MBD_system.Ai_ui_P import Ai_ui_P_vector
from MBD_system.joint.joint import Joint
from MBD_system.joint.joint_C_q_matrix import Joint_C_q_matrix
from MBD_system.joint.joint_Q_d_vector import Joint_Q_d_vector
from MBD_system.q2dtheta_i import q2dtheta_i
from MBD_system.q2theta_i import q2theta_i


class JointRevolute(Joint):
    """
    Class of a revolute joint between two rigid bodies
    """

    def __init__(self, joint_type, body_id_i, body_id_j, u_iP_CAD=np.array([0, 0], dtype=float), u_jP_CAD=np.array([0, 0], dtype=float), u_iQ=np.array([0, 0], dtype=float), properties_dict={}, parent=None):
        """

        """
        super(JointRevolute, self).__init__(joint_type, body_id_i, body_id_j, u_iP_CAD=u_iP_CAD, u_jP_CAD=u_jP_CAD, properties_dict=properties_dict, parent=parent)

        #   auto generated name by type and consecutive number
        self.revolute_joint_id = (sum(1 for joint in self.joints_list if joint.joint_type == "revolute") - 1) + 1
        self._name = self.joint_type + "_" + str(self.revolute_joint_id)

        #   joint DOF
        self.joint_DOF = 1

        #   number of constrained (nodal) coordinates
        self.n_CNC = 3 - self.joint_DOF

        #   C_q matrix dimensions [rows, cols]
        self.C_q_dim = [self.n_CNC, [3, 3]]

        #   initial evaluation of constraint equation for revolute joint
        self.C0 = self.evaluate_C(self.q0)

        #   list of markers
        self.markers = self._create_markers()

    def evaluate_C_q(self, q):
        """
        Function evaluates C_q matrix
        :return:
        """
        #    predefine empty list to store joint C_q matrix for each body in joint
        self.C_q_list = []
        #    construct C_q matrix for each body in joint
        for body_id, _u_P, _id in zip(self.body_id_list, self.u_P_LCS_list, [0, 1]):
            if body_id == "ground":
                C_q_body = Joint_C_q_matrix(self.joint_type)
            else:
                C_q_body = Joint_C_q_matrix(self.joint_type)
                #   create body Cq matrix
                C_q_body.matrix[:, -1] = Ai_theta_ui_P_vector(_u_P, q2theta_i(q, body_id), _id)

            # append joint C_q matrix to list
            self.C_q_list.append(C_q_body)

        [C_qi, C_qj] = self.C_q_list

        return C_qi.matrix, C_qj.matrix

    def evaluate_Q_d(self, q):
        """
        Function evaluates Q_d vector of a joint
        Q_d equation according to eq. 3.133 in Computational Dynamics 3rd Ed. (A.A.Shabana)
        Args:
            bodies - list of all bodies
            q - array of ppositions and velosities
            N - number of bodies
        Returns:
            Q_d - vector Q_d for joint
        """
        #    if body i and body j are not ground
        if (self.body_id_i != "ground") and (self.body_id_j != "ground"):
            self.Q_d_list = []
            for body_id in self.body_id_list:
                #   vector Q_d object for each body
                Q_d_body = Joint_Q_d_vector(self.joint_type)
                #   evaluated for each body
                Q_d_body.Q_d = Ai_ui_P_vector(self.u_P_LCS_list[Q_d_body.id], q2theta_i(q, body_id)) * (q2dtheta_i(q, body_id) ** 2)
                #    add calculated joint matrix of a body to list
                self.Q_d_list.append(Q_d_body)

            Q_d = self.Q_d_list[0].Q_d - self.Q_d_list[1].Q_d

        #   construct Q_d vector if body i is ground
        elif self.body_id_i == "ground":
            Q_d_body = Joint_Q_d_vector(self.joint_type, body_connected_to_ground=True)
            Q_d_body.Q_d = Ai_ui_P_vector(self.u_P_LCS_list[Q_d_body.id], q2theta_i(q, self.body_id_j)) * (q2dtheta_i(q, self.body_id_j) ** 2)

            Q_d = Q_d_body.Q_d

        #   construct Q_d vector if body j is ground
        elif self.body_id_j == "ground":

            Q_d_body = Joint_Q_d_vector(self.joint_type, body_connected_to_ground=True)
            Q_d_body.Q_d = Ai_ui_P_vector(self.u_P_LCS_list[Q_d_body.id], q2theta_i(q, self.body_id_i)) * (q2dtheta_i(q, self.body_id_i) ** 2)
            Q_d = Q_d_body.Q_d

        else:
            raise ValueError, "Q_d vector for revolute joint not constructed. Check calculation process."

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
        if self.body_id_i != "ground" or self.body_id_j != "ground":
            self.C = self.evaluate_rijP(q)
            return self.C
        elif self.body_id_i == "ground":
            self.C = Ai_ui_P_vector(self.u_jP, q2theta_i(q, self.body_id_j))
            if self.C0 is not None:
                self.C = self.C - self.C0
            return self.C
        elif self.body_id_j == "ground":
            self.C = Ai_ui_P_vector(self.u_jP, q2theta_i(q, self.body_id_j))
            if self.C0 is not None:
                self.C = self.C - self.C0
            return self.C
        else:
            raise ValueError, "Not right!"




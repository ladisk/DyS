"""

created by: lskrinjar
date of creation: 01/02/2016
time of creation: 10:11
"""


import numpy as np


from MBD_system.joint.joint import Joint
from MBD_system.Ai_ui_P import Ai_ui_P_vector
from MBD_system.Ai_theta_ui_P import Ai_theta_ui_P_vector
from MBD_system.joint.joint_C_q_matrix import Joint_C_q_matrix
from MBD_system.joint.joint_Q_d_vector import Joint_Q_d_vector
from MBD_system.q2R_i import q2R_i
from MBD_system.q2dR_i import q2dR_i
from MBD_system.q2theta_i import q2theta_i
from MBD_system.q2dtheta_i import q2dtheta_i


class JointFixed(Joint):
    """
    classdocs
    """
    def __init__(self, joint_type, body_id_i, body_id_j, u_iP_CAD=np.array([0, 0]), u_jP_CAD=np.array([0, 0]), u_iQ=np.array([0, 0]), parent=None):
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
        self.joint_type = "fixed"
        super(JointFixed, self).__init__(self.joint_type, body_id_i, body_id_j, u_iP_CAD=u_iP_CAD, u_jP_CAD=u_jP_CAD, parent=parent)

        #   joint DOF
        self.joint_DOF = 0

        #   initial evaluation of constraint equation for fixed joint
        self.C0 = self.evaluate_C0(self.q0)
        self.theta0 = self.C0[2]

        #   list of markers
        self.markers = self._create_markers()

    def evaluate_C_q(self, q):
        """
        Function evaluates C_q matrix
        :return:
        """
        self.C_q_list = []
        #    construct C_q matrix for fixed joint
        for body_id, _u_P in zip(self.body_id_list, self.u_P_LCS_list):
            if body_id == "ground":
                joint_C_q_matrix_body = Joint_C_q_matrix(joint_type_=self.joint_type)
            else:
                joint_C_q_matrix_body = Joint_C_q_matrix(joint_type_=self.joint_type)

                joint_C_q_matrix_body.matrix[[0, 1], -1] = Ai_theta_ui_P_vector(_u_P, q2theta_i(q, body_id), joint_C_q_matrix_body.id)

            self.C_q_list.append(joint_C_q_matrix_body)

        [C_qi, C_qj] = self.C_q_list

        return C_qi.matrix, C_qj.matrix

    def evaluate_Q_d(self, q):
        """
        Q_d equation according to eq. 3.133 in Computational Dynamics 3rd Ed. (A.A.Shabana)
        expanded with z component of moment
        """
        #    if body i and j are not ground
        if self.body_id_i != "ground" and self.body_id_j != "ground":
            self.Q_d_list = []
            for body_id in self.body_id_list:
                Q_d_body = Joint_Q_d_vector(self.joint_type)

                Q_d_body.Q_d[0:2] = Ai_ui_P_vector(self.u_P_list[Q_d_body.id], q2theta_i(q, body_id)) * (q2dtheta_i(q, body_id) ** 2)

                #    add calculated joint matrix of a body to list
                self.Q_d_list.append(Q_d_body)

            Q_d = +(self.Q_d_list[0].Q_d - self.Q_d_list[1].Q_d)


        elif self.body_id_i == "ground":
            Q_d_body = Joint_Q_d_vector(self.joint_type, body_connected_to_ground=True)

            Q_d_body.id = 1

            body_id_ = self.body_id_list[Q_d_body.id]
            #    calculate Q_d vector of only non-ground body
            Q_d_body.Q_d[0:2] = Ai_ui_P_vector(u_P=self.u_P_list[Q_d_body.id], theta=q2theta_i(q, body_id_)) * (q2dtheta_i(q, body_id_) ** 2)

            Q_d = Q_d_body.Q_d


        elif self.body_id_j == "ground":
            Q_d_body = Joint_Q_d_vector(self.joint_type, body_connected_to_ground=True)
            body_id_ = self.body_id_list[Q_d_body.id]

            #    calculate Q_d vector of only non-ground body
            Q_d_body.Q_d[0:2] = Ai_ui_P_vector(self.u_P_LCS_list[Q_d_body.id], q2theta_i(q, body_id_)) * (q2dtheta_i(q, body_id_) ** 2)

            Q_d = Q_d_body.Q_d
        else:
            raise ValueError, "Q_d vector for fixed joint not constructed. Check calculation process."

        return Q_d

    def evaluate_C0(self, q):
        """
        Function evaluates C vector for the joint
        :param q:   vector of absolute coordinates (positions and velocities)
        :
        """
        #   position (x, y)
        C = np.zeros(3)
        C[0:2] = self.evaluate_rijP(q)
        #   rotation (theta)
        C[2] = q2theta_i(q, self.body_id_i) - q2theta_i(q, self.body_id_j)
        return C

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
        C = np.zeros(3)
        if self.body_id_i != "ground" or self.body_id_j != "ground":
            C[0:2] = self.evaluate_rijP(q)
            C[2] = q2theta_i(q, self.body_id_i) - q2theta_i(q, self.body_id_j)
        elif self.body_id_i == "ground":
            C[0:2] = Ai_ui_P_vector(self.u_jP, q2theta_i(q, self.body_id_j))
            C[2] = q2theta_i(q, self.body_id_k)
        elif self.body_id_j == "ground":
            C[0:2] = Ai_ui_P_vector(self.u_iP, q2theta_i(q, self.body_id_i))
            C[2] = q2theta_i(q, self.body_id_i)
        else:
            raise ValueError, "Not right!"

        return C - self.C0

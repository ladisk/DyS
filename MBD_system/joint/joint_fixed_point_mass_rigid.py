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


class JointFixedPointMassRigid(Joint):
    """
    classdocs
    """
    def __init__(self, body_id_i, body_id_j, u_iP=np.zeros(2, dtype=float), u_jP=np.zeros(2, dtype=float), properties_dict={}, parent=None):
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
        self.joint_type = "fixed joint point mass-rigid"
        super(JointFixedPointMassRigid, self).__init__(self.joint_type, body_id_i, body_id_j, properties_dict=properties_dict, parent=parent)

        #   joint DOF
        self.joint_DOF = 0

        #   number of constrained (nodal) coordinates
        self.n_CNC = 2

        #   C_q matrix dimensions [rows, cols]
        self.C_q_dim = [self.n_CNC, [2, 3]]

        #   vectors
        self.u_iP = self.u_iP_LCS = u_iP
        self.u_jP = self.u_jP_LCS = u_jP

        self.u_P_LCS_list = [self.u_iP, self.u_jP]

        #   initial evaluation of constraint equation for fixed joint
        self.C0 = self.evaluate_C0(self.q0)

        #   list of markers
        self.markers = self._create_markers()

    def evaluate_C_q(self, q):
        """
        Function evaluates C_q matrix
        :return:
        """
        #    construct C_q matrix for fixed joint
        C_q_i = self._evaluate_C_q_i(q)

        C_q_j = self._evaluate_C_q_j(q)

        #   create list
        self.C_q_list = [C_q_i, C_q_j]

        return C_q_i, C_q_j

    def _evaluate_C_q_i(self, q):
        """

        :param q:
        :return:
        """
        C_q_i = np.zeros([self.C_q_dim[0], self.C_q_dim[1][0]])

        C_q_i[0:2, 0:2] = np.eye(2)

        return C_q_i

    def _evaluate_C_q_j(self, q):
        """

        :param q:
        :return:
        """
        C_q_j = np.zeros([self.n_CNC, self.C_q_dim[1][1]])

        C_q_j[0:2, 0:2] = -np.eye(2)

        C_q_j[:, -1] = Ai_theta_ui_P_vector(self.u_P_LCS_list[1], q2theta_i(q, self.body_id_j), 1)

        return C_q_j

    def evaluate_Q_d(self, q):
        """
        Q_d equation according to eq. 3.133 in Computational Dynamics 3rd Ed. (A.A.Shabana)
        expanded with z component of moment
        """
        Q_d = np.zeros(self.n_CNC)

        #   get uPi in LCS
        _u_P_j = self.u_P_LCS_list[self.body_id_list.index(self.body_id_j)]

        #   rigid body i data
        theta_j = q2theta_i(q, self.body_id_j)
        omega_j = q2dtheta_i(q, self.body_id_j)

        #   evaluate uPi in GCS
        u_P_i = -Ai_ui_P_vector(_u_P_j, theta_j)

        #   vector Qd
        Q_d[0:2] = u_P_i * (omega_j ** 2)

        # #    if body i and j are not ground
        # if self.body_id_i != "ground" and self.body_id_j != "ground":
        #     self.Q_d_list = []
        #     for body_id in self.body_id_list:
        #         Q_d_body = Joint_Q_d_vector(self.joint_type)
        #
        #         Q_d_body.Q_d[0:2] = Ai_ui_P_vector(self.u_P_LCS_list[Q_d_body.id], q2theta_i(q, body_id)) * (q2dtheta_i(q, body_id) ** 2)
        #
        #         #    add calculated joint matrix of a body to list
        #         self.Q_d_list.append(Q_d_body)
        #
        #     Q_d = +(self.Q_d_list[0].Q_d - self.Q_d_list[1].Q_d)
        #
        # elif self.body_id_i == "ground":
        #     Q_d_body = Joint_Q_d_vector(self.joint_type, body_connected_to_ground=True)
        #
        #     Q_d_body.id = 1
        #
        #     body_id_ = self.body_id_list[Q_d_body.id]
        #     #    calculate Q_d vector of only non-ground body
        #     Q_d_body.Q_d[0:2] = Ai_ui_P_vector(u_P=self.u_P_list[Q_d_body.id], theta=q2theta_i(q, body_id_)) * (q2dtheta_i(q, body_id_) ** 2)
        #
        #     Q_d = Q_d_body.Q_d
        #
        # elif self.body_id_j == "ground":
        #     Q_d_body = Joint_Q_d_vector(self.joint_type, body_connected_to_ground=True)
        #     body_id_ = self.body_id_list[Q_d_body.id]
        #
        #     #    calculate Q_d vector of only non-ground body
        #     Q_d_body.Q_d[0:2] = Ai_ui_P_vector(self.u_P_LCS_list[Q_d_body.id], q2theta_i(q, body_id_)) * (q2dtheta_i(q, body_id_) ** 2)
        #
        #     Q_d = Q_d_body.Q_d
        # else:
        #     raise ValueError, "Q_d vector for fixed joint not constructed. Check calculation process."

        return Q_d

    def evaluate_C0(self, q):
        """
        Function evaluates C vector for the joint
        :param q:   vector of absolute coordinates (positions and velocities)
        :
        """
        #   position (x, y)
        C = np.zeros(self.n_CNC)
        C[0:2] = self.evaluate_rijP(q)

        return C

    def evaluate_C(self, q, t=0.):
        """
        Function creates C vector for the joint
        Args:
        bodies - list of all bodies in MBD system
        q - vector of absolute coordinates (positions and velocities)
        Returns:
        C_q matrix for each body in joint
        Raises:
        """
        self.C = np.zeros(self.n_CNC)
        if self.body_id_i != "ground" or self.body_id_j != "ground":
            self.C[0:2] = self.evaluate_rijP(q)

        elif self.body_id_i == "ground":
            self.C[0:2] = Ai_ui_P_vector(self.u_jP, q2theta_i(q, self.body_id_j))

        elif self.body_id_j == "ground":
            self.C[0:2] = Ai_ui_P_vector(self.u_iP, q2theta_i(q, self.body_id_i))

        else:
            raise ValueError, "Not right!"

        self.C = self.C - self.C0
        return self.C


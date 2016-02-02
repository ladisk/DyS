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
        super(JointPrismatic, self).__init__(joint_type, body_id_i, body_id_j, u_iP_CAD=np.array([0, 0]), u_jP_CAD=np.array([0, 0]), u_iQ_CAD=u_iQ_CAD, parent=parent)

        #   joint DOF
        self.joint_DOF = 0

        #   evaluate initial angle between connected bodies
        self.theta0 = self.evaluate_C0(self.q0)

    def _evaluate_C_q(self, q):
        """
        Function evaluates C_q matrix
        :return:
        """


    def evaluate_Q_d(self, q):
        """
        Function evaluates Q_d vector of a joint
        Args:
            bodies - list of all bodies
            q - array of ppositions and velosities
            N - number of bodies
        Returns:
            Q_d - vector Q_d for joint
        """

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


    def evaluate_C0(self, q):
        """

        :param q:
        :return:
        """
        return q2theta_i(q, self.body_id_i) - q2theta_i(q, self.body_id_j)



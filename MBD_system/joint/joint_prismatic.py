# coding=utf-8


import numpy as np


from MBD_system.joint.joint import Joint


class JointPrismatic(Joint):
    """
    classdocs
    """
    def __init__(self, joint_type, body_id_i, body_id_j, u_iP_CAD=np.array([0, 0]), u_jP_CAD=np.array([0, 0]), u_iQ=np.array([0, 0]), parent=None):
        super(JointPrismatic, self).__init__(joint_type, body_id_i, body_id_j, u_iP_CAD=np.array([0, 0]), u_jP_CAD=np.array([0, 0]), u_iQ=np.array([0, 0]), parent=None)


    def create_joint_C_q_matrix(self, q, bodies):
        """
        Function evaluates C_q matrix
        :return:
        """

    def create_joint_Q_d_vector(self, bodies, q, N_b):
        """
        Function evaluates Q_d vector of a joint
        Args:
            bodies - list of all bodies
            q - array of ppositions and velosities
            N - number of bodies
        Returns:
            Q_d - vector Q_d for joint
        """

    def create_joint_C_vector(self, q):
        """
        Function creates C vector for the joint
        Args:
        bodies - list of all bodies in MBD system
        q - vector of absolute coordinates (positions and velocities)
        Returns:
        C_q matrix for each body in joint
        Raises:
        """

# coding=utf-8
"""

created by: lskrinjar
date of creation: 13/08/2016
time of creation: 11:38
"""




import numpy as np


from MBD_system.joint.joint import Joint
from global_variables import GlobalVariables


class SupportRoller(Joint):
    """
    classdocs
    """
    def __init__(self, joint_type, body_id_i, body_id_j, u_iP_CAD=np.array([0, 0]), u_jP_CAD=np.array([0, 0]), u_iQ=np.array([0, 0]), properties_dict={}, parent=None):
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
        self.joint_type = "roller support"
        super(SupportRoller, self).__init__(self.joint_type, body_id_i, body_id_j, u_iP_CAD=u_iP_CAD, u_jP_CAD=u_jP_CAD, properties_dict=properties_dict, parent=parent)

        #   joint DOF
        self.joint_DOF = 0

        #   initial evaluation of constraint equation for fixed joint
        self.C0 = self.evaluate_C0(self.q0)
        self.theta0 = self.C0[2]

    def evaluate_C_q(self, q):
        """
        Function evaluates C_q matrix for every element connected by this support (joint) object
        :param q:
        :return:
        """
        C_q_list = []


        GlobalVariables.q_i_dim[body_id]


    def evaluate_Q_d(self, q):
        """

        :param q:
        :return:
        """

    def evaluate_C0(self, q):
        """

        :param q:
        :return:
        """

    def evaluate_C(self, q, t=None):
        """

        :param q:
        :param t:
        :return:
        """
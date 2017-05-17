"""

created by: lskrinjar
date of creation: 01/12/2016
time of creation: 13:53
"""

# coding=utf-8

import numpy as np

from MBD_system.joint.joint import Joint


class Support(Joint):
    """
    classdocs
    """
    def __init__(self, support_type, body_id_i, body_id_j, node_id_i=None, node_id_j=None, properties_dict={}, parent=None):
        """

        :param support:
        :param body_id_i:
        :param body_id_j:
        :param u_iP_CAD:
        :param u_jP_CAD:
        :param u_iQ:
        :param parent:
        :return:
        """
        super(Support, self).__init__(self.joint_type, body_id_i, body_id_j, properties_dict=properties_dict, parent=parent)

        #   support DOF (defined in subclass)
        self.joint_DOF = None

        #   node data
        self.node_id_i = node_id_i
        self.node_id_j = node_id_j
        self.node_id_list = [self.node_id_i, self.node_id_j]

        #   list of constraint matrix for every body
        #   for deformable bodies this matrix C_q is constant is evaluated only once in preprocessor
        self.C_q_list = []

        #   number of constrained nodal coordinates per node of finite element
        self.n_CNC = None

        #   C_q matrix dimensions [rows, cols]
        self.C_q_dim = None

        #   initial evaluation of constraint equation for fixed joint
        self.C0 = None
        self.theta0 = None

    def evaluate_Q_d(self, q):
        """
        Function is defined in subclass
        :param q:
        :return:
        """
        Q_d = np.zeros(self.n_CNC)
        return Q_d

    def evaluate_C0(self, q):
        """
        Function is defined in subclass
        :param q:
        :return:
        """
        return None

    def evaluate_C(self, q, t=None):
        """
        Function is defined in subclass
        :param q:
        :param t:
        :return:
        """
        return None

if __name__ == "__main__":
    pass
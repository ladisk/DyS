# coding=utf-8
"""

created by: lskrinjar
date of creation: 11.08.2016
time of creation: 11:24
"""


import numpy as np


from MBD_system.joint.joint import Joint
from MBD_system.joint.joint_C_q_matrix import Joint_C_q_matrix
from global_variables import GlobalVariables
from MBD_system.joint.support import Support
from MBD_system.ancf.mesh.mesh import Mesh


class SupportFixed(Support):
    """
    classdocs
    """
    def __init__(self, body_id_i, body_id_j, node_id_i=0, node_id_j=0, properties_dict={}, parent=None):
        """

        :param body_id_i:
        :param body_id_j:
        :param u_iP_CAD:
        :param u_jP_CAD:
        :param u_iQ:
        :param parent:
        :return:
        """
        self.joint_type = "fixed support"
        super(SupportFixed, self).__init__(self.joint_type, body_id_i, body_id_j, properties_dict=properties_dict, parent=parent)

        #   joint DOF
        self.joint_DOF = 0

        #   boolean for constant C_q matrix of joint
        self.constant = True

        #   node data
        self.node_id_i = node_id_i
        self.node_id_j = node_id_j
        self.node_id_list = [self.node_id_i, self.node_id_j]

        #   number of constrained nodal coordinates of a support
        for body_id in self.body_id_list:
            if type(body_id) is int:
                if isinstance(self._parent._parent.bodies[body_id].mesh, Mesh):
                    e_n = self._parent._parent.bodies[body_id].mesh.element_e_n

                #   number of constrained nodal coordinates by support
                self.n_CNC = e_n / 2

                #   C_q matrix dimensions [rows, cols]
                self.e_n_list = [e_n, e_n]
                self.C_q_dim = [self.n_CNC, self.e_n_list]

        #   list of constraint matrix for every body
        #   for deformable bodies this matrix C_q is constant is evaluated only once in preprocessor
        self.C_q_list = []

        #   initial evaluation of constraint equation for fixed joint
        # self.C0 = self.evaluate_C0(self.q0)
        # self.theta0 = self.C0[2]

    def evaluate_C_q(self, q):
        """
        Function evaluates C_q matrix for every element connected by this support (joint) object
        :param q:
        :return:
        """
        self.C_q_list = []
        #    construct C_q matrix for fixed support (joint)
        for body_id, node_id in zip(self.body_id_list, self.node_id_list):
            if body_id == "ground":
                C_q = Joint_C_q_matrix(self.joint_type, parent=self)

            else:
                C_q = Joint_C_q_matrix(self.joint_type, body_id=body_id, parent=self)

                C_q.matrix = self._parent._parent.bodies[body_id].evaluate_C_q_fixed(node_id)

            self.C_q_list.append(C_q)

        [C_qi, C_qj] = self.C_q_list

        return C_qi.matrix, C_qj.matrix

    def evaluate_Q_d(self, q):
        """

        :param q:
        :return:
        """
        Q_d = np.zeros(self.n_CNC)
        return Q_d

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

    def evaluate_Q_c(self, L):
        """
        Function evaluates
        :param L:    vector of lagrange multipliers of a joint
        """
        return np.zeros(self.n_CNC)



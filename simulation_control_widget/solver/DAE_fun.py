"""
Created on 13. mar. 2014

@author: luka.skrinjar
"""
import time
import logging
from pprint import pprint
import numpy as np

from MBD_system import inverse_blockwise
from gaussian_elimination import gaussian_elimination
from global_variables import GlobalVariables


class DAEfun(object):
    """
    classdocs
    """

    def __init__(self, MBD_system, parent=None):
        """
        DAE - Automatic generation of Differential Algebraic Equations
        """
        #   parent
        self._parent = parent
        #   pointer to MBD system as object attribute
        self.MBD_system = MBD_system

        #    gravity
        if MBD_system is not None:
            self.gravity = self.MBD_system.gravity_magnitude * self.MBD_system.gravity_vector
        else:
            self.gravity = np.zeros(2, dtype=float)

        #   variables during integration
        #   counter attribute
        self.number_of_evaluations = 0
        self.Q_g = None
        self.Q_e = None
        self.Q_s = None
        self.Q_dis = None
        self.Q_q = None
        self.Q_dq = None
        self.M_dim = None
        self.rows_C_s = None
        self.cols_C_s = None
        self.q_dim = None
        self.cols_C_q = None
        self.rows_C_q = None
        self.C_q = None
        self.C_qT = None
        self.J_dim = None
        self.M = None
        self.M_inv = None
        self.M_created = False
        self.error = False

        #   vector of system unknowns
        self.q = []
        self.dq = []
        self.ddq = []
        self.L = []

    def _C_q_size(self, M_dim=None):
        """
        Evaluate MBD system matrix C_q size
        :return:
        """
        if M_dim is None:
            cols_C_q = self.MBD_system.evaluate_M_size()
        else:
            cols_C_q = M_dim

        rows_C_q = self.MBD_system.evaluate_C_number_of_rows()

        return cols_C_q, rows_C_q

    def preprocessing(self):
        """
        Create mass matrix
        """
        #   size of mass matrix of MBD system
        self.M_dim = self.MBD_system.evaluate_M_size()

        #   size of vector q (y, dy) of MBD system, twice the size of mass matrix
        self.q_dim = 2 * self.M_dim

        self.cols_C_q, self.rows_C_q = self._C_q_size(M_dim=self.M_dim)

        self.J_dim = self.cols_C_q + self.rows_C_q

        #    predefine mass matrix size
        M = np.zeros([self.M_dim, self.M_dim])
        _row = 0
        for body in self.MBD_system.bodies:
            M[_row:_row + body.M_size, _row:_row + body.M_size] = body.evaluate_M()
            _row += body.M_size

        self.M = M
        self.M_created = True

        #    inverse mass matrix
        self.M_inv = np.linalg.inv(self.M)

        #   evaluate vector of gravitational forces as it is constant
        self.Q_g = self.evaluate_Q_g(q=self.MBD_system.q0)

        #   check for constant stiffness matrix and evaluate if
        for body in self.MBD_system.bodies:
            if hasattr(body, "mesh"):
                if body.mesh is not None:
                    if body.mesh.K_constant_matrix:
                        body.mesh._construct_constant_matrix_K_t()

        #   set data of forces on flexible bodies - meshes
        self._set_forces_properties()

    def _set_forces_properties(self):
        """

        :return:
        """
        for force in self.MBD_system.forces:
            if force._body.body_type in ["flexible body", "finite element"]:
                force.set_element()

    def evaluate_K(self, q=None):
        """
        Evaluate stiffness matrix of the system if there are any ANCF finite elements
        :return:
        """
        for body in self.MBD_system.bodies:
            if hasattr(body, "evaluate_K"):
                body.evaluate_K(q)

    def getDOF(self, q):
        """
        Evaluate degrees of freedom
        """
        q0 = self.MBD_system.create_q0()
#         print "number of coordinates =", len(self.MBD_system.q0) / 2
        C_q, C_qT = self.evaluate_C_q(q)
#         print "DOFs =", np.linalg.matrix_rank(C_q, tol=None)
        #    num_Cq_eq - number of constraint equations
        num_Cq_eq, num_C = C_q.shape

        DOF = num_C - num_Cq_eq
        print "DOF =", DOF

        _q0 = q0[0:0.5 * len(q0)]
        gaussian_elimination(C_q, None)

        return DOF

    def evaluate_C(self, q, t=None):
        """
        Function calculates C vector for every joint - this is used with Baumgarte's Stabilization Method - BSM
        """
        C = np.zeros(self.MBD_system.evaluate_C_number_of_rows())
        for joint in self.MBD_system.joints:
            C_joint = joint.evaluate_C(q, t=t)

            if joint.joint_type in joint._supported_types:
                rows = sum(C_q_i_dim_i[0] for C_q_i_dim_i in GlobalVariables.C_q_i_dim[0:joint.joint_id])
                C[rows:rows + joint.n_CNC] = C_joint
            #
            # if joint.joint_type == "fixed":
            #     C[3 * joint.fixed_joint_id:3 * joint.fixed_joint_id + 3] = C_joint
            #
            # elif joint.joint_type == "revolute":
            #     C[3 * self.MBD_system.number_of_fixed_joints + 2 * joint.revolute_joint_id:3 * self.MBD_system.number_of_fixed_joints + 2 * joint.revolute_joint_id + 2] = C_joint
            #
            # elif joint.joint_type == "prismatic":
            #     C[3 * self.MBD_system.number_of_fixed_joints + 2 * self.MBD_system.number_of_revolute_joints + 2 * joint.prismatic_joint_id:3 * self.MBD_system.number_of_fixed_joints + 2 * self.MBD_system.number_of_revolute_joints + 2 * joint.prismatic_joint_id + 2] = C_joint
            #
            # elif joint.joint_type == "hinged support":
            #     rows = sum(C_q_i_dim_i[0] for C_q_i_dim_i in GlobalVariables.C_q_i_dim[0:joint.joint_id])
            #     C[rows:rows + joint.n_CNC] = C_joint
            #
            # elif joint.joint_type == "revolute joint rigid-flexible":
            #     rows = sum(C_q_i_dim_i[0] for C_q_i_dim_i in GlobalVariables.C_q_i_dim[0:joint.joint_id])
            #     C[rows:rows + joint.n_CNC] = C_joint
            #
            # elif joint.joint_type == "revolute joint flexible-flexible":
            #     rows = sum(C_q_i_dim_i[0] for C_q_i_dim_i in GlobalVariables.C_q_i_dim[0:joint.joint_id])
            #     C[rows:rows + joint.n_CNC] = C_joint
            #
            # elif joint.joint_type == "rigid joint flexible-flexible":
            #     rows = sum(C_q_i_dim_i[0] for C_q_i_dim_i in GlobalVariables.C_q_i_dim[0:joint.joint_id])
            #     C[rows:rows + joint.n_CNC] = C_joint
            #
            # elif joint.joint_type == "rigid joint rigid-flexible":
            #     rows = sum(C_q_i_dim_i[0] for C_q_i_dim_i in GlobalVariables.C_q_i_dim[0:joint.joint_id])
            #     C[rows:rows + joint.n_CNC] = C_joint
            #
            # elif joint.joint_type == "rigid joint point mass-flexible":
            #     rows = sum(C_q_i_dim_i[0] for C_q_i_dim_i in GlobalVariables.C_q_i_dim[0:joint.joint_id])
            #     C[rows:rows + joint.n_CNC] = C_joint

            else:
                print "joint.joint_type =", joint.joint_type
                raise ValueError, "C_q matrix not constructed. Unknown joint type!"

        for i, motion in enumerate(self.MBD_system.motions):
            C[self.MBD_system.C_q_number_of_rows + i] = motion.evaluate_C(q, t=t)

        return C

    def check_C(self, q, t):
        """
        Function evaluates C if tolerance is violated
        :param q:
        :param t:
        :return:
        """
        C = self.evaluate_C(q, t)

        if self._parent.errorControl:
            #   all elements of vector C have to be under specified tolerance defined by used
            # print "(abs(C) >= self.MBD_system.TOL_C).any() =", (abs(C) >= self.MBD_system.TOL_C).any()
            if (abs(C) >= self.MBD_system.TOL_C).any():
                for i, joint in enumerate(self.MBD_system.joints):
                    print i, joint._name, "C =", joint.C
                    print "r_P_GCS_list =", joint.r_P_GCS_list

                info = "C(q, t) check - tolerances violated!"
                logging.getLogger("DyS_logger").info(info)
                if self._parent is not None:
                    if hasattr(self._parent, "simulation_error"):
                        self._parent.simulation_error.setError(info)

                return True

        return False

    def evaluate_C_kinematic_analysis(self, q, t=None):
        """

        :param q:
        :param t:
        :return:
        """
        return -self.evaluate_C(q, t)

    def evaluate_C_q(self, q, t):
        """
        Evaluate matrix C_q
        """
        #    C_q matrix if there are no contacts in compression phase and only constraint equations
        self.cols_C_q, self.rows_C_q = self._C_q_size()

        C_q = np.zeros([self.rows_C_q, self.cols_C_q])

        for i, joint in enumerate(self.MBD_system.joints):
            #    creates a C_q matrix of a joint based on joint's properties (type, geometry, etc.)
            C_q_i, C_q_j = joint.evaluate_C_q(q)
            if joint.joint_type in joint._supported_types:
                for _C_q_i, body_id, _C_q_i_cols in zip([C_q_i, C_q_j], joint.body_id_list, joint.C_q_dim[1]):
                    if _C_q_i is not None and isinstance(body_id, int):
                        rows = sum(C_q_i_dim_i[0] for C_q_i_dim_i in GlobalVariables.C_q_i_dim[0:joint.joint_id])
                        cols = np.sum(GlobalVariables.q_i_dim[0:body_id])
                        C_q[rows:rows + joint.n_CNC, cols:cols + _C_q_i_cols] = _C_q_i

            else:
                print "joint.joint_type =", joint.joint_type
                raise ValueError, "C_q matrix not constructed. Unknown joint type!"
            #
            # #    C_q matrix for joint type: fixed
            # if joint.joint_type == "fixed":
            #     #    add matrix of body j to Cq matrix
            #     if joint.body_id_i == "ground":
            #         C_q[np.sum(GlobalVariables.q_i_dim[0:joint.body_id]):np.sum(GlobalVariables.q_i_dim[0:joint.body_id + 1])] = C_q_j
            #     #    add matrix of body i to Cq matrix
            #     elif joint.body_id_j == "ground":
            #         C_q[3 * joint.fixed_joint_id:3 * joint.fixed_joint_id + 3, 3 * joint.body_id_i:3 * joint.body_id_i + 3] = C_q_i
            #     #    add matrix of body i and j to Cq matrix
            #     elif (joint.body_id_i != "ground") and (joint.body_id_j != "ground"):
            #         C_q[3 * joint.fixed_joint_id:3 * joint.fixed_joint_id + 3, 3 * joint.body_id_i:3 * joint.body_id_i + 3] = C_q_i
            #         C_q[3 * joint.fixed_joint_id:3 * joint.fixed_joint_id + 3, 3 * joint.body_id_j:3 * joint.body_id_j + 3] = C_q_j
            #
            # #    C_q matrix for joint type: revolute
            # elif joint.joint_type == "revolute":
            #     #    add matrix of body i and j to Cq matrix
            #     if (joint.body_id_i != "ground") and (joint.body_id_j != "ground"):
            #         C_q[3 * self.MBD_system.number_of_fixed_joints + 2 * joint.revolute_joint_id:3 * self.MBD_system.number_of_fixed_joints + 2 * joint.revolute_joint_id + 2, 3 * joint.body_id_i:3 * joint.body_id_i + 3] = C_q_i
            #         C_q[3 * self.MBD_system.number_of_fixed_joints + 2 * joint.revolute_joint_id:3 * self.MBD_system.number_of_fixed_joints + 2 * joint.revolute_joint_id + 2, 3 * joint.body_id_j:3 * joint.body_id_j + 3] = C_q_j
            #     #    add matrix of body j to Cq matrix
            #     elif joint.body_id_i == "ground":
            #         C_q[3 * self.MBD_system.number_of_fixed_joints + 2 * joint.revolute_joint_id:3 * self.MBD_system.number_of_fixed_joints + 2 * joint.revolute_joint_id + 2, 3 * joint.body_id_j:3 * joint.body_id_j + 3] = C_q_j
            #     #    add matrix of body i to Cq matrix
            #     elif joint.body_id_j == "ground":
            #         C_q[3 * self.MBD_system.number_of_fixed_joints + 2 * joint.revolute_joint_id:3 * self.MBD_system.number_of_fixed_joints + 2 * joint.revolute_joint_id + 2, 3 * joint.body_id_i:3 * joint.body_id_i + 3] = C_q_i
            #
            # #    C_q matrix for joint type: prismatic
            # elif joint.joint_type == "prismatic":
            #     #    add matrix of body i and j to Cq matrix
            #     if (joint.body_id_i != "ground") and (joint.body_id_j != "ground"):
            #         C_q[3 * self.MBD_system.number_of_fixed_joints + 2 * self.MBD_system.number_of_revolute_joints + 2 * joint.prismatic_joint_id:3 * self.MBD_system.number_of_fixed_joints + 2 * self.MBD_system.number_of_revolute_joints + 2 * joint.prismatic_joint_id + 2, 3 * joint.body_id_i:3 * joint.body_id_i + 3] = C_q_i
            #         C_q[3 * self.MBD_system.number_of_fixed_joints + 2 * self.MBD_system.number_of_revolute_joints + 2 * joint.prismatic_joint_id:3 * self.MBD_system.number_of_fixed_joints + 2 * self.MBD_system.number_of_revolute_joints + 2 * joint.prismatic_joint_id + 2, 3 * joint.body_id_j:3 * joint.body_id_j + 3] = C_q_j
            #     #    add matrix of body j to Cq matrix
            #     elif joint.body_id_i == "ground":
            #         C_q[3 * self.MBD_system.number_of_fixed_joints + 2 * self.MBD_system.number_of_revolute_joints + 2 * joint.prismatic_joint_id:3 * self.MBD_system.number_of_fixed_joints + 2 * self.MBD_system.number_of_revolute_joints + 2 * joint.prismatic_joint_id + 2, 3 * joint.body_id_j:3 * joint.body_id_j + 3] = C_q_j
            #     #    add matrix of body i to Cq matrix
            #     elif joint.body_id_j == "ground":
            #         C_q[3 * self.MBD_system.number_of_fixed_joints + 2 * self.MBD_system.number_of_revolute_joints + 2 * joint.prismatic_joint_id:3 * self.MBD_system.number_of_fixed_joints + 2 * self.MBD_system.number_of_revolute_joints + 2 * joint.prismatic_joint_id + 2, 3 * joint.body_id_i:3 * joint.body_id_i + 3] = C_q_i
            #
            # elif joint.joint_type in ["fixed support", "hinged support"]:
            #     for _C_q_i, body_id in zip([C_q_i, C_q_j], joint.body_id_list):
            #         if _C_q_i is not None:
            #             C_q[sum(C_q_i_dim_i[0] for C_q_i_dim_i in GlobalVariables.C_q_i_dim[0:joint.joint_id]):sum(C_q_i_dim_i[0] for C_q_i_dim_i in GlobalVariables.C_q_i_dim[0:joint.joint_id+1]), np.sum(GlobalVariables.q_i_dim[0:body_id]):np.sum(GlobalVariables.q_i_dim[0:body_id]) + joint.C_q_dim[1]] = _C_q_i
            #
            # elif joint.joint_type in ["slope discontinuity"]:
            #     cols = np.sum(GlobalVariables.q_i_dim[0:joint.body_id]) + self.MBD_system.bodies[joint.body_id].mesh.node_dim[joint.node_id] + 0
            #     C_q[sum(C_q_i_dim_i[0] for C_q_i_dim_i in GlobalVariables.C_q_i_dim[0:joint.joint_id]):sum(C_q_i_dim_i[0] for C_q_i_dim_i in GlobalVariables.C_q_i_dim[0:joint.joint_id + 1]),
            #     cols:cols + len(np.hstack(joint.C_q_list))] = np.hstack(joint.C_q_list)
            #
            # elif joint.joint_type == "revolute joint rigid-flexible":
            #     for _C_q_i, body_id, _C_q_i_cols in zip([C_q_i, C_q_j], joint.body_id_list, joint.C_q_dim[1]):
            #         rows = sum(C_q_i_dim_i[0] for C_q_i_dim_i in GlobalVariables.C_q_i_dim[0:joint.joint_id])
            #         cols = np.sum(GlobalVariables.q_i_dim[0:body_id])
            #         C_q[rows:rows + joint.n_CNC, cols:cols + _C_q_i_cols] = _C_q_i
            #
            # elif joint.joint_type == "revolute joint flexible-flexible":
            #     for _C_q_i, body_id, _C_q_i_cols in zip([C_q_i, C_q_j], joint.body_id_list, joint.C_q_dim[1]):
            #         rows = sum(C_q_i_dim_i[0] for C_q_i_dim_i in GlobalVariables.C_q_i_dim[0:joint.joint_id])
            #         cols = np.sum(GlobalVariables.q_i_dim[0:body_id])
            #         C_q[rows:rows + joint.n_CNC, cols:cols + _C_q_i_cols] = _C_q_i
            #
            # elif joint.joint_type == "rigid joint flexible-flexible":
            #     for _C_q_i, body_id, _C_q_i_cols in zip([C_q_i, C_q_j], joint.body_id_list, joint.C_q_dim[1]):
            #         rows = sum(C_q_i_dim_i[0] for C_q_i_dim_i in GlobalVariables.C_q_i_dim[0:joint.joint_id])
            #         cols = np.sum(GlobalVariables.q_i_dim[0:body_id])
            #         C_q[rows:rows + joint.n_CNC, cols:cols + _C_q_i_cols] = _C_q_i
            #
            # elif joint.joint_type == "rigid joint rigid-flexible":
            #     for _C_q_i, body_id, _C_q_i_cols in zip([C_q_i, C_q_j], joint.body_id_list, joint.C_q_dim[1]):
            #         rows = sum(C_q_i_dim_i[0] for C_q_i_dim_i in GlobalVariables.C_q_i_dim[0:joint.joint_id])
            #         cols = np.sum(GlobalVariables.q_i_dim[0:body_id])
            #         C_q[rows:rows + joint.n_CNC, cols:cols + _C_q_i_cols] = _C_q_i
            #
            # elif joint.joint_type == "rigid joint point mass-flexible":
            #     for _C_q_i, body_id, _C_q_i_cols in zip([C_q_i, C_q_j], joint.body_id_list, joint.C_q_dim[1]):
            #         rows = sum(C_q_i_dim_i[0] for C_q_i_dim_i in GlobalVariables.C_q_i_dim[0:joint.joint_id])
            #         cols = np.sum(GlobalVariables.q_i_dim[0:body_id])
            #         C_q[rows:rows + joint.n_CNC, cols:cols + _C_q_i_cols] = _C_q_i
            #
            # else:
            #     raise ValueError, "C_q matrix not constructed. Unknown joint type!"

        for i, motion in enumerate(self.MBD_system.motions):
            C_q_motion = motion.evaluate_C_q(q)
            C_q[self.MBD_system.C_q_number_of_rows + i, :] = C_q_motion

        #    size of C_q matrix
        [self.rows_C_q, self.col_C_q] = C_q.shape

        return C_q

    def evaluate_C_t(self, q, t):
        """
        Function evaluates vector C_t of MBD system
        :param q:
        :param t:
        :return:
        """
        C_t = np.zeros(self.MBD_system.evaluate_C_number_of_rows())

        for i, motion in enumerate(self.MBD_system.motions):
            C_t[self.MBD_system.C_q_number_of_rows + i] = motion.evaluate_C_t(q, t=t)

        return C_t

    def evaluate_Q_d(self, q):
        """
        Evaluates vector of terms that are quadratic in the velocities
        """
        Q_d = np.zeros(self.rows_C_q)

        for joint in self.MBD_system.joints:

            Q_d_joint_vector = joint.evaluate_Q_d(q)


            if joint.joint_type in joint._supported_types:

                rows = sum(C_q_i_dim_i[0] for C_q_i_dim_i in GlobalVariables.C_q_i_dim[0:joint.joint_id])
                Q_d[rows:rows + joint.n_CNC] = Q_d_joint_vector

            else:
                print "joint.joint_type =", joint.joint_type
                raise ValueError, "C_q matrix not constructed. Unknown joint type!"


            #
            # if joint.joint_type == "fixed":
            #     Q_d[3 * joint.fixed_joint_id:3 * joint.fixed_joint_id + 3] = Q_d_joint_vector
            #
            # elif joint.joint_type == "revolute":
            #     Q_d[3 * self.MBD_system.number_of_fixed_joints + 2 * joint.revolute_joint_id:3 * self.MBD_system.number_of_fixed_joints + 2 * joint.revolute_joint_id + 2] = Q_d_joint_vector
            #
            # elif joint.joint_type == "prismatic":
            #     Q_d[3 * self.MBD_system.number_of_fixed_joints + 2 * self.MBD_system.number_of_revolute_joints + 2 * joint.prismatic_joint_id:3 * self.MBD_system.number_of_fixed_joints + 2 * self.MBD_system.number_of_revolute_joints + 2 * joint.prismatic_joint_id + 2] = Q_d_joint_vector
            #
            # elif joint.joint_type == "hinged support":
            #     rows = sum(C_q_i_dim_i[0] for C_q_i_dim_i in GlobalVariables.C_q_i_dim[0:joint.joint_id])
            #     Q_d[rows:rows + joint.n_CNC] = Q_d_joint_vector
            #
            # elif joint.joint_type == "revolute joint rigid-flexible":
            #     rows = sum(C_q_i_dim_i[0] for C_q_i_dim_i in GlobalVariables.C_q_i_dim[0:joint.joint_id])
            #     Q_d[rows:rows + joint.n_CNC] = Q_d_joint_vector
            #
            # elif joint.joint_type == "revolute joint flexible-flexible":
            #     rows = sum(C_q_i_dim_i[0] for C_q_i_dim_i in GlobalVariables.C_q_i_dim[0:joint.joint_id])
            #     Q_d[rows:rows + joint.n_CNC] = Q_d_joint_vector
            #
            # elif joint.joint_type == "rigid joint flexible-flexible":
            #     rows = sum(C_q_i_dim_i[0] for C_q_i_dim_i in GlobalVariables.C_q_i_dim[0:joint.joint_id])
            #     Q_d[rows:rows + joint.n_CNC] = Q_d_joint_vector
            #
            # elif joint.joint_type == "rigid joint rigid-flexible":
            #     rows = sum(C_q_i_dim_i[0] for C_q_i_dim_i in GlobalVariables.C_q_i_dim[0:joint.joint_id])
            #     Q_d[rows:rows + joint.n_CNC] = Q_d_joint_vector
            #
            # elif joint.joint_type == "rigid joint point mass-flexible":
            #     rows = sum(C_q_i_dim_i[0] for C_q_i_dim_i in GlobalVariables.C_q_i_dim[0:joint.joint_id])
            #     Q_d[rows:rows + joint.n_CNC] = Q_d_joint_vector
            #
            # else:
            #     print "Unknown joint type! Joint type is %s"%joint.joint_type

        return Q_d

    def evaluate_Q_g(self, q=None):
        """
        Evaluate vector of gravitational forces
        :argument:  q   vector of positions and velocities of MBD system
        :return:
        """
        Q_g = np.zeros(self.M_dim)
        #    add gravity forces
        if (self.gravity != np.zeros(3)).any():
            for body in self.MBD_system.bodies:
                Q_g_i = body.evaluate_Q_g(self.gravity)
                Q_g[sum(GlobalVariables.q_i_dim[0:body.body_id]):sum(GlobalVariables.q_i_dim[0:body.body_id + 1])] = Q_g_i

        return Q_g

    def evaluate_Q_e(self, t, q):
        """
        Evaluate vector of generalized external forces of the system
        """
        if self._parent.FLAG_contact == 1:
            self.evaluate_contacts(t, q)

        #    vector of generalised external forces
        Q_e = np.zeros(self.M_dim)
        #   evaluate a generalized force vector of all external forces
        #   contact forces are included
        for i, force in enumerate(self.MBD_system.forces):
            if force.active:
                Q_e_f = force.evaluate_Q_e(t, q)
                if Q_e_f is not None:
                    Q_e[np.sum(GlobalVariables.q_i_dim[0:force.body_id]):np.sum(GlobalVariables.q_i_dim[0:force.body_id + 1])] += Q_e_f

        #   evaluate a generalized force of spring element (translational and rotational)
        for i, spring in enumerate(self.MBD_system.springs):
            if spring.active:
                Q_i, Q_j = spring.evaluate_Q_e(q)
                for body_id, Q_e_ij in zip(spring.body_id_list, [Q_i, Q_j]):
                    if isinstance(body_id, int) and Q_e_ij is not None:
                        Q_e[np.sum(GlobalVariables.q_i_dim[0:body_id]):np.sum(GlobalVariables.q_i_dim[0:body_id + 1])] += Q_e_ij

        return Q_e

    def evaluate_Q_s(self, q):
        """
        Function evaluates elastic strain forces of finite element meshes of each body
        :return:
        """
        #    vector of generalised strain forces
        Q_s = np.zeros(self.M_dim)
        for body in self.MBD_system.bodies:
            if hasattr(body, "evaluate_Q_s"):
                # q_i = q[np.sum(GlobalVariables.q_i_dim[0:body.body_id]):np.sum(GlobalVariables.q_i_dim[0:body.body_id + 1])]
                Q_s_i = body.evaluate_Q_s(q)

                Q_s[sum(GlobalVariables.q_i_dim[0:body.body_id]):sum(GlobalVariables.q_i_dim[0:body.body_id + 1])] = Q_s_i

        return Q_s

    def evaluate_C_s(self, q):
        """
        Function creates C_s matrix of a contact
        """
        #    predefine C_q matrix size
        self.cols_C_s = 2
        C_s = np.zeros([self.rows_C_q_contacts, self.cols_C_s * self.active_contacts])

        #    create C_s matrix
        for i, contact in enumerate(self.MBD_system.contacts):
            if contact.contact_detected and contact.contact_distance_inside_tolerance:
                #    creates a C_q matrix of a joint based on joint's properties (type, geometry, etc.)
                contact_matrix_body_i, contact_matrix_body_j = contact.create_contact_C_s_matrix()

                C_s[i:i + 3, 2 * i] = contact_matrix_body_i.matrix[:, 0]
                C_s[i:i + 3, 2 * i + 1] = contact_matrix_body_j.matrix[:, 0]

#                 C_s[self.rows_C_q + i:self.rows_C_q + i + 3, 2 * i] = contact_matrix_body_i.matrix[:, 0]
#                 C_s[self.rows_C_q + i:self.rows_C_q + i + 3, 2 * i + 1] = contact_matrix_body_j.matrix[:, 0]

        #    size of C_q matrix
        [self.rows_C_s, self.cols_C_s] = C_s.shape

        #    matrix C_q transposed
        C_s_trans = C_s.T
        return C_s, C_s_trans

    def evaluate_dC(self, C_q, q):
        """

        """
        dq = q[0.5 * self.q_dim:]
        dC = np.dot(C_q, dq)
        return dC

    def evaluate_Q_dis(self, q):
        """

        :param q:
        :return:
        """
        #    vector of generalised strain forces
        Q_dis = np.zeros(self.M_dim)
        for body in self.MBD_system.bodies:
            if hasattr(body, "evaluate_Q_dis") and body.dissipation_coefficient is not None:
                Q_dis_i = body.evaluate_Q_dis(q)

                Q_dis[sum(GlobalVariables.q_i_dim[0:body.body_id]):sum(GlobalVariables.q_i_dim[0:body.body_id + 1])] = Q_dis_i

        return Q_dis

    def evaluate_dq(self, t, q):
        """
        Function evaluates vector dq of MBD system
        """
        #   function evaluation counter
        self.number_of_evaluations += 1

        #    construct Cq matrix
        self.C_q = self.evaluate_C_q(q, t)

        #    matrix C_q transposed
        self.C_qT = self.C_q.T

        #   size of C_s matrix
        self.rows_C_s = 0
        self.cols_C_s = 0

        #    construct augmented matrix
        #    size of predefined empty matrix
        rows = self.M_dim + self.rows_C_q
        cols = self.M_dim + self.rows_C_q

        #    construct vector of external and constraint forces
        _vector = np.empty([rows])

        #   vector of gravitational forces
        self.Q_g = self.evaluate_Q_g(q=q)
        # print "Q_g =", self.Q_g

        #    vector of external and contact forces
        self.Q_e = self.evaluate_Q_e(t, q)
        # print "Q_e =", self.Q_e

        #   vector of elastic strain forces
        self.Q_s = self.evaluate_Q_s(q)
        # print "self.Q_s =", self.Q_s

        #   vector of dissipative forces
        self.Q_dis = self.evaluate_Q_dis(q)

        _vector[0:self.M_dim] = self.Q_g + self.Q_e - self.Q_s - self.Q_dis

        #    vector of elements, that are quadratic in velocities
        if self.MBD_system.use_BSM:
            C = self.evaluate_C(q, t)
            dC = self.evaluate_dC(self.C_q, q)

            self.alpha = self._parent.h
            self.beta = self.alpha #* np.sqrt(2)

            Q_d = self.evaluate_Q_d(q) - 2*self.alpha*dC - (self.beta**2)*C

        else:
            Q_d = self.evaluate_Q_d(q)

        _vector[self.M_dim:rows] = Q_d
        # print "_vector ="
        # print _vector

        B = np.zeros([self.col_C_q, self.rows_C_q])
        B = self.C_qT

        C = np.zeros([self.rows_C_q, self.col_C_q])
        C = self.C_q

        D = np.zeros([self.rows_C_q, self.rows_C_q])

        #    inverse matrix block-wise
        # print "self.M_inv =", self.M_inv
        # print "B =", B
        # print "C =", C
        # print "D =", D
        # print "self.M_dim =", self.M_dim
        if B != [] and C != [] and D != []:
            _matrix_inverse = inverse_blockwise.inverse_blockwise(self.M_inv, B, C, D, self.M_dim)
        else:
            _matrix_inverse = self.M_inv

        #    solve system of equations
        _sol_vector = np.dot(_matrix_inverse, _vector)

        #    accelerations vector
        ddq = _sol_vector[0:self.M_dim]

        #    vector of lagrange multipliers
        self.L = _sol_vector[self.M_dim + self.cols_C_s:rows]
        
        #    evaluate vector of constraint forces for the system
        self._evaluate_Q_c()

        #    the dq vector in constructed from vector ddq and the second part of the input vector q
        dq = np.append(q[0.5 * self.q_dim:], ddq)

        #   check C
        self.error = self.check_C(q=q, t=t)
        return dq

    def _evaluate_Q_c(self):
        """
        Function evaluates vector of constraint reaction forces and assigns vector to each joint
        """
        #    generalized reaction forces associated with C_q matrix
        for joint in self.MBD_system.joints:
            if joint.joint_type in ["fixed", "revolute", "prismatic"]:
                if joint.joint_type == "fixed":
                    L = self.L[3 * joint.fixed_joint_id:3 * joint.fixed_joint_id + 3]

                if joint.joint_type == "revolute":
                    L = self.L[3 * self.MBD_system.number_of_fixed_joints + 2 * joint.revolute_joint_id:3 * self.MBD_system.number_of_fixed_joints + 2 * joint.revolute_joint_id + 2]

                if joint.joint_type == "prismatic":
                    L = self.L[3 * self.MBD_system.number_of_fixed_joints + 2 * self.MBD_system.number_of_revolute_joints + 2 * joint.prismatic_joint_id:3 * self.MBD_system.number_of_fixed_joints + 2 * self.MBD_system.number_of_revolute_joints + 2 * joint.prismatic_joint_id + 2]

                joint.evaluate_Q_c(L)
            
    def evaluate_contacts(self, t, q):
        """
        Solve all contact that are detected at integration time t
        """
        #   loop through all contacts
        for contact in self.MBD_system.contacts:
            #   check if contact is detected and calculated contact distance is inside the tolerance, then contact is present
            if contact.contact_detected and contact.status == 1 and contact.active:#and contact.contact_distance_inside_tolerance
                #   get contact forces
                contact.evaluate_contact(t, q)

    def evaluate_Jacobian(self, t, q):
        """

        :param t:
        :param q:
        :return:
        """
        J = np.zeros([self.J_dim, self.J_dim])

        #   evaluate jacobian matrix of constraint equations
        self.C_q = self.evaluate_C_q(q, t)

        #    matrix C_q transposed
        self.C_qT = self.C_q.T

        J[self.rows_C_q::, 0:self.cols_C_q] = self.C_qT
        J[0:self.cols_C_q, self.rows_C_q::] = self.C_q

        #   numerical differentiation
        J[0:self.rows_C_q, 0:self.cols_C_q] = None
        return J

    def evaluate_Q_q(self, q):
        """

        :return:
        """
        Q_q = np.zeros_like(self.M)

        for i, spring in enumerate(self.MBD_system.springs):
            pass

        for i, force in enumerate(self.MBD_system.forces):
            pass

    def evaluate_Q_dq(self, q):
        """

        :return:
        """
        Q_dq = np.zeros_like(self.M)

        for i, spring in enumerate(self.MBD_system.springs):
            pass

        for i, force in enumerate(self.MBD_system.forces):
            pass

    def evaluate_M_ddq(self):
        """

        :return:
        """
        M_ddq = np.zeros_like(self.M)
        return M_ddq

    def evaluate_C_q_L_q(self, q, L):
        """

        :return:
        """




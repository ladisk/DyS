"""
Created on 13. mar. 2014

@author: luka.skrinjar
"""
import logging
from pprint import pprint
import sys
import time
import cProfile

import numpy as np
import scipy as sp
from scipy.sparse import linalg
from scipy.sparse import find
from scipy.sparse import issparse
from scipy.sparse import csr_matrix
from scipy.sparse import csc_matrix


from MBD_system import inverse_blockwise
from gaussian_elimination import gaussian_elimination
from MBD_system.force.force import Force


class DAEfun(object):
    """
    classdocs
    """

    def __init__(self, MBD_system, parent=None):
        """
        DAE - Automatic generation of Differential Algebraic Equations
        """
        # super(DAEfun, self).__init__(MBD_system, parent)

        #   parent
        self._parent = parent
        #   pointer to MBD system as object attribute
        self.MBD_system = MBD_system

        #    gravity
        self.gravity = MBD_system.gravity * MBD_system.gravity_vector

    def _C_q_size(self):
        """
        Evaluate MBD system matrix C_q size
        :return:
        """
        self.cols_C_q = 3 * self.MBD_system.number_of_bodies
        self.rows_C_q = self.MBD_system.evaluate_C_number_of_rows()

    def evaluate_M(self):
        """
        Create mass matrix
        """
        self.M_dim = 3 * self.MBD_system.number_of_bodies
        self.q_dim = 6 * self.MBD_system.number_of_bodies

        self._C_q_size()

        #    predefine mass matrix size
        M = np.zeros([self.M_dim, self.M_dim])

        for i, body in enumerate(self.MBD_system.bodies):
            #    body mass matrix for 2D motion
            body_matrix = np.array([[body.mass, 0, 0],
                                    [0, body.mass, 0],
                                    [0, 0, body.J_zz]])

            M[3 * i:3 * i + 3, 3 * i:3 * i + 3] = body_matrix

        self.M = M
        self.M_created = True

        #    inverse mass matrix
        self.M_inv = np.linalg.inv(self.M)

    def getDOF(self):
        """
        Evaluate degrees of freedom
        """
        q0 = self.MBD_system.create_q0()
#         print "number of coordinates =", len(self.MBD_system.q0) / 2
        C_q, C_qT = self.create_C_q(q)
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
            C_joint = joint.evaluate_C(q)

            if joint.joint_type == "fixed":
                C[3 * joint.fixed_joint_id:3 * joint.fixed_joint_id + 3] = C_joint
            if joint.joint_type == "revolute":
                C[3 * self.MBD_system.number_of_fixed_joints + 2 * joint.revolute_joint_id:3 * self.MBD_system.number_of_fixed_joints + 2 * joint.revolute_joint_id + 2] = C_joint
            if joint.joint_type == "prismatic":
                C[3 * self.MBD_system.number_of_fixed_joints + 2 * self.MBD_system.number_of_revolute_joints + 2 * joint.prismatic_joint_id:3 * self.MBD_system.number_of_fixed_joints + 2 * self.MBD_system.number_of_revolute_joints + 2 * joint.prismatic_joint_id + 2] = C_joint

        for i, motion in enumerate(self.MBD_system.motions):
            C[self.MBD_system.C_q_number_of_rows + i] = motion.evaluate_C(q, t=t)
        # print "C =", C
        return C

    def check_C(self, q, t):
        """
        Function evaluates C if tolerance is violated
        :param q:
        :param t:
        :return:
        """
        C = self.evaluate_C(q, t)

        #   all elements of vector C have to be under specified tolerance defined by used
        #   and number of elementc of C has to be equal to number of absolute coordinates of MBD system
        if (abs(C) < self.MBD_system.TOL_C).all() and (len(C) == 3*self.MBD_system.number_of_bodies):
            logging.getLogger("DyS_logger").info("C(q, t) check - proceed with simulation run!")
            return False
        else:
            logging.getLogger("DyS_logger").info("C(q, t) check - tolerances violated!")
            return True

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
        #    predefine C_q matrix size

        # #    C_q matrix if there are contacts in compression phase
        # if self._parent.FLAG_contact and any(contact._phase == "C" for contact in self.MBD_system.contacts):
        #
        #     #    count number of active contacts
        #     self.active_contacts = 0
        #     self.active_contacts = sum(1 for contact in self.MBD_system.contacts if contact.contact_detected and contact.contact_distance_inside_tolerance and contact.status == 2 and contact._phase == "C")
        #
        #     self.rows_C_q_contacts = 3 * self.active_contacts
        #
        #     C_q = np.zeros([self.rows_C_q_contacts, self.cols_C_q])
        #
        #
        #     for i, contact in enumerate(self.MBD_system.contacts):
        #
        #         if contact.status == 2 and contact._phase == "C":
        #             #    creates a C_q matrix of a contact based on contact's properties (type, geometry, etc.)
        #             contact_matrix_body_i, contact_matrix_body_j = contact.create_contact_C_q_matrix(t_, q_, self.MBD_system.bodies)
        #
        #             C_q[i:self.rows_C_q + i + 3, 3 * contact.body_id_i:3 * contact.body_id_i + 3] = contact_matrix_body_i.matrix
        #             C_q[i:self.rows_C_q + i + 3, 3 * contact.body_id_j:3 * contact.body_id_j + 3] = contact_matrix_body_j.matrix


        #    C_q matrix if there are no contacts in compression phase and only constraint equations
        self._C_q_size()
        C_q = np.zeros([self.rows_C_q, self.cols_C_q])
        for joint in self.MBD_system.joints:
            #    creates a C_q matrix of a joint based on joint's properties (type, geometry, etc.)
            C_q_i, C_q_j = joint.evaluate_C_q(q)

            #    C_q matrix for joint type: fixed
            if joint.joint_type == "fixed":
                #    add matrix of body j to Cq matrix
                if joint.body_id_i == "ground":
                    C_q[3 * joint.fixed_joint_id:3 * joint.fixed_joint_id + 3, 3 * joint.body_id_j:3 * joint.body_id_j + 3] = C_q_j
                #    add matrix of body i to Cq matrix
                elif joint.body_id_j == "ground":
                    C_q[3 * joint.fixed_joint_id:3 * joint.fixed_joint_id + 3, 3 * joint.body_id_i:3 * joint.body_id_i + 3] = C_q_i
                #    add matrix of body i and j to Cq matrix
                elif (joint.body_id_i != "ground") and (joint.body_id_j != "ground"):
                    C_q[3 * joint.fixed_joint_id:3 * joint.fixed_joint_id + 3, 3 * joint.body_id_i:3 * joint.body_id_i + 3] = C_q_i
                    C_q[3 * joint.fixed_joint_id:3 * joint.fixed_joint_id + 3, 3 * joint.body_id_j:3 * joint.body_id_j + 3] = C_q_j

            #    C_q matrix for joint type: revolute
            elif joint.joint_type == "revolute":
                #    add matrix of body i and j to Cq matrix
                if (joint.body_id_i != "ground") and (joint.body_id_j != "ground"):
                    C_q[3 * self.MBD_system.number_of_fixed_joints + 2 * joint.revolute_joint_id:3 * self.MBD_system.number_of_fixed_joints + 2 * joint.revolute_joint_id + 2, 3 * joint.body_id_i:3 * joint.body_id_i + 3] = C_q_i
                    C_q[3 * self.MBD_system.number_of_fixed_joints + 2 * joint.revolute_joint_id:3 * self.MBD_system.number_of_fixed_joints + 2 * joint.revolute_joint_id + 2, 3 * joint.body_id_j:3 * joint.body_id_j + 3] = C_q_j
                #    add matrix of body j to Cq matrix
                elif joint.body_id_i == "ground":
                    C_q[3 * self.MBD_system.number_of_fixed_joints + 2 * joint.revolute_joint_id:3 * self.MBD_system.number_of_fixed_joints + 2 * joint.revolute_joint_id + 2, 3 * joint.body_id_j:3 * joint.body_id_j + 3] = C_q_j
                #    add matrix of body i to Cq matrix
                elif joint.body_id_j == "ground":
                    C_q[3 * self.MBD_system.number_of_fixed_joints + 2 * joint.revolute_joint_id:3 * self.MBD_system.number_of_fixed_joints + 2 * joint.revolute_joint_id + 2, 3 * joint.body_id_i:3 * joint.body_id_i + 3] = C_q_i

            #    C_q matrix for joint type: prismatic
            elif joint.joint_type == "prismatic":
                #    add matrix of body i and j to Cq matrix
                if (joint.body_id_i != "ground") and (joint.body_id_j != "ground"):
                    C_q[3 * self.MBD_system.number_of_fixed_joints + 2 * self.MBD_system.number_of_revolute_joints + 2 * joint.prismatic_joint_id:3 * self.MBD_system.number_of_fixed_joints + 2 * self.MBD_system.number_of_revolute_joints + 2 * joint.prismatic_joint_id + 2, 3 * joint.body_id_i:3 * joint.body_id_i + 3] = C_q_i
                    C_q[3 * self.MBD_system.number_of_fixed_joints + 2 * self.MBD_system.number_of_revolute_joints + 2 * joint.prismatic_joint_id:3 * self.MBD_system.number_of_fixed_joints + 2 * self.MBD_system.number_of_revolute_joints + 2 * joint.prismatic_joint_id + 2, 3 * joint.body_id_j:3 * joint.body_id_j + 3] = C_q_j
                #    add matrix of body j to Cq matrix
                elif joint.body_id_i == "ground":
                    C_q[3 * self.MBD_system.number_of_fixed_joints + 2 * self.MBD_system.number_of_revolute_joints + 2 * joint.prismatic_joint_id:3 * self.MBD_system.number_of_fixed_joints + 2 * self.MBD_system.number_of_revolute_joints + 2 * joint.prismatic_joint_id + 2, 3 * joint.body_id_j:3 * joint.body_id_j + 3] = C_q_j
                #    add matrix of body i to Cq matrix
                elif joint.body_id_j == "ground":
                    C_q[3 * self.MBD_system.number_of_fixed_joints + 2 * self.MBD_system.number_of_revolute_joints + 2 * joint.prismatic_joint_id:3 * self.MBD_system.number_of_fixed_joints + 2 * self.MBD_system.number_of_revolute_joints + 2 * joint.prismatic_joint_id + 2, 3 * joint.body_id_i:3 * joint.body_id_i + 3] = C_q_i

            else:
                raise Exception, "C_q matrix not constructed. Unknown joint type!"
            
        
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
#         if self._parent.FLAG_contact == 2 and :
        #    create Q_d vector with contact present
        # if any(contact._phase == "C" for contact in self.MBD_system.contacts):
        #     Q_d = np.zeros(self.rows_C_q + self.cols_C_s)
        #     for contact in self.MBD_system.contacts:
        #
        #         if contact.contact_detected and contact.contact_distance_inside_tolerance and contact._phase == "C":
        #             Q_d_contact_vector = contact.create_contact_Q_d_vector(_q=q_)
        #
        #             Q_d = Q_d_contact_vector
        #
        #     if len(Q_d) == self.rows_C_q:
        #         pass
        #     else:
        #         print "Q_d =", Q_d
        #         raise ValueError, "Vector Q_d not correct size."
        #     return Q_d

        #    create Q_d vector as no contact is present
        # else:
        Q_d = np.zeros(self.rows_C_q)

        for joint in self.MBD_system.joints:

            Q_d_joint_vector = joint.evaluate_Q_d(q)

            if joint.joint_type == "fixed":
                Q_d[3 * joint.fixed_joint_id:3 * joint.fixed_joint_id + 3] = Q_d_joint_vector

            if joint.joint_type == "revolute":
                Q_d[3 * self.MBD_system.number_of_fixed_joints + 2 * joint.revolute_joint_id:3 * self.MBD_system.number_of_fixed_joints + 2 * joint.revolute_joint_id + 2] = Q_d_joint_vector

            if joint.joint_type == "prismatic":
                Q_d[3 * self.MBD_system.number_of_fixed_joints + 2 * self.MBD_system.number_of_revolute_joints + 2 * joint.prismatic_joint_id:3 * self.MBD_system.number_of_fixed_joints + 2 * self.MBD_system.number_of_revolute_joints + 2 * joint.prismatic_joint_id + 2] = Q_d_joint_vector

        return Q_d

    def create_Q_e(self, t, q_):
        """
        Evaluate vector of generalized external forces of the system
        """
        # print "create_Q_e()"
        # print "t =", t
        if self._parent.FLAG_contact == 1:
            self.solve_contacts(t, q_)
        # print "********************************************"
        # print "t =", t
        # print "q(in) =", q_
        #    vector of generalised external forces
        Q_e = np.zeros(3 * self.MBD_system.number_of_bodies)

        #    add gravity forces
        if (self.gravity != np.zeros(3)).any():
            for body in self.MBD_system.bodies:
                M_b = np.array([body.mass, body.mass, body.J_zz])
                Q_e_b = M_b * self.gravity
                Q_e[3 * body.body_id:3 * body.body_id + 3] = Q_e_b

        for force in self.MBD_system.forces:
            Q_e_f = force.create_force_Q_e_vector(t, q_)

            Q_e[3 * force.body_id:3 * force.body_id + 3] = Q_e[3 * force.body_id:3 * force.body_id + 3] + Q_e_f

        for spring in self.MBD_system.springs:
            Q_i, Q_j = spring.create_spring_force_Q_e_vector(q=q_)

#             self.Q_list = [Q_i, Q_j]
#             print "---------------------------------------"
#             for body_id in spring.body_id_list:
#                 if body_id == "ground":
#                     None
#                 else:
#                     Q_e[3 * int(body_id):3 * int(body_id) + 3] = Q_e[3 * int(body_id):3 * int(body_id) + 3] + self.Q_list[spring.body_id_list.index(body_id)]
#         print "*********************************"
            if spring.body_id_i != "ground":
                Q_e[3 * spring.body_id_i:3 * spring.body_id_i + 3] = Q_e[3 * spring.body_id_i:3 * spring.body_id_i + 3] + Q_i
            if spring.body_id_j != "ground":
                Q_e[3 * spring.body_id_j:3 * spring.body_id_j + 3] = Q_e[3 * spring.body_id_j:3 * spring.body_id_j + 3] + Q_j
        # print "Q_e =", Q_e
        return Q_e

    def create_C_s(self, q):
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
#                 print "contact_matrix_body_i ="
#                 print contact_matrix_body_i.matrix
#
#                 print "contact_matrix_body_j ="
#                 print contact_matrix_body_j.matrix
#                 print "==============================================="
                C_s[i:i + 3, 2 * i] = contact_matrix_body_i.matrix[:, 0]
                C_s[i:i + 3, 2 * i + 1] = contact_matrix_body_j.matrix[:, 0]

#                 C_s[self.rows_C_q + i:self.rows_C_q + i + 3, 2 * i] = contact_matrix_body_i.matrix[:, 0]
#                 C_s[self.rows_C_q + i:self.rows_C_q + i + 3, 2 * i + 1] = contact_matrix_body_j.matrix[:, 0]


        #    size of C_q matrix
        [self.rows_C_s, self.cols_C_s] = C_s.shape

        #    matrix C_q transposed
        C_s_trans = C_s.T
        return C_s, C_s_trans

    def create_dC(self, C_q, q):
        """

        """
        dq = q[0.5 * self.q_dim:]
        dC = np.dot(C_q, dq)
        return dC

    def evaluate_dq(self, h, t, q):
        """
        Function evaluates vector dq of MBD system
        """
        #    construct Cq matrix
        self.C_q = self.evaluate_C_q(q)
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

        #    vector of external and contact forces
        Q_e = self.create_Q_e(t, q)
        _vector[0:self.M_dim] = Q_e
        # print "%2.10f"%t, Q_e, "FLAG =", self._parent.FLAG_contact

        #    vector of elements, that are quadratic in velocities
        if self.MBD_system.use_BSM:
            C = self.evaluate_C(q, t)
            dC = self.create_dC(self.C_q, q)

            self.alpha = h#1./h
            self.beta = self.alpha #* np.sqrt(2)

            Q_d = self.evaluate_Q_d(q) - 2*self.alpha*dC - (self.beta**2)*C
        else:
            Q_d = self.evaluate_Q_d(q)

        _vector[self.M_dim:rows] = Q_d

        B = np.zeros([self.col_C_q, self.rows_C_q])
        B = self.C_qT

        C = np.zeros([self.rows_C_q, self.col_C_q])
        C = self.C_q

        D = np.zeros([self.rows_C_q, self.rows_C_q])

        #    inverse matrix block-wise
        _matrix_inverse = inverse_blockwise.inverse_blockwise(self.M_inv, B, C, D, self.M_dim)

        #    solve system of equations
        _sol_vector = np.dot(_matrix_inverse, _vector)

        #    accelerations vector
        ddq = _sol_vector[0:self.M_dim]

        #    vector of lagrange multipliers
        L = _sol_vector[self.M_dim + self.cols_C_s:rows]

        #    generalized reaction forces associated with C_q matrix
        Q_c = -np.dot(self.C_qT, L)

        #    the dq vector in constructed from vector ddq and the second part of the input vector q
        dq = np.append(q[0.5 * self.q_dim:], ddq)
        return dq

    def solve_contacts(self, t, q):
        """
        Solve all contact that are detected at integration time t
        """
        #   loop through all contacts
        # print "q =", q
        for contact in self.MBD_system.contacts:
            # print "contact.contact_detected =", contact.contact_detected
            # print "contact.contact_distance_inside_tolerance =", contact.contact_distance_inside_tolerance
            #   check if contact is detected and calculated contact distance is inside the tolerance, then contact is present
            if contact.contact_detected:#and contact.contact_distance_inside_tolerance
                #   get contact forces
                contact.solve(t, q)

                # if not contact.list_of_contact_force_objects_constructed:
                #     for body_id, u_P in zip(contact.body_id_list, contact.u_P_list):
                #         _force = Force(force_name=contact._name + "_on_body_" + str(body_id), body_id=body_id, Fx=contact_force[0], Fy=contact_force[1], Mz=0, u_iP_f=u_P)
                #
                #         #   add pair of contact forces to forces list of MBD system
                #         self.MBD_system.forces.append(_force)
                #         #   add pair of contact forces to forces list of contact
                #         contact.contact_force_list.append(_force)
                #
                #         # contact.contact_force_list.append(_force)
                #     contact.list_of_contact_force_objects_constructed = True
                    
                # else:
                #     for _force in contact.contact_force_list:
                #         _force.update_force_vector(Fx=contact_force[0], Fy=contact_force[1])
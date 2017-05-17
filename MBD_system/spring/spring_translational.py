"""
Created on 21. feb. 2014

@author: lskrinjar
"""
import itertools
from pprint import pprint
import numpy as np

from MBD_system.Ai_ui_P import Ai_ui_P_vector
from MBD_system.dr_ij_P_dq import dr_ij_P_dq
from MBD_system.force.force_vector import Force_Q_e_vector
from MBD_system.q2R_i import q2R_i
from MBD_system.q2dR_i import q2dR_i
from MBD_system.q2dtheta_i import q2dtheta_i
from MBD_system.q2theta_i import q2theta_i
from MBD_system.r_ij_P import r_ij_P
from MBD_system.u_P_cad2cm_lcs import u_P_cad2cm_lcs
from MBD_system.u_P_lcs2gcs import u_P_lcs2gcs
from simulation_control_widget.vtk_widget.marker.marker import Marker
from MBD_system.spring.spring import Spring
from global_variables import GlobalVariables


class SpringTranslational(Spring):
    """
    classdocs
    """
    def __init__(self, _name, body_id_i, body_id_j, u_iP_CAD=np.zeros(2, dtype=float), u_jP_CAD=np.zeros(2, dtype=float), u_iP_LCS=np.zeros(2, dtype=float), u_jP_LCS=np.zeros(2, dtype=float), k=0., c=0., l_0=None, properties_dict={}, parent=None):
        """

        :return:
        """
        super(SpringTranslational, self).__init__(_name, body_id_i, body_id_j, u_iP_CAD=u_iP_CAD, u_jP_CAD=u_jP_CAD, k=k, c=c, properties_dict=properties_dict, parent=parent)

        #    spring type
        self.spring_type = "translational"

        #   undeformed length
        self.l_0 = l_0

        #   add additional properties
        self.add_attributes_from_dict(self.properties)

        #   list of markers
        self.markers = self._create_markers()

        #   create pair of forces list
        self._Fn_list = self._create_Fn_forces()

        #   create list of force objects
        self.Q_e_list = self._create_Q_e()

    def evaluate_Q_e(self, q):
        """
        Function evaluates a generalized external force vector
        :param q:   vector of MBD system
        :return:
        """
        #   predefined list for spring force matrix
        self.spring_force_matrix_list = self._evaluate_r_ij_P_q_list(q)

        #    vector between points of spring element
        self.r_ij_P = self.evaluate_rijP(q)

        #    length of vector - spring length
        self.l = np.linalg.norm(self.r_ij_P, ord=2)

        #    unit vector in direction of vector r_ij_P
        self.I_r_ij_P = self._evaluate_I_r_ij_P()

        #    velocity of deformation of spring length
        self.dl = self._evaluate_dl(q, self.I_r_ij_P)

        #    force value (amplitude) of spring element
        self.f_s = self._evaluate_F(self.l, self.dl)

        #    generalized force on each body that are connected with a spring
        self.Q_e_list = []

        for body_id, force_matrix in zip(self.body_id_list, self.spring_force_matrix_list):
            Q_e = Force_Q_e_vector(self.f_s, force_matrix.matrix, self.I_r_ij_P)
            #    append spring force vector object to list
            self.Q_e_list.append(Q_e)

        [Q_i, Q_j] = self.Q_e_list

        #   update force values in object (for data tracking and visualization)
        self._update_F(Q_e_list=self.Q_e_list)

        return Q_i.Q_e, Q_j.Q_e

    def _evaluate_dl(self, q, I_r_ij_P):
        """

        :return:
        """
        #   vector of dq_i, dq_j
        dq = np.array([np.append(q2dR_i(q, self.body_id_i), q2dtheta_i(q, self.body_id_i)),
                       np.append(q2dR_i(q, self.body_id_j), q2dtheta_i(q, self.body_id_j))]).flatten()

        #    dr_ij_P_dq
        dr_ij_P_dq_ = dr_ij_P_dq(body_id_i=self.body_id_i,
                                 theta_i=q2theta_i(q, self.body_id_i),
                                 u_iP=self.u_iP_LCS,
                                 body_id_j=self.body_id_j,
                                 theta_j=q2theta_i(q, self.body_id_j),
                                 u_jP=self.u_jP_LCS)

        #    velocity of deformation of spring length
        dl = np.dot(I_r_ij_P, np.dot(dr_ij_P_dq_, dq))

        return dl

    def _evaluate_I_r_ij_P(self):
        """

        :return:
        """
        if self.I_r_ij_0 is None:
            I_r_ij_P = self.r_ij_P / self.l
        else:
            I_r_ij_P = self.I_r_ij_0

        return I_r_ij_P

    def evaluate_L(self, q):
        """

        :param q:
        :return:
        """
        self.r_ij_P = self.evaluate_rijP(q)

        l = np.linalg.norm(self.r_ij_P, ord=2)

        if self.l_0 is None:
            self.l_0 = l

        return l
    
    def evaluate_l0l(self, q):
        """
        Function evaluates deformation of spring length
        """
        #    length of vector - spring length
        self.l = self.evaluate_rijP(q)
        
        self.l0l = self.l - self.l_0
        print self.l0l

    def evaluate_potential_energy(self, q):
        """
        Function evaluates mechanical energy of linear spring element
        :param q:
        :return:
        """
        #    length of vector - spring length
        self.l = self.evaluate_L(q)
        
        self.energy = 0.5 * self.k * (self.l_0 - self.l)**2

        return self.energy

    def evaluate_Q_q(self, q):
        """

        :param q:
        :return Q_q_list:
        """
        #   if list is empty evaluate elements
        if not self.spring_force_matrix_list:
            self.spring_force_matrix_list = self._evaluate_r_ij_P_q_list(q)

        if self.l is None:
            self.l = self.evaluate_L(q)

        if self.I_r_ij_P is None:
            self.I_r_ij_P = self._evaluate_I_r_ij_P()

        if self.dl is None:
            self.dl = self._evaluate_dl(q, self.I_r_ij_P)

        if self.f_s is None:
            self.f_s = self._evaluate_F(self.dl)

        #    partial derivative of generalized spring force vector with respect to generalized coordinates
        Q_q_list = []
        for body_id, r_ij_q in zip(self.body_id_list, self.spring_force_matrix_list):
            Q_q_i = np.zeros([3, 3], dtype=float)

            print self.f_s / self.l * r_ij_q.matrix




            Q_q_list.append(Q_q_i)

        return Q_q_list

    def evaluate_Q_dq(self, q):
        """

        :return:
        """

    def _evaluate_f_s_l_q(self, q):
        """
        Evaluate a partial derivative of expression (f_s / l) with respect to generalized coordinates q
        :return:
        """

    def print_F(self):
        """

        :return:
        """
        if self.l is None:
            self.l = self.evaluate_L(self.q0)

        if self.dl is None:
            self.dl = 0.
        print "F ="
        print self._evaluate_F(self.l, self.dl)


if __name__ == "__main__":
    GlobalVariables.q_i_dim = [3, 3]
    spring = SpringTranslational("test", 0, 1, u_iP_CAD=np.array([-1., 0.]), u_jP_CAD=np.array([1., 0.]))
    spring.k = 1.
    spring.l_0 = 2.
    q = np.array([0., 0., 0., 1., 1., 0., 0., 0., 0., 0., 0., 0.])

    spring.evaluate_Q_q(q)
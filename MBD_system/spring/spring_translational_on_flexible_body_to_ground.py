"""
Created on 21. feb. 2014

@author: lskrinjar
"""
import itertools
from pprint import pprint
import numpy as np

from MBD_system.Ai_ui_P import Ai_ui_P_vector
from MBD_system.MBD_system_items import SpringItem
from MBD_system.dr_ij_P_dq import dr_ij_P_dq
from MBD_system.force.force import Force
from MBD_system.force.force_matrix import Force_Q_e_matrix
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


class SpringTranslationalOnFlexibleBodyToGround(Spring):
    """
    classdocs
    """
    def __init__(self, name, body_id_i, node_id=None, element_id=None, element_ksi=None, u_jP_CAD=np.zeros(2), L_i=None, k=0, c=0, l_0=None, properties_dict={}, parent=None):
        """

        :return:
        """
        super(SpringTranslationalOnFlexibleBodyToGround, self).__init__(name, body_id_i, "ground", u_iP_CAD=None, u_jP_CAD=u_jP_CAD, k=k, c=c, properties_dict=properties_dict, parent=parent)

        #    spring type
        self.spring_type = "translational on flexible body to ground"

        #   spring properties
        self.l_0 = l_0

        #   flexible body data
        self.node_id = node_id
        self.element_id = element_id
        self.element_ksi = element_ksi
        self.L_i = L_i

        #   update list
        self.node_id_list[0] = self.node_id

        #   add additional properties
        self.add_attributes_from_dict(self.properties)

        #   flexible body i
        self.body_i = self._parent._parent.bodies[self.body_id_i]

        self.element_ksi, self.element_id, self.node_id = self._parent._parent.bodies[self.body_id_i].evaluate_mesh_data_at_position(self.L_i)

        #   pointer to element
        self._element = self.body_i.mesh.elements[self.element_id]

        #   list of markers
        self.markers = self._create_markers()

        #   create pair of forces list
        self._Fn_list = self._create_Fn_forces()

        #   length
        #   get q of MBD system
        q = self._parent._parent.evaluate_q()
        rijP = self.evaluate_rijP(q)
        self.l = np.linalg.norm(rijP)
        self.dl = 0.

        #   create list of force objects
        self.Q_e_list = self._create_Q_e()

    def add_attributes_from_dict(self, dict):
        """
        Function adds attributes to object from dictionary
        :param dict:
        :return:
        """
        for key in dict:
            val = dict[key]
            if hasattr(self, key):
                setattr(self, key, val)

    def _create_markers(self):
        """
        Function create markers
        :return:
        """
        markers = []

        #   get q of MBD system
        q = self._parent._parent.evaluate_q()

        #   marker on flexible body i
        rP_i = self._parent._parent.bodies[self.body_id_i].evaluate_r(q, element_id=self.element_id, ksi=self.element_ksi)

        #   marker on ground body j
        self.rP_j = self.u_jP_CAD
        self.r_P_list = [rP_i, self.rP_j]

        for body_id, rP in zip(self.body_id_list, self.r_P_list):
            marker = Marker(rP, body_id=body_id, parent=None)
            markers.append(marker)

        return markers

    def _create_Fn_forces(self):
        """

        :return:
        """
        Fn_list = []
        for body_id, rP in zip(self.body_id_list, self.r_P_list):
            #   create force object
            force = Force(body_id, force_name=self._name + "_on_body_" + str(body_id), r_P_GCS=rP, parent=self)

            #   add pair of contact forces to forces list of MBD system
            # self._parent._parent.forces.append(force)

            #   add pair of contact forces to forces list of spring
            Fn_list.append(force)

        return Fn_list

    def evaluate_Q_e(self, q):
        """
        Function evaluates a generalized external force vector
        :param q:   vector of MBD system
        :param q:
        :return:
        """
        #   positions of spring endpoints in GCS
        #   distance vector
        r_ij_P = self.evaluate_rijP(q)

        #    length of vector - spring length
        self.l = np.linalg.norm(r_ij_P, ord=2)

        I_r_ij = self._evaluate_I_r_ij(r_ij_P, self.l)

        #    velocity of deformation of spring length
        dq_ = self.body_i.evaluate_dr(q, element_id=self.element_id, ksi=self.element_ksi)
        dl = np.dot(I_r_ij, dq_)

        #    force value (amplitude) of spring element
        self.F_s = self._evaluate_F(self.l, self.dl)
        if self.direction == "compression":
            if self.F_s > 0.:
                self.F_s = 0.

        if self.direction == "tension":
            if self.F_s < 0.:
                self.F_s = 0.

        F = -self.F_s * I_r_ij

        #   force on flexible body
        S = self._element._evaluate_S(self.element_ksi)
        Q_e_i_element = np.dot(S.T, F)

        Q_e_i = reduce(np.dot, [self._element.B.T, self._element.T.T, Q_e_i_element])

        #   force on ground body
        Q_e_j = None

        return Q_e_i, Q_e_j

    def evaluate_rijP(self, q):
        """
        Length of spring at time when vector of absolute coordinates is equal to q
        """
        rP_i = self._parent._parent.bodies[self.body_id_i].evaluate_r(q, element_id=self.element_id, ksi=self.element_ksi)

        self.r_P_list[0] = rP_i

        #   distance vector
        r_ij_P = rP_i - self.rP_j

        return r_ij_P

    def evaluate_F(self, **kwargs):
        """

        :param q:
        :return:
        """
        if "q" not in kwargs:
            q = self._parent._parent.evaluate_q()

        self.f_s = self._evaluate_F(self.l, self.dl)

        r_ij_P = self.evaluate_rijP(q)
        l = np.linalg.norm(r_ij_P, ord=2)

        I_r_ij = self._evaluate_I_r_ij(r_ij_P, l)

        F = self.f_s * I_r_ij

        return F

    def _evaluate_I_r_ij(self, r_ij_P, l):
        """

        :return:
        """
        #    unit vector in direction of vector r_ij_P
        if self.I_r_ij_0 is None:
            I_r_ij = r_ij_P / l
        else:
            I_r_ij = self.I_r_ij_0

        return I_r_ij

    def reset(self, q):
        """

        :return:
        """
        self.r_P_list[0] = self._parent._parent.bodies[self.body_id_i].evaluate_r(q, element_id=self.element_id, ksi=self.element_ksi)


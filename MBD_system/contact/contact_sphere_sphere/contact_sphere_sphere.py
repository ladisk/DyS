__author__ = 'lskrinjar'

import time
from pprint import pprint
import numpy as np


from MBD_system.contact.contact import Contact
from MBD_system.q2R_i import q2R_i
from MBD_system.contact.distance.distance import Distance
from MBD_system.A import A_matrix
from MBD_system.dr_contact_point_uP import dr_contact_point_uP
from MBD_system.q2dR_i import q2dR_i
from MBD_system.q2theta_i import q2theta_i
from MBD_system.q2dtheta_i import q2dtheta_i
from MBD_system.transform_cs import gcs2cm_lcs
from MBD_system.transform_cs import uP_gcs2lcs
from MBD_system.n2t import n2t


class ContactSphereSphere(Contact):
    """
    classdocs
    """


    def __init__(self, _type, body_id_i, body_id_j, R0_i, R0_j, properties_dict={}, parent=None):
        """
        
        """
        #    this has to be after attributes contact_id and _name as name is constructed from contact_id and _name
        super(ContactSphereSphere, self).__init__(_type, body_id_i, body_id_j, name=None, properties_dict=properties_dict, parent=parent)

        #    name as string
        self._name = "Contact_Sphere_Sphere_" + str(self.contact_id)

        #    contact type
        self._type = "Contact Sphere-Sphere"
        
        #    radius of spheres
        self.R0_i = R0_i
        self.R0_j = R0_j
        self.R0_list = [self.R0_i, self.R0_j]
        
        #    minimum distance value at contact
        self._distance0 = self.R0_i + self.R0_j

        #   contact model
        self._create_contact_model(self.properties)

        #   friction model
        self._create_friction_model(self.properties)

        #   list of markers
        self.markers = self._create_markers()

    def _create_markers(self):
        """

        :return:
        """
        return []

    def check_for_contact_started_condition(self, delta, sign):
        """

        :return:
        """
        if (sign == -1) and (delta <= 0):
            return True
        else:
            return False

    def check_for_contact_continued_condition(self, delta, dq_n):
        """

        :return:
        """
        return (delta <= self._delta0)

    def check_for_contact(self, step, t, q):
        """
        Function check for contact between two spheres
        returns:
                -1 - first contact indentation is too large, return to previous time step and continue with smaller step size
                0 - no contact
                +1 - contact detected
        """
        # print "check_for_contact()"
        #   step
        self._step = step
        #   time
        self._t = t

        #   evaluate normal
        # self._n = self._distance_obj.get_normal_2D()
        self._distance, self._delta, self._n_GCS, self._t_GCS = self._contact_geometry_GCS(q)
        # print "self._delta =", self._delta
        #   calculate current contact geometry
        self._get_contact_geometry_data(q)

        #   calculate relative contact velocity in normal direction at time t
        # self._dq_n, self._dq_t = self._contact_velocity(q)
        # print "self._dq_n =", self._dq_n
        #   add calculated data to container
        # if self._step_solution_accepted:
        #     self._distance_solution_container = np.append(self._distance_solution_container, self._delta)  # _distance_sign, _distance
        # else:
        #     self._distance_solution_container[-1] = self._delta
        
        #   check sign of product between last to values of distance
        self._sign_check = np.sign(self._distance_solution_container[-1] * self._delta)
        # print "self._distance =", self._delta
        #    contact
        if self.check_for_contact_started_condition(self._delta, self._sign_check):#(self._delta <= 0) and (self._sign_check == -1):# and self._dq_n < 0:
            #   condition to find first contact
            if  abs(self._delta) <= self.distance_TOL:
                self._delta0 = self._delta
                self.contact_detected = True
                self.contact_distance_inside_tolerance = True
                self.status = 1
                print "self._n_GCS =", self._n_GCS
                # print "contact detected! @check_for_contact()"
                # print "step =", self._step
                # print "self._delta =", self._delta
                # self._status_container = np.append(self._status_container, self.status)
            #   condition to check if penetration depth of first contact is too large
            else:# self._sign_check == -1 and self._distance > -self.distance_TOL:
                self.status = -1
            # #   no contact
            # else:
            #     self.contact_distance_inside_tolerance = False
            #     self.status = 1
            #     self._status_container = np.append(self._status_container, self.status)
        #    no contact
        else:
            self.contact_detected = False
            self.contact_distance_inside_tolerance = False
            self.status = 0
            # self._status_container = np.append(self._status_container, self.status)
            self.no_contact()

        return self.status

    def contact_update(self, step, t, q):
        """
        Function updates status only for contact type - 2 speheres contact
        :return:
            status - status of contact 0-no contact, 2-contact
        """
        # print "contact_update()"
        #   step value
        self._step = step

        #   calculate relative contact velocity in normal direction at time t
        self._dq_n, self._dq_t = self._contact_velocity(q)

        #   calculate distance: node-edge
        # self._distance, _inside = evaluate_distance_2D(self.node_GCS, self.edge_GCS[0], self._n, self._t)
        self._delta = self._evaluate_contact_distance(q)

        if self._delta <= self._delta0:#self._delta0:# and self._dq_n <= abs(self._dq0_n):#_inside and self._distance < -1*self.distance_TOL:
            self.status = 1
        else:
            # print "contact OVER@contact_update()"
            self.status = 0
            self._contact_point_found = False
            self.initial_contact_velocity_calculated = False
            self.contact_distance_inside_tolerance = False
            self.contact_detected = False
            self._status_container = np.append(self._status_container, self.status)

            #   reset values of contact force object to zero if no contact
            _F0 = np.zeros(2)
            _u_P0 = np.zeros(2)
            for force_n, force_t in zip(self._Fn_list, self._Ft_list):
                force_n.update(self._step, F=_F0, u_P=_u_P0)
                force_t.update(self._step, F=_F0, u_P=_u_P0)
        return self.status

    def _get_contact_geometry_data(self, q):
        """
        Function calculates normal and tangent vector of body i, j in GCS
        """
        #   normal in GCS of body i, j
        self._n_GCS_list = [+self._n_GCS, -self._n_GCS]

        #   tangent in GCS of body i, j
        self._t_GCS_list = [+self._t_GCS, -self._t_GCS]

    def _contact_velocity(self, q):
        """
        Function evaluates relative contact velocity vectors in normal and tangent direction
        """
        dr_P = []
        for body_id, u_P in zip(self.body_id_list, self.u_P_LCS_list):
            #   body velocity, R, theta
            _dR = q2dR_i(q, body_id)
            _theta = q2theta_i(q, body_id)

            #    dtheta - omega
            _dtheta = q2dtheta_i(q, body_id)

            #    velocity at contact point (in GCS) on a body
            _dr_P = dr_contact_point_uP(_dR, _theta, _dtheta, u_P)
            
            #    appned to list
            dr_P.append(_dr_P)
        
        #    relative contact velocity vector
        _dq = dr_P[1] - dr_P[0]
        
        #   relative contact velocity (scalar value)
        #   normal direction
        _dq_n = np.dot(_dq, self._n_GCS)

        #   tangent direction
        _dq_t = np.dot(_dq, self._t_GCS)
        return _dq_n, _dq_t

    def _contact_geometry_GCS(self, q):
        """
        Function calculates distance between center points os spheres and indentation based on radius of sphere i and j
        Evaluated parameters:
        :param _distance:   distance between centers of sphere i and j
        :param _delta:      penetration depth at contact point
        """
        #    create distance object
        # print "q =", q
        # print "i =", self.body_id_i
        # print "j =", self.body_id_j
        # print "i =", q2R_i(q, self.body_id_i)
        # print "j =", q2R_i(q, self.body_id_j)
        self._distance_obj = Distance(q2R_i(q, self.body_id_i), q2R_i(q, self.body_id_j))

        #   distance
        _d = self._distance_obj._distance

        #   delta
        _delta = self._distance_obj._distance - self._distance0

        #   normal in GCS
        _n_GCS = -self._distance_obj._distance_vector / self._distance_obj._distance

        #   tangent in GCS
        _t_GCS = n2t(_n_GCS)

        return _d, _delta, _n_GCS, _t_GCS

    def _evaluate_contact_distance(self, q):
        """
        Function evaluates contact distance for this type of contact
        :return:
        """
        #    create distance object
        self._distance_obj = Distance(q2R_i(q, self.body_id_i), q2R_i(q, self.body_id_j))

        #    distance value
        _distance = self._distance_obj._distance - self._distance0
        return _distance

    def _contact_geometry_LCS(self, q):
        """
        Function calculates the contact geometry from GCS to LCS
        Especially uPi, uPj
        :param q:
        :return:
        """
        #   evaluate contact points on each body in body LCS
        for i, (body_id, n_i, R0_i) in enumerate(zip(self.body_id_list, self._n_GCS_list, self.R0_list)):
            #   R
            R_i = q2R_i(q, body_id)

            #   theta
            theta_i = q2theta_i(q, body_id)

            #   normal in LCS
            self._n_LCS_list[i] = uP_gcs2lcs(u_P=-n_i, theta=theta_i)

            #   tangent in LCS
            self._t_LCS_list[i] = n2t(self._n_LCS_list[i])

            #   calculate actual contact point in LCS on a body surface
            self.u_P_LCS_list[i] = R0_i * self._n_LCS_list[i]

            #   calculate contact point on each body in GCS
            self.r_P_GCS_list[i] = R_i + uP_gcs2lcs(self.u_P_LCS_list[i], theta_i)

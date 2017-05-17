"""
Created on 27. jan. 2016

@author: lskrinjar
"""
import time
from __builtin__ import enumerate
from operator import attrgetter
from pprint import pprint
import itertools
import numpy as np
import matplotlib.pyplot as plt


from MBD_system.contact.contact import Contact
from MBD_system.q2R_i import q2R_i
from MBD_system.q2theta_i import q2theta_i
from MBD_system.transform_cs import gcs2cm_lcs
from MBD_system.Ai_ui_P import Ai_ui_P_vector
from MBD_system.q2dR_i import q2dR_i
from MBD_system.q2dtheta_i import q2dtheta_i
from MBD_system.dr_contact_point_uP import dr_contact_point_uP
from MBD_system.contact.distance.contact_point import ContactPoint


class GeneralContact(Contact):
    """
    classdocs
    """
    __id = itertools.count(0)

    def __init__(self, _type, body_id_i, body_id_j, u_iP, u_jP, properties_dict={}, parent=None):
        """
        Constructor of class contact of revolute clearance joint
        :param _type:       type of clearance joint or contact
        :param body_id_i:   id of hole body
        :param body_id_j:   id of pin body
        :param properties_dict: additioanl parameters to override default values or add new parameters
        """
        super(GeneralContact, self).__init__( _type, body_id_i, body_id_j, properties_dict=properties_dict, parent=parent)
        #   parent
        self._parent = parent

        #    type of contact
        self._contact_type = "general"

        #   vector of axis on revolute joint in LCS of a body i, j
        self.u_iP = u_iP
        self.u_jP = u_jP
        #   predefined empty list of center point or clearance joint (to center of pin/hole) in body LCS
        self.u_CP_LCS_list = [self.u_iP, self.u_jP]

        #   centers of revolute clearance joint in GCS
        self.u_iCP_GCS = np.zeros(2)
        self.u_jCP_GCS = np.zeros(2)
        self.r_CP_GCS_list = [self.u_iCP_GCS, self.u_jCP_GCS]

        #    joint body ids
        self.body_id_i = body_id_i
        self.body_id_j = body_id_j
        self.body_id_list = [self.body_id_i, self.body_id_j]

        #   contact model
        self._create_contact_model(self.properties)

        #   friction model
        self._create_friction_model(self.properties)

        #   list of markers
        self.markers = self._create_markers()

        #   AABB - Axis Aligned Bounding Box properties
        #   initialize a list of AABB
        self.AABB_i = None
        self.AABB_j = None
        self.AABB_list = []
        self.roughness_profile_list = []
        #   initialize a list of contact pairs
        self.AABB_list_of_overlap_pairs = []
        #   contact status attributes
        self.AABB_overlap_detected = False

        #   reset to initial value
        self.reset()

    def set_vtk_data(self, interactor=None):
        """

        :return:
        """
        # print "exe in contact_general.py"
        for AABB, contact_profile in zip(self.AABB_list, self.roughness_profile_list):
            # pprint(vars(contact_profile.geometry))
            AABB.set_vtk_data(interactor, contact_profile.geometry.polygon)#points, actor, polygon

    def contact_velocity(self, q):
        """
        Function evaluates relative contact velocity vectors in normal and tangent direction
        :param q:
        :return:
        """
        dq_n = []
        dq_t = []
        for contact_point in self._contact_point_obj_list:
            _dq_n, _dq_t = contact_point._contact_velocity(q)
            dq_n.append(_dq_n)
            dq_t.append(_dq_t)

        return dq_n, dq_t

    def _no_contact(self):
        """

        :return:
        """
        self.body_id_edge = None
        self.body_id_node = None

        for contact_point in self._contact_point_obj_list:
            if not contact_point.active:
                contact_point.deactivate_forces()

        for contact_point in self._contact_point_obj_list:
            if not contact_point.active:
                self._contact_point_obj_list.remove(contact_point)
                # print "contact point removed!"
            # for Fn, Ft in zip(contact_point._Fn_list, contact_point._Ft_list):
            #     Fn._update_F(np.zeros(2, dtype=float))#step=self._step)
            #     Ft._update_F(np.zeros(2, dtype=float))#step=self._step)
                #
                # del Fn
                # del Ft

        self._distance_obj_list = []
        # self._contact_point_obj_list = []

    def _solution_containers_additional(self):
        """
        Function creates additional solution containers that are specific for type of contact
        :return:
        """
        #   eccentricity vector
        self._e_solution_container = [0, 0]

    def _create_markers(self):
        """
        Function creates markers
        :return:
        """
        return []

    def check_for_contact_started_condition(self, delta, sign):
        """
        Function checks condition for contact
        """
        # print delta
        # print sign
        if type(delta) is not list:
            if (sign == -1) and (delta <= 0.):
                return True
            else:
                return False
        else:
            status_list = []
            for _delta, _sign in zip(delta, sign):
                # print '_sign =', _sign, "_delta =", _delta
                if (_sign == -1) and (_delta <= 0.):
                    status_list.append(True)
                else:
                    status_list.append(False)
            # print "status_list =", status_list
            if any(status_list) or all(status_list):
                return True
            else:
                return False

    def check_for_contact(self, step, t, q):
        """
        Function check for contact between contact pair of AABBs
        returns:
                -1 - first contact indentation is too large, return to previous time step and continue with smaller step size
                0 - no contact
                +1 - contact detected
        """
        # print "check_for_contact(), step =", step
        self._step = step

        #   if distance object list is empty, get distance objects from pair of overlapping AABBs
        #   and when distanc elist object is already created
        if self._distance_obj_list==[]:
            #    loop through list of overlapped AABB object pairs
            for i, overlap_pair_obj in enumerate(self.AABB_list_of_overlap_pairs):
                #   check for contact of overlapped AABB pair (object)
                #   list of distance objects is returned from every pair of AABBs
                _distance_list = overlap_pair_obj.check_for_contact_2D(q)
                #   extend list with list
                self._distance_obj_list.extend(_distance_list)

            for i, _distance_obj in enumerate(self._distance_obj_list):
                if _distance_obj.inside():
                    _distance_obj.geometry_LCS(q)

        #   distance object has already been created and only needs to be updated
        else:
            for i, _distance_obj in enumerate(self._distance_obj_list):
                if _distance_obj._inside:
                    _distance_obj.update_contact_geometry_GCS(q)

        #   contact geometry  in GCS
        self._distance, self._delta, self._n_GCS, self._t_GCS = self._contact_geometry_GCS(q)

        if type(self._distance_solution_container[-1]) is float:
            _distance_solution_container_list = self._distance_solution_container[-1] * np.ones(len(self._delta))
        else:
            _distance_solution_container_list = self._distance_solution_container[-1]

        #   list of sign checks
        self._sign_check = np.sign(np.multiply(np.array(self._delta), _distance_solution_container_list))
        # print "self._sign_check =", self._sign_check
        # print "self._delta =", self._delta
        if self.check_for_contact_started_condition(self._delta, self._sign_check):
            if np.any(abs(np.array(self._delta)) <= self.distance_TOL):
                # print "contact detected, step =", step
                self._delta0 = self._delta
                self.contact_detected = True
                self.contact_distance_inside_tolerance = True
                self._contact_point_found = True
                self.status = 1

                #   body ids
                # fig = plt.figure(num=1, figsize=(6, 4), dpi=100, facecolor='w', edgecolor='k')
                # ax = plt.subplot(111, aspect="equal")

                #   remove duplicate objects
                self._distance_obj_list = list(set(self._distance_obj_list))
                for i, _distance_obj in enumerate(self._distance_obj_list):
                    # _distance_obj.save_plot()
                    if abs(_distance_obj._distance_sign) <= self.distance_TOL:
                        # print "i =", i, "id =", _distance_obj.id, _distance_obj
                        contact_point_obj = ContactPoint(_distance_obj, parent=self)
                        contact_point_obj.evaluate_contact_points(q)
                        contact_point_obj._contact_point_found = True
                        contact_point_obj._delta0 = contact_point_obj._distance_sign

                        #   check and if contact point object with value of r_iP is already in list it does not append it
                        if next((_contact_point_obj for _contact_point_obj in self._contact_point_obj_list if (_contact_point_obj.r_iP == contact_point_obj.r_iP).all()), None) is None:
                            self._contact_point_obj_list.append(contact_point_obj)
                            # contact_point_obj.plot()
                # plt.xlim([-0.025, +0.009])
                # plt.ylim([-0.009, +0.009])
                # filename = "step_" + str(step).zfill(4) + "_contact_started" + ".png"
                # plt.savefig(filename)
                # print "Plot saved to file: ", filename
                # plt.clf()

                #   clear distance object list when contact point object list is created
                delattr(self, "_distance_obj_list")

                print "contact detected!", "self._step =", self._step#, self._contact_point_found#, self.F

                return self.status

            # step back
            if np.any(abs(np.array(self._delta)) > self.distance_TOL):
                if self._step > 66:
                    print "self._delta ="
                    print self._delta
                self.contact_detected = True
                self.status = -1
                return self.status

        # all calculated distances are greater than tolerance and bodies are not in contact
        self.status = 0
        self.no_contact()

        return self.status

    def check_for_contact_continued_condition(self, delta, dq_n, q):
        """
        Function checks condition for contact if contact is still present
        """
        TOL = 0.
        # print 'step=', self._step
        if type(delta) is float:
            return [(delta <= TOL)] # self._delta0and (self._dq_n < abs(self._dq0_n))
        else:
            _delta = []
            for contact_point in self._contact_point_obj_list:

                contact_point.update_contact_geometry_GCS(q)
                # contact_point._evaluate_distance()

                # print "fi1, fi2 =", contact_point.fi1, np.rad2deg(contact_point.fi1), contact_point.fi2, np.rad2deg(contact_point.fi2)
                # print contact_point.fi1 <= np.pi/2.
                # pprint(vars(contact_point))
                # print 'contact_point._distance_sign=', contact_point._distance_sign
                if contact_point._distance_sign <= contact_point._delta0:
                    _delta.append(True)
                else:
                    contact_point._contact_point_found = False
                    contact_point.active = False
                    _delta.append(False)

            if any(_delta):
                return True
            else:
                return False

    def contact_update(self, step, t, q):
        """
        Function evaluates contact distance while, contact is present
        :return:
            status - status of contact 0-no contact, 2-contact
        """
        self._step = step

        #   current contact velocity at time t
        self._dq_n, self._dq_t = self.contact_velocity(q)

        #   calculate distance between joint centers and delta (penetration)
        self._distance, self._delta, self._n_GCS, self._t_GCS = self._contact_geometry_GCS(q)

        #   if distance is greater than radial clearance, contact is present
        if self.check_for_contact_continued_condition(self._delta, self._dq_n, q):
            self.status = 1

        else:
            self.status = 0
            self.contact_detected = False
            self.contact_distance_inside_tolerance = False
            # self.list_of_contact_force_objects_constructed = False
            self._contact_point_found = False
            self.initial_contact_velocity_calculated = False

            self.no_contact()

        return self.status

    def _get_contact_geometry_data(self, q):
        """
        Function calculates normal and tangent vector of body i, j in GCS
        """
        for contact_point in self._contact_point_obj_list:
            contact_point.contact_geometry(q)

    def _contact_geometry_LCS(self, q):
        """
        Function evaluates contact geometry parameters in body LCS based on contact geometry in GCS
        Function evaluates:
        self._n_LCS_list:   normal of contact in body LCS:
        self._t_LCS_list:   tangent of contact in body LCS
        self.u_P_LCS_list:  contact point in body LCS
        :param q:
        :return:
        """
        for contact_point in self._contact_point_obj_list:
            contact_point._contact_geometry_LCS(q)

    def _contact_geometry_GCS(self, q):
        """

        :param q:
        :return:
        """
        _distance = self.distance_TOL
        _distance_list = []
        _delta_list = []
        _n_GCS_list = []
        _t_GCS_list = []

        #   penetration depth
        for contact_point in self._contact_point_obj_list:
            if contact_point._contact_point_found:
                distance, delta, n_GCS, t_GCS = contact_point.update_contact_points_GCS(q)

                _distance_list.append(distance)
                _delta_list.append(delta)
                _n_GCS_list.append(n_GCS)
                _t_GCS_list.append(t_GCS)

        if hasattr(self, "_distance_obj_list"):
            for distance_obj in self._distance_obj_list:
                _delta = distance_obj._distance_sign
                _delta_list.append(_delta)
                _distance_list.append(_delta)

                #   contact normal in GCS
                _n_GCS = distance_obj.normal
                _n_GCS_list.append(_n_GCS)

                #   contact tangent in GCS
                _t_GCS = distance_obj.tangent
                _t_GCS_list.append(_t_GCS)

        return _distance_list, _delta_list, _n_GCS_list, _t_GCS_list

    def testing(self):
        """

        :return:
        """
        print "testing@", __name__
        print "BREAK CONTACT!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print "contact_point_obj =", self._contact_point_obj
        print "_contact_point_obj_list =", self._contact_point_obj_list

        for contact_point in self._contact_point_obj_list:
            if contact_point.active:
                contact_point.deactivate_forces()

        # for contact_point in self._contact_point_obj_list:
        #     if not contact_point.active:
        #         self._contact_point_obj_list.remove(contact_point)


if __name__ == "__main__":
    contact = None
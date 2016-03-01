__author__ = 'lskrinjar'
import itertools
from operator import attrgetter
from pprint import pprint

import numpy as np

from MBD_system.Ai_ui_P import Ai_ui_P_vector
from MBD_system.cad2cm_lcs import cad2cm_lcs
from MBD_system.contact.contact import Contact
from MBD_system.contact.distance.distance_line_node import DistanceLineNode
from MBD_system.contact.distance.distance_revolute_clearance_joint import DistanceRevoluteClearanceJoint
from MBD_system.n2t import n2t
from MBD_system.q2R_i import q2R_i
from MBD_system.q2theta_i import q2theta_i
from MBD_system.t2n import t2n
from MBD_system.transform_cs import gcs2cm_lcs, uP_gcs2lcs
from MBD_system.u_P_gcs2lcs import u_P_gcs2lcs
from MBD_system.u_P_lcs2gcs import u_P_lcs2gcs
from simulation_control_widget.opengl_widget.marker.marker import Marker


class PinSlotClearanceJointLinear(Contact):
    """
    classdocs
    """
    __id = itertools.count(0)

    def __init__(self, _type, body_id_i, body_id_j, u_iP, u_jP, u_jR, R0_i, h0_jP, L=None, h0_jR=None, properties_dict=[], parent=None):
        """

        :param _type:
        :param body_id_i:   id of a body - pin
        :param body_id_j:   id of a body - slot
        :param u_iP:        point vector to center of a pin on body i (in body CAD CS)
        :param u_jP:        point vector to start point center of the slot on body j (in body CAD CS)
        :param u_jR:        point vector to end point center of the slot on body j (in body CAD CS)
        :param R0_i:        radius of a pin on body i
        :param h0_jP:       width of a slot at start point center on body i
        :param h0_jR:       width of a slot at end point center on body j
        :return:
        """
        #    number
        self.contact_id = self.__id.next()
        self._name = "Pin_Slot_Clearance_Joint_"+str(self.contact_id)

        super(PinSlotClearanceJointLinear, self).__init__(_type, body_id_i, body_id_j, name=self._name, properties_dict=properties_dict, parent=parent)

        #   parent
        self._parent = parent

        #   type of contact
        self._type = "pin-slot clearance joint linear"

        #   contact model

        #   add all to list of CP - center points
        #   1 - body i (pin) center
        #   2 - body j (slot) center uP
        #   3 - body j (slot) center uR
        self.u_CP_LCS_list = []

        #   body i
        #   body i (pin) geometry parameters
        if self._parent is not None:
            if hasattr(self._parent._parent.bodies[body_id_i], "CM_CAD_LCS"):
                self.u_iP = cad2cm_lcs(u_iP, self._parent._parent.bodies[body_id_i].CM_CAD_LCS, 0)
        else:
            self.u_iP = u_iP
        self.u_CP_LCS_list.append(self.u_iP)
        #   radius of a pin
        self.R0_i = R0_i

        #   list of radii of pin and of slot
        self.R0_list = [self.R0_i, h0_jP / 2.]

        #   body j
        #   body j (slot) geometry data
        #   predefine empty list of slot centers
        self.u_CP_slot_LCS = []

        if self._parent is not None:
            if hasattr(self._parent._parent.bodies[body_id_j], "CM_CAD_LCS"):
                self.u_jP = cad2cm_lcs(u_jP, self._parent._parent.bodies[body_id_j].CM_CAD_LCS, 0)
                self.u_jR = cad2cm_lcs(u_jR, self._parent._parent.bodies[body_id_j].CM_CAD_LCS, 0)
        else:
            #   vector to first point of center of a slot in
            self.u_jP = u_jP
            #   vector to second point of center of a slot
            self.u_jR = u_jR
        #   append to centers list of a slot
        self.u_CP_slot_LCS.append(self.u_jP)
        self.u_CP_slot_LCS.append(self.u_jR)
        #   append to centers list
        self.u_CP_LCS_list.append(self.u_jP)
        self.u_CP_LCS_list.append(self.u_jR)
        #   in GCS
        self.u_CP_GCS_list = [np.zeros(2), np.zeros(2), np.zeros(2)]

        #   body ids for three point (1 point pin, 2 points slot)
        self.body_id_extended_list = [self.body_id_i, self.body_id_j, self.body_id_j]
        #   width of a slot at first point of center of a slot
        self.h0_jP = h0_jP
        #   width of a slot at second point of center of a slot
        if h0_jR is None:
            self.h0_jR = h0_jP
        else:
            self.h0_jR = h0_jR

        #   clearance
        self._radial_clearance = (self.h0_jR / 2.) - self.R0_i

        #   tangent of slot edge
        self.t_LCS = self.u_jR - self.u_jP
        #   length of flat edge of slot
        self.L = np.linalg.norm(self.t_LCS, ord=2)
        #   unit vector of direction of edge
        self.t_LCS_unit_vector = self.t_LCS / self.L
        #   list of tangents in LCS or edge vectors
        self.t_edge_slot_LCS = []
        self.t_edge_slot_LCS.append(self.t_LCS_unit_vector)
        self.t_edge_slot_LCS.append(-self.t_LCS_unit_vector)

        #   normal unit vector
        self.n_LCS_unit_vector = t2n(self.t_LCS_unit_vector)

        #   list of normals in LCS of flat edges in slot- section jPjR
        self.n_edge_slot_LCS = []
        self.n_edge_slot_LCS.append(self.n_LCS_unit_vector)
        self.n_edge_slot_LCS.append(-self.n_LCS_unit_vector)

        #   create slot frame geometry for flat surface
        self.frame_nodes_LCS = self._create_slot_frame_flat_surface(self.u_jP, self.u_jR, self.n_LCS_unit_vector, self.t_LCS_unit_vector, self.h0_jP)

        #   list of markers
        self.markers = self._create_markers()

        #   reset additional parameters before simulation
        self.reset()

    def _create_markers(self):
        """
        Function creates markers
        :return:
        """
        markers = []
        #   create markers for joint centers
        for i, u_P in enumerate(self.u_CP_LCS_list):
            #   node coordinates
            _node = np.array(np.append(u_P, self.z_dim), dtype="float32")

            #   get correct body object
            if i == 0:
                body = self.body_list[0]
            if i == 1 or i == 2:
                body = self.body_list[1]

            #   create marker object
            _marker = Marker(_node, visible=True, parent=body)

            #   append marker to list of body markers
            body.markers.append(_marker)
            #   append marker to list of contact markers
            markers.append(_marker)

        return markers

    def reset(self):
        """
        Function is used to reset simulation parameters before new simulation run
        :return:
        """
        #   contact status
        #   left half cylinder
        self.pin_in_section_jP = False
        #   flat section
        self.pin_in_section_jPjR = False
        #   right half cylinder
        self.pin_in_section_jR = False

        #   predefined list of contact point on each body (one contact point on a body)
        self.u_P_list_LCS = np.array([np.zeros(2), np.zeros(2)])

    def check_for_contact(self, q):
        """
        Check for contact
        :param q:
        :return:
        """
        # print "check_for_contact()"
        self._distance, self._delta = self._contact_geometry_GCS(q)
        # print "self._delta =", self._delta

        # print "self._delta@check_for_contact() =", self._delta

        #   add distance value to container

        # self._distance_solution_container = np.append(self._distance_solution_container, self._delta)

        #   check sign
        self._sign_check = np.sign(self._delta * self._distance_solution_container[-1])

        #    contact has happened, but time step has to be reduced as initial penetration depth is too large
        # if (np.sign(self._dqn_solution_container[-1]) == +1) or (self._dqn_solution_container[-1] == 0) and (self._sign_check == -1) and (self._distance >= self._radial_clearance):
        # print "(self._sign_check == -1) and (self._delta <= 0) =", (self._sign_check == -1) and (self._delta <= 0)
        if (self._sign_check == -1) and (self._delta <= 0):
            # print "sign and delta"
            #    beginning of contact detected, all parameters are within limits
            if abs(self._delta) <= self.distance_TOL:
            # else:
                self._delta0 = self._delta
                self.contact_detected = True
                self.contact_distance_inside_tolerance = True
                self.status = 1
                # print "---------------------------------------------"
                # print "contact detected", "L =", self.L
                # print "self._delta =", self._delta
                # print "DA =", self._pin_distance_to_DA_obj._distance
                # print "BC =", self._pin_distance_to_BC_obj._distance
                # print "section jPjR =", self._pin_distance_to_DA_obj._distance < self.L, self._pin_distance_to_BC_obj._distance < self.L
                # print "section jP =", self._pin_distance_to_DA_obj._distance < self.L, self._pin_distance_to_BC_obj._distance >= self.L
                # print "section jR =", self._pin_distance_to_DA_obj._distance >= self.L, self._pin_distance_to_BC_obj._distance < self.L
                # print "t =", self.t, "step =", self._step, "pin in: ",
                # if self.pin_in_section_jPjR:
                #     print "flat section"
                # if self.pin_in_section_jP:
                #     print "section jP"
                # if self.pin_in_section_jR:
                #     print "section jR"


                # print "contact point =",
                # print int(self._step_num_solution_container[-1]), self.status, t, self._distance,
                return 1

            #   step back
            if abs(self._delta) > self.distance_TOL:
                self.contact_detected = True
                self.status = -1
                # self._distance_solution_container = np.delete(self._distance_solution_container, -1)
                return -1

        #    all calculated distances are greater than tolerance and bodies are not in contact
        self.status = 0
        # self._status_container = np.append(self._status_container, self.status)
        self.no_contact()

        return 0

    def contact_update(self, step, t, q):
        """

        :param step:
        :param t:
        :param q:
        :return:
        """
        # print "contact_update ="
        #   update time step
        self._step = step

        #   current contact velocity at time t
        self._dq_n, self._dq_t = self._contact_velocity(q)

        #   calculate distance between joint centers and delta (penetration)
        self._distance, self._delta = self._contact_geometry_GCS(q)

        # print "self._delta@ =", self._delta,
        # print "self._delta0 =", self._delta0
        # print self._delta <= self._delta0
        # print "self._contact_point_found =", self._contact_point_found
        # time.sleep(1)
        if self._delta <= self._delta0:#(self._distance >= self._radial_clearance) and (abs(self._delta) >= self.distance_TOL):
            # print "exe??"
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
        Function locates and creates contact point vector in GCS based on position of both bodies in contact
        :param q:
        :return:
        """
        # print "_get_contact_geometry_data()"
        if self.pin_in_section_jPjR:
            #   contact point on body j - slot
            #   GCS
            self.u_P_GCS_list[1] = self.u_jP_GCS = self._distance_obj.contact_point_on_line_GCS()
            #   LCS
            self.u_jP_LCS = u_P_gcs2lcs(self.u_jP_GCS, q, self.body_id_j)
            # print "self.u_jP_LCS =", self.u_jP_LCS

            #   get normal in GCS
            self._n_GCS = self._distance_obj.normal

            # print "self._n_GCS =", self._n_GCS

        elif self.pin_in_section_jP or self.pin_in_section_jR:
            #   get normal in GCS
            self._n_GCS = self._distance_obj._distance_vector / self._distance_obj._distance

        else:
            raise ValueError, "Contact in slot section not detected!"

        #   get tangent in CS
        self._t_GCS = n2t(self._n_GCS)
        self._contact_point_found = True

    def _contact_geometry_GCS(self, q):
        """
        Function evaluates contact geometry in GCS
        Function evaluates:
        self._n_GCS:    normal of contact in GCS:
        self._t_GCS:    tangent of contact in GCS
        self.node_GCS:  contact point in GCS (on body i and on body j - unverified)
        :param q:
        :return delta:
        :return distance:
        """
        # print "_contact_geometry_GCS()"
        #   evaluate coordinates in GCS of mid frame section and pin center
        pin_CP_GCS, slot_CP_GCS, slot_frame_nodes_GCS, slot_frame_normals_GCS, slot_frame_tangents_GCS = self._contact_geometry_CP_GCS(q)

        #   center points in GCS
        self.u_CP_GCS_list = [pin_CP_GCS, slot_CP_GCS[0], slot_CP_GCS[1]]
        # print "self.u_CP_GCS_list =", self.u_CP_GCS_list
        # time.sleep(100)
        #   evaluate position of pin in slot:
        #   section jP
        #   section jPjR
        #   section jR
        self.pin_in_section_jP, self.pin_in_section_jPjR, self.pin_in_section_jR = self._check_pin_position(slot_frame_nodes_GCS, pin_CP_GCS)


        #   evaluate contact distance and penetration depth based on in what part of the slot is the pin
        if self.pin_in_section_jPjR:
            #   define empty list to store 2 distance object
            #   distance of pin CP to edge AB
            #   distance of pin CP to edge CD
            distance_obj_list = []
            # print "--------------------"
            for i, (_n, _t) in enumerate(zip(slot_frame_normals_GCS, slot_frame_tangents_GCS)):
                _distance_obj = DistanceLineNode(u_iP=slot_frame_nodes_GCS[2*i, :], u_jP=pin_CP_GCS, normal=_n, u_iR=slot_frame_nodes_GCS[2*i+1, :], tangent=_t)
                distance_obj_list.append(_distance_obj)

                # print "i =", i, "_distance =", _distance_obj._distance
            #   get distance object with min distance attribute value
            self._distance_obj = min(distance_obj_list, key=attrgetter('_distance'))

            #   add edge id to object to have info with which edge of flat slot surface is contact with a pin
            self._distance_obj.edge_id = distance_obj_list.index(self._distance_obj)

            #   get distance attribute from distance object
            _distance = self._distance_obj._distance
            #   evaluated penetration depth
            _delta = self._evaluate_delta_flat_section(_distance)
            # print "_delta =", _delta
            #   tangent in GCS
            # pprint(vars(self._distance_obj))
            # time.sleep(100)

        else:
            if self.pin_in_section_jP:
                #   calculate distance between axis of both bodies in half revolute clearance joint at point Pj
                _slot_CP_GCS = slot_CP_GCS[0]

            elif self.pin_in_section_jR:
                #   calculate distance between axis of both bodies in half revolute clearance joint at point Rj
                _slot_CP_GCS = slot_CP_GCS[1]

            else:
                raise ValueError, "Contact point not in section!"

            # print "pin_CP_GCS =", pin_CP_GCS
            # print "_slot_CP_GCS =", _slot_CP_GCS
            self._distance_obj = DistanceRevoluteClearanceJoint(_slot_CP_GCS, pin_CP_GCS, parent=self)

            #   penetration depth is difference between nominal radial clearance and actual calculated clearance at time t
            _distance = self._distance_obj._distance
            # print "_distance =", _distance
            _delta = self._evaluate_delta_cylindrical_section(_distance)
            # print "_delta =", _delta
        #   if penetration depth is below small value, contact is present and contact geometry parameters are evaluated
        # if _delta <= self._delta0 and not self._contact_point_found:
        #     if self.pin_in_section_jPjR:
        #         #   contact normal in GCS
        #         self._n_GCS = self._distance_obj.normal

        return _distance, _delta

    def _contact_geometry_LCS(self, q):
        """
        Function evaluates contact geometry parameters in body LCS based on contact geometry in GCS
        Function evaluates:
        self._n_LCS_list:   normal of contact in body LCS:
        self._t_LCS_list:   tangent of contact in body LCS
        self.u_P_LCS_list:  contact point in body LCS
        """
        # print "_contact_geometry_LCS()"
        if self.pin_in_section_jPjR:
            #   contact point on body i - pin
            #   in LCS
            self.u_P_LCS_list[0] = self.u_iP_LCS = Ai_ui_P_vector(-self._n_GCS * self.R0_i, q2theta_i(q, self.body_id_i))
            #   in GCS
            self.u_P_GCS_list[0] = self.u_iP_GCS = q2R_i(q, self.body_id_i) + self.u_iP_LCS

            #   contact point on body j in GCS
            self.u_P_GCS = self.u_P_GCS_list[1] = self.u_jP_GCS = self._distance_obj.contact_point_on_line_GCS()
            # time.sleep(10)
            #   contact point on body j in LCS
            self.u_P_LCS_list[1] = self.u_jP_LCS = gcs2cm_lcs(self.u_P_GCS, q2R_i(q, self.body_id_j), q2theta_i(q, self.body_id_j))

            #   in GCS
            # self.u_iP_GCS = u_P_lcs2gcs(self.u_iP_LCS, q, self.body_id_i) + (-self._n_GCS) * self.R0_i
            # print "self.u_iP_GCS =", self.u_iP_GCS

            #   get normal and tangent of each body
            for i, (sign, body_id) in enumerate(zip([-1, +1], self.body_id_list)):
                theta = q2theta_i(q, body_id)
                #   evaluate body normal in body LCS
                self._n_LCS_list[i] = gcs2cm_lcs(sign * self._n_GCS, theta=theta)
                #   evaluate body tangent in body LCS
                self._t_LCS_list[i] = gcs2cm_lcs(sign * self._t_GCS, theta=theta)

        elif self.pin_in_section_jP or self.pin_in_section_jR:
            #   create normal list in LCS
            self._n_GCS_list = [self._n_GCS, -self._n_GCS]
            for i, (body_id, _normal) in enumerate(zip(self.body_id_list, self._n_GCS_list)):
                #   normal in LCS
                _theta = q2theta_i(q, body_id)
                # _normal_LCS = Ai_ui_P_vector(_normal, _theta)
                _normal_LCS = uP_gcs2lcs(_normal, _theta)
                #   append normal to list
                self._n_LCS_list[i] = _normal_LCS

            if self.pin_in_section_jP:
                _u_CP_LCS_list = [self.u_CP_LCS_list[0], self.u_CP_LCS_list[1]]
                _u_CP_GCS_list = [self.u_CP_GCS_list[0], self.u_CP_GCS_list[1]]

            if self.pin_in_section_jR:
                _u_CP_LCS_list = [self.u_CP_LCS_list[0], self.u_CP_LCS_list[2]]
                _u_CP_GCS_list = [self.u_CP_GCS_list[0], self.u_CP_GCS_list[2]]

            #   evaluate actual contact point in LCS of each body and in GCS
            for i, (body_id, _u_CP_LCS, _R0) in enumerate(zip(self.body_id_list, _u_CP_LCS_list, self.R0_list)):
                #   R of body
                R = q2R_i(q, body_id)
                #   theta of body
                theta = q2theta_i(q, body_id)

                #   contact point in body LCS
                _u_P_LCS = _u_CP_LCS
                self.u_P_LCS_list[i] = _u_P_LCS# + _R0 * self._n_GCS
                #   contact point in GCS
                self.u_P_GCS_list[i] = R + Ai_ui_P_vector(_u_CP_LCS, theta) + _R0 * self._n_GCS

            # print "self.u_P_GCS_list =", self.u_P_GCS_list


    def _evaluate_delta_flat_section(self, distance):
        """
        Evaluate penetration depth - delta for contact between pin and flat surface of the slot
        :param distamce:
        :return:
        """
        delta = distance - self.R0_i#((self.h0_jP - 2 * ) / 2.) -
        return delta

    def _evaluate_delta_cylindrical_section(self, distance):
        """
        Evaluate penetration depth - delta for contact between pin and cylindrical surface of the slot
        :param distamce:
        :return:
        """
        # print "_evaluate_delta_cylindrical_section()"
        # print "distance =", distance
        # print "self._radial_clearance =", self._radial_clearance
        delta = self._radial_clearance - distance
        return delta

    def _check_pin_position(self, frame_nodes_GCS, pin_CP_GCS):
        """
        Function eveluates where in slot is pin located to use correct equation for evaluation of penetration depth
        :param frame_nodes_GCS:     position of slot frame nodes in GCS
        :param pin_CP_GCS:          position of pin center in GCS
        :return:
        """
        # print "_check_pin_position()"
        #   set default values
        pin_in_section_jP = pin_in_section_jPjR = pin_in_section_jR = False

        #   distance of pin center to line DA and BC
        #   node A
        u_Aj = frame_nodes_GCS[0, :]
        #   node B
        u_Bj = frame_nodes_GCS[1, :]
        #   node C
        u_Cj = frame_nodes_GCS[2, :]
        #   node D
        u_Dj = frame_nodes_GCS[3, :]
        #   DA
        self._pin_distance_to_DA_obj = DistanceLineNode(u_Aj, pin_CP_GCS, u_iR=u_Dj, parent=self)
        # print "_pin_distance_to_DA_obj =", _pin_distance_to_DA_obj._distance
        #   BC
        self._pin_distance_to_BC_obj = DistanceLineNode(u_Cj, pin_CP_GCS, u_iR=u_Bj, parent=self)
        # print "_pin_distance_to_BC_obj =", _pin_distance_to_BC_obj._distance

        if self._pin_distance_to_DA_obj._distance < self.L and self._pin_distance_to_BC_obj._distance < self.L:
            pin_in_section_jPjR = True
            # print "pin in jPjR"

        if self._pin_distance_to_DA_obj._distance < self.L and self._pin_distance_to_BC_obj._distance >= self.L:
            pin_in_section_jP = True
            # print "pin in jP"

        if self._pin_distance_to_DA_obj._distance >= self.L and self._pin_distance_to_BC_obj._distance < self.L:
            pin_in_section_jR = True
            # print "pin in jR"
        #   section from (u_)jP to (u_)jR
        # if (frame_nodes_GCS[0, 0] < pin_CP_GCS[0] < frame_nodes_GCS[1, 0]) and (frame_nodes_GCS[0, 1] < pin_CP_GCS[1] < frame_nodes_GCS[2, 1]):
        #     pin_in_section_jPjR = True
        #     # print "pin in jPjR"
        #
        # #   half-circlurar section at (u_)jP
        # if (pin_CP_GCS[0] <= max(frame_nodes_GCS[0, 0], frame_nodes_GCS[3, 0])) and (pin_CP_GCS[1] <= max(frame_nodes_GCS[0, 1], frame_nodes_GCS[3, 1])):
        #     # print "pin_CP_GCS[0] =", pin_CP_GCS[0]
        #     # print "max x =", max(frame_nodes_GCS[0, 0], frame_nodes_GCS[3, 0])
        #     # print "pin_CP_GCS[1] =", pin_CP_GCS[1]
        #     # print "max y =", max(frame_nodes_GCS[0, 1], frame_nodes_GCS[3, 1])
        #     pin_in_section_jP = True
        #     # print "pin in jP"
        #
        # #   half-circlurar section at (u_)jR
        # if (pin_CP_GCS[0] >= min(frame_nodes_GCS[1, 0], frame_nodes_GCS[2, 0])) and (pin_CP_GCS[1] >= min(frame_nodes_GCS[1, 1], frame_nodes_GCS[2, 1])):
        #     pin_in_section_jR = True
        #     # print "pin in jR"

        return pin_in_section_jP, pin_in_section_jPjR, pin_in_section_jR

    def _contact_geometry_CP_GCS(self, q):
        """

        :param q:
        :return:
        """
        # print "_contact_geometry_CP_GCS()"
        #   calculate position of pin/hole centers of each body in GCS
        #   this is a list of CP - center points
        #   1 -
        u_CP_list_GCS = []

        #   center points in GCS
        # for uP, id in zip(self.u_CP_LCS_list, self.body_id_extended_list):
        #     uP_gcs = q2R_i(q, id) + Ai_ui_P_vector(uP, q2theta_i(q, id))
        #     print "uP_gcs =", uP_gcs

        # print "_update_slot_frame_flat_surface"
        #   slot frame nodes in GCS
        frame_nodes_GCS, frame_normals_GCS, frame_tangents_GCS = self._slot_frame_flat_surface_GCS(q)
        # print "frame_nodes_GCS =", frame_nodes_GCS

        #   pin center in GCS
        pin_CP_GCS = self._pin_GCS(q)
        # print "pin_CP_GCS =", pin_CP_GCS

        #   slot center points in GCS
        slot_CP_GCS = self._slot_GCS(q)

        return pin_CP_GCS, slot_CP_GCS, frame_nodes_GCS, frame_normals_GCS, frame_tangents_GCS

    def _pin_GCS(self, q):
        """
        Function updates and calculates position of a pin center in global coordinates - GCS
        :param q:
        :return:
        """
        # CP_GCS = q2R_i(q, self.body_id_i) + Ai_ui_P_vector(self.u_iP, q2theta_i(q, self.body_id_i))
        CP_GCS = u_P_lcs2gcs(self.u_iP, q, self.body_id_i)
        return CP_GCS

    def _slot_GCS(self, q):
        """

        :param q:
        :return:
        """
        slot_CP_GCS = []
        for CP in self.u_CP_slot_LCS:
            _uP_gcs = u_P_lcs2gcs(CP, q, self.body_id_j)
            slot_CP_GCS.append(_uP_gcs)
        return slot_CP_GCS

    def _create_slot_frame_flat_surface(self, uPj, uRj, n, t, hPj, hRj=None):
        """
        Direction of tangent uPj->hRj
        Direction of normal uPj->T1
        Frame shape:
        A------------------>B
        |                   |
        uPj                 uRj
        |                   |
        D<------------------C
        :param q:                   a vector (numpy array) of positions
        :return frame_nodes_LCS:    a matrix (numpy array) of frame nodes
        """
        #   4 frame nodes in a matrix
        frame_nodes_LCS = np.array([uPj - n*hPj*.5,     #0   A
                                    uRj - n*hPj*.5,     #1   B
                                    uRj + n*hPj*.5,     #2   C
                                    uPj + n*hPj*.5],    #3   D
                                    dtype='float32')

        return frame_nodes_LCS

    def _slot_frame_flat_surface_GCS(self, q):
        """
        Function updates and calculates position of slot frame in global coordinates - GCS
        :param q:
        :return:
        """
        #   nodes in GCS
        frame_nodes_GCS = q2R_i(q, self.body_id_j) + Ai_ui_P_vector(self.frame_nodes_LCS, q2theta_i(q, self.body_id_j))

        #   normals in GCS
        frame_normals_GCS = Ai_ui_P_vector(np.array(self.n_edge_slot_LCS), q2theta_i(q, self.body_id_j))

        #   tangents in GCS
        frame_tangents_GCS = Ai_ui_P_vector(np.array(self.t_edge_slot_LCS), q2theta_i(q, self.body_id_j))

        return frame_nodes_GCS, frame_normals_GCS, frame_tangents_GCS

if __name__ == "__main__":
    _type = "pin-slot clearance joint linear"
    body_id_i = 0
    body_id_j = 1
    u_iP = np.array([0, 0])*1E-3
    u_jP = np.array([10, 20])*1E-3
    u_jR = np.array([70, 20])*1E-3
    R0_i = 4E-3
    h0_jP = 10E-3
    _joint = PinSlotClearanceJointLinear(_type, body_id_i, body_id_j, u_iP, u_jP, u_jR, R0_i, h0_jP)
    pprint(vars(_joint))
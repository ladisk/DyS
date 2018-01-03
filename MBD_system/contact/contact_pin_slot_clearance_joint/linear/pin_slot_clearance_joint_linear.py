__author__ = 'lskrinjar'
import os
import time
import itertools
import inspect
from operator import attrgetter
from pprint import pprint
import numpy as np
from matplotlib import pyplot as plt


from MBD_system.Ai_ui_P import Ai_ui_P_vector
from MBD_system.cad2cm_lcs import cad2cm_lcs
from MBD_system.contact.contact import Contact
from MBD_system.n2t import n2t
from MBD_system.q2R_i import q2R_i
from MBD_system.q2theta_i import q2theta_i
from MBD_system.t2n import t2n
from MBD_system.transform_cs import gcs2cm_lcs, uP_gcs2lcs
from MBD_system.u_P_gcs2lcs import u_P_gcs2lcs
from MBD_system.u_P_lcs2gcs import u_P_lcs2gcs
from MBD_system.q2dR_i import q2dR_i
from MBD_system.q2dtheta_i import q2dtheta_i
from MBD_system.dr_contact_point_uP import dr_contact_point_uP
from simulation_control_widget.vtk_widget.marker.marker import Marker
from MBD_system.contact.distance.distance_line_node import DistanceLineNode
from MBD_system.contact.distance.distance_RCJ import DistanceRCJ
from MBD_system.contact.distance.distance_PSCJ import DistancePSCJ
from MBD_system.contact.distance.contact_point import ContactPoint
from MBD_system.contact.distance.contact_point_PSCJ import ContactPointPSCJ
from MBD_system.contact.distance.contact_point_RCJ import ContactPointRCJ


class PinSlotClearanceJointLinear(Contact):
    """
    classdocs
    """

    def __init__(self, _type, body_id_i, body_id_j, u_iP, u_iR, u_jP, h0_iP, R0_j, L=None, h0_iR=None, properties_dict={}, parent=None):
        """
        contact_models[1] - on cylindrical section
        contact_models[0] - on flat section
        :param _type:       type of a clearance joint or contact
        :param body_id_i:   id of a body - slot
        :param body_id_j:   id of a body - pin
        :param u_iP:        point vector to center of a pin on body i (in body CAD CS)
        :param u_jP:        point vector to start point center of the slot on body j (in body CAD CS)
        :param u_jR:        point vector to end point center of the slot on body j (in body CAD CS)
        :param R0_i:        radius of a pin on body i
        :param h0_jP:       width of a slot at start point center on body i
        :param h0_jR:       width of a slot at end point center on body j
        :return:
        """
        #    name
        self._name = "Pin_Slot_Clearance_Joint_"

        super(PinSlotClearanceJointLinear, self).__init__(_type, body_id_i, body_id_j, name=self._name, properties_dict=properties_dict, parent=parent)

        #   parent
        self._parent = parent

        #   type of contact
        self._contact_type = "pin-slot clearance joint linear"

        #   contact model
        #   because of different contact geometries, different contact models can be used
        #   1 contact model for circular section
        #   2 contact model for flat section
        #   contact models can be the same with different contact properties, or completely different contact modeles for 1 and 2
        #   can be used (they are constructed in contact.py)

        #   add all to list of CP - center points (index-position per element)
        #   0 - body i (slot) center uP
        #   1 - body i (slot) center uR
        #   2 - body j (pin) center
        self.u_CP_LCS_list = []

        #   body j
        #   body j (slot) geometry data
        #   predefine empty list of slot centers
        self.u_CP_slot_LCS = []

        if self._parent is not None:
            if hasattr(self._parent._parent.bodies[body_id_i], "u_CAD"):
                self.u_iP = cad2cm_lcs(u_iP, self._parent._parent.bodies[body_id_i].u_CAD, 0)
                self.u_iR = cad2cm_lcs(u_iR, self._parent._parent.bodies[body_id_i].u_CAD, 0)
        else:
            #   vector to first point of center of a slot in
            self.u_iP = u_iP
            #   vector to second point of center of a slot
            self.u_iR = u_iR

        #   append to centers list of a slot
        self.u_CP_slot_LCS.append(self.u_iP)
        self.u_CP_slot_LCS.append(self.u_iR)

        #   append to centers list
        self.u_CP_LCS_list.append(self.u_iP)
        self.u_CP_LCS_list.append(self.u_iR)

        #   in GCS
        self.r_CP_GCS_list = [np.zeros(2), np.zeros(2), np.zeros(2)]

        #   vectors in GCS
        #   slot
        self.r_iP = None
        self.r_iR = None
        #   pin
        self.r_jP = None

        #   body ids for three point (1 point pin, 2 points slot)
        self.body_id_extended_list = [self.body_id_i, self.body_id_i, self.body_id_j]
        #   width of a slot at first point of center of a slot
        self.h0_iP = h0_iP
        self.R0_i = h0_iP / 2.
        #   width of a slot at second point of center of a slot
        if h0_iR is None:
            self.h0_iR = h0_iP
        else:
            self.h0_iR = h0_iR

        #   radius of a pin
        self.R0_j = R0_j

        #   list of radii of pin and of slot
        self.R0_list = [self.R0_i, self.R0_j]

        #   clearance
        self._radial_clearance = (self.h0_iP / 2.) - self.R0_j

        #   body j
        #   body j (pin) geometry parameters
        self.u_jP = u_jP
        if self._parent is not None:
            if hasattr(self._parent._parent.bodies[body_id_j], "u_LCS"):
                self.u_jP = cad2cm_lcs(u_jP, self._parent._parent.bodies[body_id_j].u_LCS, 0)
        self.u_CP_LCS_list.append(self.u_jP)

        #   tangent of slot edge
        self.t_LCS = self.u_iR - self.u_iP
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
        self.n_edge_slot_LCS = [+self.n_LCS_unit_vector,
                                -self.n_LCS_unit_vector]

        #   create slot frame geometry for flat surface
        #   LCS
        self.frame_nodes_LCS = self._create_slot_frame_flat_surface(self.u_iP, self.u_iR, self.n_LCS_unit_vector, self.t_LCS_unit_vector, self.h0_iP)
        #   GCS
        if self._parent is not None:
            self.slot_frame_nodes_GCS, self.frame_normals_GCS, self.frame_tangents_GCS = self._slot_frame_flat_surface_GCS(self._parent._parent.evaluate_q())

        #   create contact model
        self._create_contact_model(self.properties)

        #   create friction model
        self._create_friction_model(self.properties)

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
            #   get correct body object
            body_id = body = None
            if i == 0 or i == 1:
                body = self.body_list[0]
                body_id = self.body_id_list[0]
            if i == 2:
                body = self.body_list[1]
                body_id = self.body_id_list[1]

            #   node coordinates in GCS
            rP = u_P_lcs2gcs(u_P, self.q0, body_id)

            #   add 3rd dimension for visualization
            rP = np.array(np.append(rP, self.z_dim), dtype="float32")

            if body is not None:
                #   create marker object
                marker = Marker(rP, parent=self)

                #   append marker to list of contact markers
                markers.append(marker)

        return markers

    def _reset(self):
        """
        Function is used to reset simulation parameters before new simulation run
        :return:
        """
        #   simulation parameters
        self._step = 0
        self._t = 0

        #   contact status
        #   left half cylinder
        self.pin_in_section_iP = False
        #   flat section
        self.pin_in_section_iPiR = False
        #   right half cylinder
        self.pin_in_section_iR = False

        self.pin_in_section = ""

        self._contact_point_obj = None
        self._distance_obj = None

        #   predefined list of contact point on each body (one contact point on a body)
        self.u_P_list_LCS = np.array([np.zeros(2), np.zeros(2)])

    def _solve_ECF_N(self, t, q, _delta, _dq_n, _dq_t):
        """
        Function evaluates contact forces based on selected contact and friction model
        :param t:       time (np.float)
        :param q:       vector of displacements and velocities (np.array)
        :param _delta:  penetration depth (np.float)
        :param _dq_n:   relative normal contact velocity (np.float)
        :param _dq_t:   relative tangent contact velocity (np.float)
        :returns None:
        """
        # print "_solve_ECF_N()"
        #   check if contact is finished
        #   contact
        # if self.check_for_contact_continued_condition(_delta, _dq_n, q):
        # for i, (contact_point, delta) in enumerate(zip(self._contact_point_obj_list, _delta)):

            # if delta <= self._delta0[i]:
        #   normal contact force
        Fn = self.contact_model.contact_force(_delta, self._contact_point_obj._dq_n, dq0_n=self._contact_point_obj._dq0_n)
        #   tangent contact force
        Ft = self.friction_model.friction_force(self._contact_point_obj.Fn, self._contact_point_obj._dq_t)
            # else:
            #     Fn = 0.
            #     Ft = 0.
        if self._step > 2164:
            print "delta =", _delta
            print "Fn =", Fn
        # print "delta =", delta, "Fn =", Fn#, "Ft =", Ft, "_dq_n =", _dq_n
        self._contact_point_obj.evaluate_F(Fn, Ft)

        #   no contact
        # else:
        #     print "CONTACT FINISHED - self._step =", self._step
        #
        #     self._contact_point_found = False
        #     self.initial_contact_velocity_calculated = False
        #     self.contact_distance_inside_tolerance = False
        #     self.contact_detected = False
        #     self.status = 0
        #     self._delta_n = self._delta
        #     print "self._delta_n =", self._delta_n
        #
        #     #   delete force objects from contact force lists (normal and tangent direction)
        #     self.no_contact()

        #   update contact forces
        self._update_contact_forces(q)

        # #   create, update or delete force object of contact forces at contact points
        # for i, contact_point in enumerate(self._contact_point_obj_list):
        #     if contact_point._contact_point_found and contact_point.active:
        #         for body_id, Fn_i, Ft_i, u_P_i, n_i, t_i in zip(contact_point.body_id_list, contact_point._Fn_list, contact_point._Ft_list, contact_point.u_P_LCS_list, contact_point._n_GCS_list, contact_point._t_GCS_list):
        #
        #             #   print only contact point of body j - pin
        #             #   normal force
        #             Fn_i.update_F(q=q, step=self._step, F=contact_point.Fn * n_i, u_P=u_P_i)
        #             Fn_i._visible = True
        #             #   tangent force
        #             Ft_i.update_F(q=q, step=self._step, F=contact_point.Ft * t_i, u_P=u_P_i)
        #             Ft_i._visible = True

    def _no_contact(self):
        """

        :return:
        """
        # print "_no_contact()"
        # print "self._delta =", self._delta
        for contact_point in self._contact_point_obj_list:
            if contact_point.active:
                contact_point.deactivate_forces()
            else:
                contact_point.remove_forces()

                self._contact_point_obj_list.remove(contact_point)

        # self._distance_obj_list = []
        # self._contact_point_obj_list = []
        self._contact_point_obj = None

        # if self._delta[0] < 0.:
        #     self._delta = self._delta
        #
        # else:
        #     self._delta[0] = self.distance_TOL

    def check_for_contact_started_condition(self, delta, sign):
        """
        Function checks condition for contact
        """
        # print "check_for_contact_started_condition()"
        # contact_started = False
        # print "delta =", delta
        # print "self._distance_obj._contact =", self._distance_obj._contact
        # print "sign =", sign
        # if :
        #     print "contact_started!", delta, sign
        if (sign == -1) and (delta <= 0):
            # print "1"
            contact_started = True
        elif (sign == +1) and (delta <= 0) and (np.sign(self._distance_solution_container[-1]) == -1) and (self._distance_solution_container[-1] > delta):
            # print "2"
            contact_started = True
        else:
            # print "3"
            contact_started = False

        return contact_started

    def check_for_contact_continued_condition(self, delta, dq_n, q):
        """
        Function checks condition for contact to continue
        :param delta:
        :return:
        """
        # print "check_for_contact_continued_condition()"
        # print "PIN IN SECTION @check_for_contact_continued_condition() =", self._contact_point_obj_list[0].pin_in_section
        # print "self._delta0 =", self._delta0
        # print "self._dq0_n =", self._dq0_n
        # print "dq_n =", dq_n

        if delta <= self._delta0:#0., self._delta0
            continued = True
        # elif abs(dq_n) <= self._dq0_n and self._dq0_n > 0. and (delta[0] < self._delta0):  # and (dq_n < self._dq0_n)
        elif abs(dq_n) <= self._dq0_n and (delta < self._delta0):# and (dq_n < self._dq0_n)
            continued = True
        # elif abs(dq_n) <= self._dq0_n and delta[0] < 0.:#self._delta0)# and (dq_n < self._dq0_n)
        #     continued = True

        else:
            continued = False
            #   delta at break of contact
            self._delta_n = self._delta

            print "CONTACT FINISHED - self._step =", self._step
            print "self._delta0 =", self._delta0
            print "self._delta_n =", self._delta_n


        # print "continued =", continued, "self._delta0 =", self._delta0, "delta[0] =", delta[0]
        return continued

    def check_for_contact(self, step, t, q):
        """
        Check for contact
        :param q:
        :return:
        """
        #   step
        self._step = step

        #   time
        self._t = t

        #   evaluate distance, delta
        self._distance, self._delta, self._n_GCS, self._t_GCS = self._contact_geometry_GCS(q)
        # print "self._distance, self._delta, self._n_GCS, self._t_GCS =", self._distance, self._delta, self._n_GCS, self._t_GCS
        #   check sign
        # print "self._delta =", self._delta
        # print "self._distance_solution_container[-1] =", self._distance_solution_container[-1]
        self._sign_check = np.sign(self._delta * self._distance_solution_container[-1])

        #    contact has happened, but time step has to be reduced as initial penetration depth is too large
        if self.check_for_contact_started_condition(self._delta, self._sign_check):
            # print "self._sign_check =", self._sign_check
            #    beginning of contact detected, all parameters are within limits
            if (abs(np.array(self._delta)) <= self.distance_TOL).all() or (self._delta < self._distance_solution_container[-1]):
                print "CONTACT DETECTED - step =", self._step
                print "PIN IN SECTION =", self.pin_in_section
                # print "normal ="
                self._delta0 = self._delta
                print  "self._delta0 =", self._delta0
                # pprint(vars(self._distance_obj))
                self.contact_detected = True
                self.contact_distance_inside_tolerance = True
                self.status = 1
                self._distance_obj_list = [self._distance_obj]
                self._contact_point_obj = ContactPointPSCJ(self._distance_obj, scale=self.scale, parent=self)
                # self._distance_obj = None
                # print "_n_GCS_list =", self._contact_point_obj._n_GCS_list
                # pprint(vars(self._contact_point_obj))_n_GCS_list
                # time.sleep(100)
                # #   create contact point object
                # if isinstance(self._distance_obj, DistanceRCJ):
                #     self._contact_point_obj = ContactPointRCJ(self._distance_obj, scale=self.scale, parent=self)
                #
                # if isinstance(self._distance_obj, DistancePSCJ):
                #     self._contact_point_obj = ContactPointPSCJ(self._distance_obj, scale=self.scale, parent=self)
                    # pprint(vars(self._contact_point_obj))

                self._contact_point_obj._contact_point_found = True
                # time.sleep(100)
                self._contact_point_obj_list = [self._contact_point_obj]
                return self.status

            #   step back
            if (abs(np.array(self._delta)) > self.distance_TOL).all():
            # if abs(self._delta) > self.distance_TOL:
                self.contact_detected = True
                self.status = -1
                return self.status

        #    all calculated distances are greater than tolerance and bodies are not in contact
        self.status = 0
        # self._status_container = np.append(self._status_container, self.status)
        self.no_contact()

        return self.status

    def contact_update(self, step, t, q):
        """

        :param step:
        :param t:
        :param q:
        :return:
        """
        #   update time step
        self._step = step

        #   calculate distance between joint centers and delta (penetration)
        self._distance, self._delta, self._n_GCS, self._t_GCS = self._contact_geometry_GCS(q)

        #   current contact velocity at time t
        self._dq_n, self._dq_t = self.contact_velocity(q)

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
        Function locates and creates contact point vector in GCS based on position of both bodies in contact
        :param q:
        :return:
        """
        #   evaluate coordinates in GCS of mid frame section and pin center
        pin_CP_GCS, slot_CP_GCS, self.slot_frame_nodes_GCS, slot_frame_normals_GCS, slot_frame_tangents_GCS = self._contact_geometry_CP_GCS(q)

        #   if contact point is already found only evaluation of its position is done
        if self._contact_point_found:
            self._distance, self._delta, self._n_GCS, self._t_GCS = self._contact_geometry_GCS(q)
            self._get_contact_geometry_data(q)

        self._contact_geometry_LCS(q)

    def _contact_geometry_GCS(self, q):
        """
        Function evaluates contact geometry in GCS
        Function evaluates:
        self._n_GCS:    normal of contact in GCS:
        self._t_GCS:    tangent of contact in GCS
        self.node_GCS:  contact point in GCS (on body i and on body j - unverified)
        :param q:               vector of states of MBD system
        :return _distance:      distance between center of pin and hole
        :paramm _delta:         depth of deformation (indentation, penetration)
        :return _n_GCS:         normal of contact in GCS
        _return _t_GCS:         tangent of contact in GCS
        """
        # print "_contact_geometry_GCS()@",__name__
        # num = self._parent._parent.solver.DAE_fun.number_of_evaluations

        #   evaluate coordinates in GCS of mid frame section and pin center
        pin_CP_GCS, slot_CP_GCS, self.slot_frame_nodes_GCS, slot_frame_normals_GCS, slot_frame_tangents_GCS = self._contact_geometry_CP_GCS(q)

        #   vectors in GCS
        self.r_iP = slot_CP_GCS[0]
        self.r_iR = slot_CP_GCS[1]
        self.r_jP = pin_CP_GCS

        #   center points in GCS
        # self.r_CP_GCS_list = [slot_CP_GCS[0], slot_CP_GCS[1], pin_CP_GCS]

        #   evaluate position of pin in slot:
        #   section jP
        #   section jPjR
        #   section jR
        #   distance object of pin slot clearance joint has 4 distance values as array
        if self._contact_point_obj is not None:
            self._contact_point_obj.update_contact_point(self.r_iP, self.r_iR, self.r_jP)
            pin_in_section = self._contact_point_obj.pin_in_section
            self._contact_point_obj.contact_points_GCS(frame_nodes_GCS=self.slot_frame_nodes_GCS)

        else:
            self._distance_obj = DistancePSCJ(self.r_iP, self.r_iR, self.r_jP, parent=self)
            pin_in_section = self._distance_obj.pin_in_section


        if pin_in_section != self.pin_in_section:
            self.section_changed = True
        else:
            self.section_changed = False


        self.pin_in_section = pin_in_section


        if self.pin_in_section == "iPiR":
            self.contact_model = self.contact_models[1]

            # self.slot_frame_nodes_GCS

        if self.pin_in_section in ["iP", "iR"]:
            self.contact_model = self.contact_models[0]

        if self._distance_obj is not None and self._contact_point_obj is None:
            # print "NO CONTACT OBJ"
            distance, delta, n_GCS, t_GCS = self._distance_obj.contact_geometry_GCS()

        else:
            # print "CONTACT OBJ"
            distance, delta, n_GCS, t_GCS = self._contact_point_obj.contact_geometry_GCS()

        # print "self.pin_in_section =", self.pin_in_section, "delta =", delta,

        # print "delta =", delta
        # print "n_GCS =", np.rad2deg(np.arctan2(n_GCS[1], n_GCS[0]))
        return distance, delta, n_GCS, t_GCS

    def _contact_geometry_LCS(self, q):
        """
        Function evaluates contact geometry parameters in body LCS based on contact geometry in GCS
        Function evaluates:
        self._n_LCS_list:   normal of contact in body LCS:
        self._t_LCS_list:   tangent of contact in body LCS
        self.u_P_LCS_list:  contact point in body LCS
        """
        self._contact_point_obj.contact_geometry_LCS()
        self._contact_point_obj.contact_points_LCS(q)

    def contact_velocity(self, q):
        """
        Function evaluates relative contact velocity vectors in normal and tangent direction
        """
        dq_n, dq_t = self._contact_point_obj_list[0].contact_velocity(q)
        # dr_P = []
        #
        # for i, (body_id, u_P) in enumerate(zip(self.body_id_list, self.u_P_LCS_list)):
        #     #   body velocity, R, theta
        #     dR = q2dR_i(q, body_id)
        #     theta = q2theta_i(q, body_id)
        #     #    dtheta - omega
        #     dtheta = q2dtheta_i(q, body_id)
        #     #    point velocity
        #     dr_P_body = dr_contact_point_uP(dR, theta, dtheta, u_P)
        #     #    add to list
        #     dr_P.append(dr_P_body)
        #
        # #    relative contact velocity vector
        # _dq = dr_P[1] - dr_P[0]
        # #   relative contact velocity
        # #   normal direction
        # dq_n = np.dot(_dq, self._n_GCS[0])
        #
        # #   tangent direction
        # dq_t = np.dot(_dq, self._t_GCS[0])
        #
        # #   set relative contact velocities to contact point object
        # self._contact_point_obj_list[0].set_dq(dq_n=dq_n, dq_t=dq_t)

        # self._contact_point_obj_list[0].contact_velocity(q)

        return dq_n, dq_t

    def set_dq0(self, dq_n, dq_t):
        """

        :return:
        """
        self._contact_point_obj_list[0].set_dq0(dq0_n=dq_n, dq0_t=dq_t)

    def _evaluate_delta_flat_section(self, distance):
        """
        Evaluate penetration depth - delta for contact between pin and flat surface of the slot
        :param distamce:
        :return:
        """
        delta = distance - self.R0_j
        return delta

    def _evaluate_delta_cylindrical_section(self, distance):
        """
        Evaluate penetration depth - delta for contact between pin and cylindrical surface of the slot
        :param distamce:
        :return:
        """
        delta = self._radial_clearance - distance
        return delta

    def _check_pin_position(self, frame_nodes_GCS, pin_CP_GCS):
        """
        Function eveluates where in slot is pin located to use correct equation for evaluation of penetration depth
        :param frame_nodes_GCS:     position of slot frame nodes in GCS
        :param pin_CP_GCS:          position of pin center in GCS
        :return:
        """
        #   set default values
        pin_in_section_jP = pin_in_section_jPjR = pin_in_section_jR = False

        #   distance of pin center to line DA and BC
        #   node A
        u_Ai = frame_nodes_GCS[0, :]
        #   node B
        u_Bi = frame_nodes_GCS[1, :]
        #   node C
        u_Ci = frame_nodes_GCS[2, :]
        #   node D
        u_Di = frame_nodes_GCS[3, :]
        #   DA
        # self._pin_distance_to_DA_obj = DistanceLineNode(pin_CP_GCS, u_Ai, r_jR=u_Di, parent=self)
        self._pin_distance_to_DA_obj = DistancePSCJ(u_Ai, u_Di, pin_CP_GCS, parent=self)
        # print "_pin_distance_to_DA_obj =", self._pin_distance_to_DA_obj._distance
        #   BC
        # self._pin_distance_to_BC_obj = DistanceLineNode(pin_CP_GCS, u_Ci, r_jR=u_Bi, parent=self)
        self._pin_distance_to_BC_obj = DistancePSCJ(u_Ci, u_Bi, pin_CP_GCS, parent=self)
        # print "_pin_distance_to_BC_obj =", self._pin_distance_to_BC_obj._distance
        # time.sleep(100)
        if self._pin_distance_to_DA_obj._distance < self.L and self._pin_distance_to_BC_obj._distance < self.L:
            pin_in_section_jPjR = True

        if self._pin_distance_to_DA_obj._distance < self.L and self._pin_distance_to_BC_obj._distance >= self.L:
            pin_in_section_jP = True

        if self._pin_distance_to_DA_obj._distance >= self.L and self._pin_distance_to_BC_obj._distance < self.L:
            pin_in_section_jR = True

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

    def evaluate_pin_position(self, q):
        """

        :param q:
        :return:
        """
        #   evaluate coordinates in GCS of mid frame section and pin center
        pin_CP_GCS, slot_CP_GCS, self.slot_frame_nodes_GCS, slot_frame_normals_GCS, slot_frame_tangents_GCS = self._contact_geometry_CP_GCS(q)

        #   vectors in GCS
        self.r_iP = slot_CP_GCS[0]
        self.r_iR = slot_CP_GCS[1]
        self.r_jP = pin_CP_GCS

        d = DistancePSCJ(self.r_iP, self.r_iR, self.r_jP, parent=self)

        d.evaluate_pin_in_section()

        print "Pin in section =", d.pin_in_section

    def evaluate_delta(self, q):
        """

        :param q:
        :return:
        """
        #   evaluate coordinates in GCS of mid frame section and pin center
        pin_CP_GCS, slot_CP_GCS, self.slot_frame_nodes_GCS, slot_frame_normals_GCS, slot_frame_tangents_GCS = self._contact_geometry_CP_GCS(q)

        #   vectors in GCS
        self.r_iP = slot_CP_GCS[0]
        self.r_iR = slot_CP_GCS[1]
        self.r_jP = pin_CP_GCS

        if self._contact_point_obj is not None and any(self._contact_point_obj._contact):
            self._contact_point_obj.update_contact_point(self.r_iP, self.r_iR, self.r_jP)

            delta = self._contact_point_obj._distance[self._contact_point_obj._contact.index(True)]

        else:
            delta = None
            print "No contact detected!"

        return delta

    def evaluate_CP_GCS(self, q):
        """

        :param q:
        :return:
        """
        rPj = self._pin_GCS(q)

        rPRi = self._slot_GCS(q)

        self.r_CP_GCS_list[2] = rPj

        self.r_CP_GCS_list[0:2] = rPRi

    def _contact_geometry_CP_GCS(self, q):
        """

        :param q:
        :return:
        """
        #   slot frame nodes in GCS
        frame_nodes_GCS, frame_normals_GCS, frame_tangents_GCS = self._slot_frame_flat_surface_GCS(q)

        #   pin center in GCS
        pin_CP_GCS = self._pin_GCS(q)

        #   slot center points in GCS
        slot_CP_GCS = self._slot_GCS(q)

        return pin_CP_GCS, slot_CP_GCS, frame_nodes_GCS, frame_normals_GCS, frame_tangents_GCS

    def _pin_GCS(self, q):
        """
        Function updates and calculates position of a pin center in global coordinates - GCS
        :param q:
        :return:
        """
        CP_GCS = u_P_lcs2gcs(self.u_jP, q, self.body_id_j)
        return CP_GCS

    def _slot_GCS(self, q):
        """

        :param q:
        :return:
        """
        slot_CP_GCS = []
        for CP in self.u_CP_slot_LCS:
            _uP_gcs = u_P_lcs2gcs(CP, q, self.body_id_i)
            slot_CP_GCS.append(_uP_gcs)
        return slot_CP_GCS

    def _create_slot_frame_flat_surface(self, uPi, uRi, n, t, hPi, hRi=None):
        """
        Direction of tangent uPj->uRj
        Direction of normal A->D for edge AB
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
        frame_nodes_LCS = np.array([uPi - n*hPi*.5,     #0   A
                                    uRi - n*hPi*.5,     #1   B
                                    uRi + n*hPi*.5,     #2   C
                                    uPi + n*hPi*.5],    #3   D
                                    dtype='float32')

        return frame_nodes_LCS

    def _slot_frame_flat_surface_GCS(self, q):
        """
        Function updates and calculates position of slot frame in global coordinates - GCS
        :param q:
        :return:
        """
        #   nodes in GCS
        frame_nodes_GCS = q2R_i(q, self.body_id_i) + Ai_ui_P_vector(self.frame_nodes_LCS, q2theta_i(q, self.body_id_i))

        #   normals in GCS
        frame_normals_GCS = Ai_ui_P_vector(np.array(self.n_edge_slot_LCS), q2theta_i(q, self.body_id_i))

        #   tangents in GCS
        frame_tangents_GCS = Ai_ui_P_vector(np.array(self.t_edge_slot_LCS), q2theta_i(q, self.body_id_i))

        return frame_nodes_GCS, frame_normals_GCS, frame_tangents_GCS

    def plot(self, q):
        """

        :param q:
        :return:
        """
        plt.clf()
        fig = plt.figure(num=1, figsize=(6, 5), dpi=72, facecolor='w', edgecolor='g')
        ax = plt.subplot(111, aspect="equal")

        pin_CP_GCS, slot_CP_GCS, frame_nodes_GCS, frame_normals_GCS, frame_tangents_GCS = self._contact_geometry_CP_GCS(q)

        #   pin - body j
        plt.plot(pin_CP_GCS[0], pin_CP_GCS[1], color='r')
        circle_j = plt.Circle((pin_CP_GCS[0], pin_CP_GCS[1]), self.R0_j, color='r', fill=False)
        ax.add_artist(circle_j)

        #   slot - body i
        plt.plot(np.append(frame_nodes_GCS[:, 0], frame_nodes_GCS[0, 0]), np.append(frame_nodes_GCS[:, 1], frame_nodes_GCS[0, 1]), color='b')
        circle_iP = plt.Circle((slot_CP_GCS[0][0], slot_CP_GCS[0][1]), self.R0_i, color='b', fill=False)
        ax.add_artist(circle_iP)
        circle_iR = plt.Circle((slot_CP_GCS[1][0], slot_CP_GCS[1][1]), self.R0_i, color='b', fill=False)
        ax.add_artist(circle_iR)

        # slot_CP_GCS, self.slot_frame_nodes_GCS

        # pin_circle = plt.Circle(pin_CP_GCS, self.R0_j, color='r', fill=False)
        # ax.add_artist(pin_circle)
        # plt.plot(pin_CP_GCS[0] + x_i, pin_CP_GCS[1] + y_i)

        # iP_circle = plt.Circle(slot_CP_GCS[0], self.h0_iP/2., color='b', fill=False)
        # ax.add_artist(iP_circle)
        # plt.plot(slot_CP_GCS[0][0] + x_j, slot_CP_GCS[0][1] + y_j, color='b')
        # plt.plot(slot_CP_GCS[0][1] + x_j, slot_CP_GCS[1][1] + y_j, color='b')
        # iR_circle = plt.Circle(slot_CP_GCS[1], self.h0_iP/2., color='b', fill=False)
        # ax.add_artist(iR_circle)
        plt.show()

    def print_r_CP_list(self, q):
        """
        Print list of center points - CP of pin slot clearance joint
        :return:
        """
        self.evaluate_CP_GCS(q)
        print "List of Center Points"

        for index, point_name, r_CP in zip(["i", "i", "j"], ["P", "R", "P"], self.r_CP_GCS_list):
            print "r_" + point_name + "_" + index + "=", r_CP

    def print_r_slot_frame(self, q):
        """

        :param q:
        :return:
        """
        print "List of Slot Frame Points"

        frame_nodes_GCS, frame_normals_GCS, frame_tangents_GCS = self._slot_frame_flat_surface_GCS(q)

        for i, rP in enumerate(frame_nodes_GCS):
            print i, rP

    def testing(self):
        """
        Only for testing
        :return:
        """
        print "TESTING@",__name__

        for contact_point_obj in self._contact_point_obj_list:
            for Fn, Ft in zip(contact_point_obj._Fn_list, contact_point_obj._Ft_list):
                Fn.active = False
                Ft.active = False

if __name__ == "__main__":
    _type = "pin-slot clearance joint linear"
    body_id_i = 0
    body_id_j = 1
    u_iP = np.array([0, 0])*1E-3
    u_iR = np.array([2, 0])*1E-3
    u_jP = np.array([0, 0])*1E-3
    R0_j = 2.45E-3
    h0_iP = 5E-3
    _pscj = PinSlotClearanceJointLinear(_type, body_id_i, body_id_j, u_iP, u_iR, u_jP, h0_iP, R0_j, parent=None)
    # pprint(vars(_joint))

    pin_CP_GCS = np.array([5.200000014E-04, 3.999999899E-05])
    frame_nodes_GCS = np.array([[4.606226617e+00, -2.417889370e+00],
                                [4.607809150e+00, -2.419112314e+00],
                                [4.604751791e+00, -2.423068644e+00],
                                [4.603169258e+00, -2.421845701e+00]])

    print _pscj._check_pin_position(frame_nodes_GCS, pin_CP_GCS)

    fig = plt.figure()
    ax = fig.add_subplot(111,aspect='equal')
    pin = plt.Circle(pin_CP_GCS,R0_j, color='b', fill=False)

    fig.gca().add_artist(pin)
    plt.plot(pin_CP_GCS[0], pin_CP_GCS[1],
             "x",
             color='b')

    plt.plot(np.append(frame_nodes_GCS[:, 0], frame_nodes_GCS[0,0]),
             np.append(frame_nodes_GCS[:, 1], frame_nodes_GCS[0,1]),
             linestyle="-",
             color="red")
    plt.show()
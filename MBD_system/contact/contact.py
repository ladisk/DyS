"""
Created on 21. feb. 2014

@author: lskrinjar
"""
import time
import copy
from __builtin__ import enumerate
from operator import attrgetter
from pprint import pprint
from matplotlib import pyplot as plt
from scipy.spatial import _distance_wrap


from MBD_system.MBD_system_items import ContactItem
from MBD_system.check_filename import check_filename
from MBD_system.contact.bounding_box.bounding_box import AABB
from MBD_system.contact.bounding_box.bounding_box_2D import AABB2D
from MBD_system.contact.contact_geometry.contact_geometry import ContactGeometry
from MBD_system.contact.evaluate_distance import evaluate_distance_2D
from MBD_system.contact.overlap_pair.overlap_pair import OverlapPair
from MBD_system.contact_model.contact_model import ContactModel
from MBD_system.contact_model.contact_model_cylinder import ContactModelCylinder
from MBD_system.dr_contact_point_uP import dr_contact_point_uP
from MBD_system.extract_from_dictionary_by_string_in_key import extract_from_dictionary_by_string_in_key
from MBD_system.force.force import Force
from MBD_system.friction_model.friction_model import FrictionModel
from MBD_system.q2dR_i import q2dR_i
from MBD_system.q2dtheta_i import q2dtheta_i
from MBD_system.q2q_body import q2q_body
from MBD_system.q2theta_i import q2theta_i
from MBD_system.solution_data.solution_data_contact import SolutionDataContact
from MBD_system.transform_cs import cm_lcs2gcs, gcs2cm_lcs
from contact_C_matrix import *


class Contact(ContactItem):
    """
    Contact object to solve contacts between a pair of bodies, ie. body i and body j
    Main logic of contact class:

    class method _contact_geometry_GCS() is evaluated at every time step (contact detection) this
    method is run in class method check_for_contact()

    class method _contact_geometry_LCS() is evaluated only when contact is detected

    class method _get_contact_geometry_data() is evaluated only when contact is detected and before class method _contact_geometry_LCS()

    class attribute _delta is depth of deformation at contact point and is:
    negative (delta) - impact
    positive (delta) - no impact

    relative normal contact velocity (dq_n)
    positive - bodies are approaching
    negative - bodies are separating
    """
    __id = itertools.count(0)

    def __init__(self, _type, body_id_i, body_id_j, name=None, contact_model_type=None, friction_model_type=None, properties_dict={}, parent=None):
        """
        Constructor of class contact
        :param  _type               supported types:
                                    1 general
                                    2 sphere-sphere
                                    3 plane-sphere
                                    4 revolute clearance joint
                                    5 pin-slot (linear) clearance joint
                                    6 prismatic (translational) clearance joint
        :param  body_id_i           body i id
        :param  body_id_j           body j id
        :param  name                name of the contact as str
        :param  contact_model_type  type of contact model (str)
        :param  friction_model_type type of friction model (str)
        :param  properties_dict     dictionary of additional object properties
        :param  parent              pointer to the parent object
        """
        #    number
        self.contact_id = self.__id.next()

        #    name as string
        if name is None:
            self._name = "Contact_" + str(self.contact_id)
        else:
            self._name = name + str(self.contact_id)

        # this has to be after attributes contact_id and _name as name is constructed from contact_id and _name
        super(Contact, self).__init__(self._name, parent)

        #   parent
        self._parent = parent

        #   bodies
        if self._parent is not None:
            self._bodies = self._parent._parent.bodies
        else:
            self._bodies = []

        #   MBD item type
        self._type = "contact"

        #    contact type
        self._contact_type = _type

        #   supported types
        self._contact_types = ["general",
                               "revolute clearance joint",
                               "prismatic clearance joint",
                               "contact sphere-sphere",
                               "contact plane-sphere",
                               "pin-slot clearance joint linear",
                               "pin-slot clearance joint radial",
                               "contact point-line",
                               "roughness profile"]

        #   contact geometry (profile) - list of object pairs
        self.contact_geometry_list = [None, None]

        #   position of contact in z direction
        self.z_dim = 0

        #   vector of initial positions
        self.q0 = self._parent._parent.evaluate_q0()

        #   time dependent variables
        self._step = 0
        self._t = 0.
        self._distance = None
        self._sign_check = None

        #   solution attributes
        #   save solution - properties
        self.save_to_file = False
        self._solution_filename = None
        #   a solutions list to store solution data of contact object of each simulation
        self.solutions = []
        self.solution_data = None
        #   solution options:
        #   Discard
        #   Overwrite (existing file)
        #   Save to new (next available) file
        self._solution_save_options = "Discard"
        #   solution file type, optional: .dat, .xlsx, .csv
        self._solution_filetype = ".xlsx"
        self._solution_filetypes = [".sol", ".dat", ".xlsx", ".csv"]

        #   visualization properties (default)
        self.color = np.array([0.4, 0.2, 0.2], dtype="float32")
        self.scale = 1.
        self._visible_F_list = [True, True]

        #   markers
        self.markers = []

        #    contact body ids
        self.body_id_i = body_id_i
        self.body_id_j = body_id_j
        self.body_id_list = [self.body_id_i, self.body_id_j]
        self.body_id_edge = None
        self.body_id_node = None
        self.body_ids_edge_list = []
        self.body_ids_node_list = []

        self.body_list = [None, None]
        if self._bodies != []:
            self.body_list = [self._bodies[self.body_id_i], self._bodies[self.body_id_j]]
        self.contact_geometry_list = [None, None]
        #   type of contact geometry of each body
        self.body_i_contact_geometry_type = ""
        self.body_j_contact_geometry_type = ""

        #    set properties dictionary as object property
        self.properties = copy.copy(properties_dict)
        #   check if extra properties are defined in dict and if not, assign default values
        if self._parent is not None:
            _body_i_contact_geometry_type = "body_i.contact_geometry_type"
            if _body_i_contact_geometry_type not in self.properties:
                self.properties[_body_i_contact_geometry_type] = None

            _body_j_contact_geometry_type = "body_j.contact_geometry_type"
            if _body_j_contact_geometry_type not in self.properties:
                self.properties[_body_j_contact_geometry_type] = None

        #   general contact default properties
        self.AABB_list_of_overlap_pairs = []
        self.AABB_list = []

        #    status of contact
        #    0 - no contact
        #    -1 - contact is detected, two bodies have already collided and are in collision, but the
        #    current initial penetration depth is too large (halve the time step)
        #    1 - contact is detected
        self.status = 0
        #   numerical error when calculating distance (point to line), in m
        self.distance_TOL = 1.E-7

        #   contact point(s) in GCS
        self.r_P_GCS = np.array(np.zeros(2))
        self.r_iP_GCS = np.zeros(2, dtype="float32")
        self.r_jP_GCS = np.zeros(2, dtype="float32")
        self.r_P_GCS_list = [self.r_iP_GCS, self.r_jP_GCS]
        self.r_P_GCS_0_list = copy.copy(self.r_P_GCS_list)

        #   contact point in body LCS
        self.u_P_LCS = None
        self.u_iP_LCS = np.zeros(2)
        self.u_jP_LCS = np.zeros(2)
        self.u_P_LCS_list = [self.u_iP_LCS, self.u_jP_LCS]
        self.u_P_LCS_0_list = copy.copy(self.u_P_LCS_list)

        #   list of normals and tangents in LCS of each body in contact
        #   normals
        self._n_GCS = np.zeros(2)
        self._n_GCS_MCP_list = []
        self._n_i_LCS = np.zeros(2)
        self._n_j_LCS = np.zeros(2)
        self._n_LCS_list = [self._n_i_LCS, self._n_j_LCS]
        self._n_GCS_list = []
        #   tangents
        self._t_GCS = np.zeros(2)
        self._t_GCS_MCP_list = []
        self._t_i_LCS = np.zeros(2)
        self._t_j_LCS = np.zeros(2)
        self._t_LCS_list = [self._t_i_LCS, self._t_j_LCS]
        self._t_GCS_list = []

        #   contact forces
        self.contact_bodies_added_to_list = False
        self.list_of_contact_force_objects_constructed = False
        #   create normal contact forces
        self._Fn_list = []#self._create_Fn_forces()
        #   create tangential (friction) contact forces
        self._Ft_list = []#self._create_Ft_forces()

        #   contact force (type:ndarray)
        #   force consists of normal and tangent force vector of contact
        self.F = np.zeros(2)
        self.Fn = 0.
        self.Ft = 0.

        #   list of generalized external forces, pair of generalized force on each body
        self.Q_e_list = []

        #   relative contact velocity
        #   define relative contact velocity components as object attributes
        self._dq_n = 0.
        self._dq_t = 0.
        self.initial_contact_velocity_calculated = False

        #   create solution containers
        self._solution_containers()

        #   add additional properties
        if self.properties != {}:
            self.add_additional_parameters(self.properties)

        #    contact distance - penetration depth
        self._delta = copy.copy(self.distance_TOL)
        self._delta0 = copy.copy(self.distance_TOL)
        self._delta_n = copy.copy(self.distance_TOL)
        self.contact_detected = False  # True or False
        self.contact_distance_inside_tolerance = False
        self._contact_point_found = False
        self._distance_obj = None
        self._contact_point_obj = None
        self._distance_obj_list = []
        self._distance_list = []
        #   list of actual contact point objects at time
        self._contact_point_obj_list = []

        #   contact node-edge points to calculate contact properties
        self.node_LCS = None
        self.node_GCS = None
        self.edge_LCS = [None, None]
        self.edge_GCS = [None, None]

    def set_vtk_data(self, interactor=None):
        """
        Function is overridden in subclass
        :return:
        """

    def update_vtk_data(self, q):
        """
        Function is overriden in subclass
        :param q:
        :return:
        """

    def _create_contact_model(self, properties_dict):
        """
        Function creates contact model object that is used to evaluate normal contact force at contact
        """
        #   default properties
        self.contact_model_type = None
        self.contact_model_types = []

        #   number of contact models defined
        _n = self._count_contact_models(properties_dict)
        self.contact_models = []

        self.properties_contact_model = extract_from_dictionary_by_string_in_key(properties_dict, "contact_model.")

        if self.contact_model_type is None:
            self.contact_model_type = "hertz"

        if "type" in self.properties_contact_model:
            self.contact_model_type = self.properties_contact_model["type"]

        if _n == 1 or _n == 0:
            _properties = extract_from_dictionary_by_string_in_key(properties_dict, "contact_model.")
            if _properties["type"] in ContactModel._supported_types():
                self.contact_model = ContactModel(self.contact_model_type, properties_dict=_properties, parent=self)
            elif _properties["type"] in ContactModelCylinder._supported_types():
                self.contact_model = ContactModelCylinder(self.contact_model_type, properties_dict=_properties, parent=self)
            else:
                self.contact_model = ContactModel(self.contact_model_type, properties_dict=_properties, parent=self)
                raise Warning, "Contact model type not correct: %s", _properties["_type"]
        else:
            self.contact_model = None
            for i in range(0, _n):
                _substring = "contact_model" + "[" + str(i) + "]" + "."
                _properties = extract_from_dictionary_by_string_in_key(properties_dict, _substring)

                if _properties["_type"] in ContactModel._supported_types():
                    _contact_model = ContactModel(_properties["_type"], properties_dict=_properties, _substring=_substring, parent=self)
                else:
                    _contact_model = ContactModelCylinder(_properties["_type"], properties_dict=_properties, _substring=_substring, parent=self)

                self.contact_models.append(_contact_model)

    def _create_friction_model(self, properties_dict):
        """
        Function creates friction model object that is used to evaluate friction force
        """
        #   default properties
        self.friction_model_type = None
        self.friction_model_types = []

        #    friction model
        self._dq_t_TOL = 1E-3
        self.coef_of_friction_dynamic = 0
        self.coef_of_friction_static = 0
        self.properties_friction_model = extract_from_dictionary_by_string_in_key(properties_dict, "friction_model.")

        if self.friction_model_type is None:
            self.friction_model_type = "ideal"
        if "type" in self.properties_friction_model:
            self.friction_model_type = self.properties_friction_model["type"]

        self.friction_model = FrictionModel(self.friction_model_type, coef_of_friction_dynamic=self.coef_of_friction_dynamic, coef_of_friction_static=self.coef_of_friction_static,
                                            properties_dict=self.properties_friction_model, parent=self)

    def _create_Fn_forces(self):
        """
        Function creates list of contact forces as objects
        :return:
        """
        self.Fn = 0.
        Fn_list = []
        if self._parent is not None:
            for body_id in self.body_id_list:
                _Fn = Force(body_id, force_name=self._name + "_Fn_on_body_" + str(body_id), parent=self)
                #   add pair of contact forces to forces list of MBD system
                self._parent._parent.forces.append(_Fn)
                #   add pair of contact forces to forces list of contact
                Fn_list.append(_Fn)
                self._parent._parent.bodies[body_id].forces.append(_Fn)

                #    set color
                _Fn._color = self._parent._parent.bodies[body_id].color

            Fn_list[0]._visible = False
            self.list_of_contact_force_objects_constructed = True

        return Fn_list

    def _create_Ft_forces(self):
        """
        Function creates list of friction (tangential) contact forces as objects
        :return:
        """
        self.Ft = 0.
        Ft_list = []
        if self._parent is not None:
            for body_id in self.body_id_list:
                _Ft = Force(body_id, force_name=self._name + "_Ft_on_body_" + str(body_id), parent=self)
                #   add pair of contact forces to forces list of MBD system
                self._parent._parent.forces.append(_Ft)
                #   add pair of contact forces to forces list of contact
                Ft_list.append(_Ft)
                self._parent._parent.bodies[body_id].forces.append(_Ft)

                #    set color
                _Ft._color = self._parent._parent.bodies[body_id].color

            Ft_list[0]._visible = False
            self.list_of_contact_force_objects_constructed = True

        return Ft_list

    def change_visibility_Fi(self):
        """

        :return:
        """
        indx = 0
        self._change_F_visibility(indx)

    def change_visibility_Fj(self):
        """

        :return:
        """
        indx = 1
        self._change_F_visibility(indx)

    def _change_F_visibility(self, indx):
        """

        :return:
        """
        for contact_point in self._contact_point_obj_list:
            for F in [contact_point._Fn_list[indx], contact_point._Ft_list[indx]]:
                if hasattr(F.vtk_actor, "GetVisibility"):
                    if F.vtk_actor.GetVisibility():
                        F.vtk_actor.VisibilityOff()
                        self._visible_F_list[indx] = False
                    else:
                        F.vtk_actor.VisibilityOn()
                        self._visible_F_list[indx] = True

    def add_additional_parameters(self, d):
        """
        Add additional attributes from dictionary self.properties
        """
        _dict = copy.copy(d)
        #    set object properties from dictionary
        for key in _dict:
            #    set Contact object additional properties
            if "body" not in key and "contact " not in key:
                if key == "name":
                    val = _dict[key]
                    key = "_name"
                    _dict[key] = val
                setattr(self, key, _dict[key])

                # _dict.pop(key, None)

            if "[" in key and "]" in key and "model" not in key:
                _attr = key[0:key.find("[")]
                if hasattr(self, _attr):
                    #   index of item (object) in a list
                    _index = int(key[key.find("[") + 1:key.find("]")])

                    #   name of the property (attribute)
                    _prop = key[key.find("]") + 2:]

                    #   get object from list by index if list is not empty
                    if not _attr:
                        _obj = getattr(self, _attr)[_index]

                        #   set attribute of object
                        setattr(_obj, _prop, _dict.get(key))

            # remove key from dictionary after added to object attributes
            elif "." not in key:
                pass
                # _dict.pop(key, None)
            else:
                None

            # _dict.pop(key, None)

    def _count_contact_models(self, _dict):
        """
        Method counts number of "contact_model._type" string in _dict keys
        :param _dict:
        :return:
        """
        _substring_1 = "contact_model"
        _substring_2 = "type"

        counter = 0
        _counter = itertools.count(1)
        for key in _dict:
            if _substring_1 in key and _substring_2 in key:
                counter = _counter.next()

        return counter

    def _solution_containers(self):
        """

        :return:
        """
        self.solution_data = SolutionDataContact(MBD_item=self)
        #    solution containers over integration times
        self._step_num_solution_container = [0]

        self._status_container = [0]

        self._step_size_solution_container = [0]

        self._t_solution_container = [0]

        self._distance_solution_container = [self.distance_TOL]

        self._dqn_solution_container = [0]

        self._dqt_solution_container = [0]

        self._Fn_solution_container = [0]

        self._Ft_solution_container = [0]

        self._Fx_solution_container = [0]

        self._Fy_solution_container = [0]

        self._u_P_solution_container = [[np.zeros(2), np.zeros(2)]]

        self._r_P_solution_container = [[np.zeros(2), np.zeros(2)]]

        self.__solution_containers_list = [self._step_num_solution_container,
                                           self._status_container,
                                           self._step_size_solution_container,
                                           self._t_solution_container,
                                           self._distance_solution_container,
                                           self._dqn_solution_container,
                                           self._dqt_solution_container,
                                           self._Fn_solution_container,
                                           self._Ft_solution_container]

        self.__solution_container_names_list = ["_step_num_solution_container",
                                                "_status_container",
                                                "_step_size_solution_container",
                                                "_t_solution_container",
                                                "_distance_solution_container",
                                                "_dqn_solution_container",
                                                "_dqt_solution_container",
                                                "_Fn_solution_container",
                                                "_Ft_solution_container"]

        self.__solution_container_dict = {"delta": self._distance_solution_container,
                                          "Fn": self._Fn_solution_container,
                                          "Ft": self._Ft_solution_container,
                                          "dqn": self._dqn_solution_container,
                                          "dqt": self._dqt_solution_container}

    def data_tracker(self, t, step):
        """
        Track contact object data during simulation. Track:
        simulation time
        step number
        step size
        """
        self.t = t
        self._step = step

        #    append last value if time t is greater than last element in data container
        if t > self._t_solution_container[-1]:
            self._step_solution_accepted = True

        # replace last value if new value is smaller than last value
        else:
            self._step_solution_accepted = False

    def _track_data(self, step, h, t):
        """
        Function saves current calculated data at
        :return:
        """
        self._step = step
        self._h = h
        self._t = t

        self.solution_data._track_data()
        #   step number
        self._step_num_solution_container.append(step)

        #   contact status
        self._status_container.append(self.status)

        #   step size
        self._step_size_solution_container.append(h)

        #   time
        self.t = t
        self._t_solution_container.append(t)

        #   delta
        if type(self._delta) is list:
            self._distance_solution_container.extend(self._delta)
        else:
            self._distance_solution_container.append(self._delta)

        #   dq
        self._dqn_solution_container.append(self._dq_n)
        self._dqt_solution_container.append(self._dq_t)

        #   F
        #   normal and tangent component of contact force
        self._Fn_solution_container.append(self.Fn)
        self._Ft_solution_container.append(self.Ft)
        self._Fn_solution_container.append(self.F)

        #   uPi, uPj
        self._u_P_solution_container.append(self.u_P_LCS_list)

        #   rPi, rPj
        self._r_P_solution_container.append(self.r_P_GCS_list)

    def evaluate_Q_e(self, t, q):
        """
        Function evaluates vector of generalized external force on each body in contact
        :param t:
        :param q:
        :return:
        """
        self.Q_e_list = []
        for Fn, Ft in zip(self._Fn_list, self._Ft_list):
            Q_e = Fn.evaluate_Q_e(t, q) + Ft.evaluate_Q_e(t, q)
            self.Q_e_list.append(Q_e)

        [Q_e_i, Q_e_j] = self.Q_e_list
        return Q_e_i, Q_e_j

    def _reset(self):
        """
        Defined in subclass
        :return:
        """

    def reset(self):
        """
        Reset main object attributes to initial value
        """
        self._contact_point_found = False
        self.initial_contact_velocity_calculated = False
        self.status = 0

        self._solution_containers()

        self.solution_data.reset()

        #    reset force object solution data containers
        for force in self._Fn_list:
            force.reset()

        #   reset action of subclass variables and properties
        self._reset()

    def contact_update(self, step, t, q):
        """
        Function updates status only for (general) contact
        :return:
            status - status of contact 0-no contact, 2-contact
        """
        # print "contact_update()"
        self._step = step

        #   calculate relative contact velocity in normal direction at time t
        self._dq_n, self._dq_t = self.contact_velocity(q)

        #   calculate contact geometry in GCS to calculate contact distance
        # self._contact_geometry_GCS(q)

        #   calculate distance: node-edge
        _distance, _inside = evaluate_distance_2D(self.node_GCS, self.edge_GCS[0], self._n_GCS, self._t_GCS)

        #   contact is present
        if _inside and self._distance < self.distance_TOL:
            self._distance = -_distance
            self.status = 1

        # no contact
        else:
            self._distance = _distance
            self.status = 0
            self.contact_detected = False
            self.contact_distance_inside_tolerance = False
            self._contact_point_found = False
            self.initial_contact_velocity_calculated = False

            for _Fn, _Ft in zip(self._Fn_list, self._Ft_list):
                _Fn._update_force_vector(np.zeros(2))
                _Ft._update_force_vector(np.zeros(2))

            self._no_contact()

        # print "self.status - contact_update() =", self.status
        return self.status

    def evaluate_contact(self, t, q):
        """
        Calculate contact and evaluate:
        1. contact penetration depth
        2. contact forces in normal and tangent direction
        returns None:
        """
        # print "evaluate_contact()"
        #    calculate coordinates of contact point from global coordinates in local coordinates of each body in contact
        if not self._contact_point_found:
            self._get_contact_geometry_data(q)
            self._contact_point_found = True

        self._distance, self._delta, self._n_GCS, self._t_GCS = self._contact_geometry_GCS(q)

        self._contact_geometry_LCS(q)

        #   kinematic properties of contact point
        self._dq_n, self._dq_t = self.contact_velocity(q)

        #   initial contact velocity
        if not self.initial_contact_velocity_calculated or self.section_changed:
            self._dq0_n, self._dq0_t = self.contact_velocity(q)
            self.contact_model.set_dq0(self._dq0_n, self._dq0_t)

            for contact_model in self.contact_models:
                contact_model.set_dq0(self._dq0_n, self._dq0_t)

            for contact_point in self._contact_point_obj_list:
                contact_point.set_dq0()

            self.initial_contact_velocity_calculated = True

        if self._contact_type.lower() in self._contact_types:  # == "general" or self._contact_type.lower() == "revolute clearance joint" or self._contact_type.lower() == "contact sphere-sphere":#ECF-N
            # print "self._dq_n =", self._dq_n
            self._solve_ECF_N(t, q, self._delta, self._dq_n, self._dq_t)  # self._delta
        else:
            raise ValueError, "Contact type not correct!"

    def check_for_contact_continued_condition(self, delta, dq_n, q):
        """
        Method must be overriden in subclass
        :param delta:
        :param dq_n:
        :return:
        """

    def set_dq0(self, dq_n, dq_t):
        """

        :return:
        """
        for contact_point_obj in self._contact_point_obj_list:
            contact_point_obj.set_dq0(dq0_n=dq_n, dq0_t=dq_t)

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
        if self.check_for_contact_continued_condition(_delta, _dq_n, q):
            for i, (contact_point, delta) in enumerate(zip(self._contact_point_obj_list, _delta)):
                #   normal contact force
                Fn = self.contact_model.contact_force(delta, contact_point._dq_n, dq0_n=contact_point._dq0_n)
                #   tangent contact force
                Ft = self.friction_model.friction_force(contact_point.Fn, contact_point._dq_t)
                # print contact_point, "Fn =", Fn, "Ft =", Ft, "_delta =", _delta
                contact_point.evaluate_F(Fn, Ft)

        #   no contact
        else:
            print "CONTACT FINISHED - self._step =", self._step
            print "delta =", _delta

            self._contact_point_found = False
            self.initial_contact_velocity_calculated = False
            self.contact_distance_inside_tolerance = False
            self.contact_detected = False
            self.status = 0

            # fig = plt.figure(num=1, figsize=(6, 4), dpi=100, facecolor='w', edgecolor='k')
            # ax = plt.subplot(111, aspect="equal")
            # for i, contact_point in enumerate(self._contact_point_obj_list):
            #     contact_point.plot()
            # filename = "step_" + str(self._step).zfill(4) + "_contact_finished" + ".png"
            # plt.savefig(filename)
            # print "Plot saved to file: ", filename
            # plt.clf()

            #   delete force objects from contact force lists (normal and tangent direction)
            self.no_contact()

        #   update contact forces
        self._update_contact_forces(q)

    def _update_contact_forces(self, q):

        #   create, update or delete force object of contact forces at contact points
        for i, contact_point in enumerate(self._contact_point_obj_list):
            if contact_point._contact_point_found and contact_point.active:
                for body_id, Fn_i, Ft_i, u_P_i, n_i, t_i in zip(contact_point.body_id_list, contact_point._Fn_list, contact_point._Ft_list, contact_point.u_P_LCS_list, contact_point._n_GCS_list, contact_point._t_GCS_list):

                    #   print only contact point of body j - pin
                    #   normal force
                    Fn_i.update_F(q=q, step=self._step, F=contact_point.Fn * n_i, u_P=u_P_i)
                    Fn_i._visible = True
                    #   tangent force
                    Ft_i.update_F(q=q, step=self._step, F=contact_point.Ft * t_i, u_P=u_P_i)
                    Ft_i._visible = True

    def _no_contact(self):
        """
        Defined in a subclass
        :return:
        """

    def no_contact(self):
        """
        Method sets values to contact attributes when there is no contact
        """
        self._no_contact()

        self._contact_point_found = False
        # self._distance_obj = None
        self.Fn = 0.
        self.Ft = 0.
        self.F = np.zeros(2, dtype="float32")

        if self._delta >= self._delta0:
            self._delta0 = None

        # self._delta = self.distance_TOL

        self.initial_contact_velocity_calculated = False

        #   uPi and uPj - position vector of contact force in body LCS
        self.u_P_LCS_list = [np.zeros(2), np.zeros(2)]
        self.r_P_GCS_list = [np.zeros(2), np.zeros(2)]

        self._dq_n = 0.
        self._dq_t = 0.

        #   update contact force values with zeros as no contact is present
        for _force_n, _force_t in zip(self._Fn_list, self._Ft_list):
            _force_n.update_step(self._step)
            _force_t.update_step(self._step)

        for _distance_obj in self._distance_obj_list:
            for Fn, Ft in zip(_distance_obj._Fn_list, _distance_obj._Ft_list):
                Fn.update_step(self._step)
                Ft.update_step(self._step)
            # self._distance_obj_list.remove(_distance_obj)

        # for contact_point_obj in self._contact_point_obj_list:
        #     contact_point_obj.remove_forces()

    def no_overlap(self, step, t, q):
        """
        Function runs at time step when there is no overlap of pair of (main) AABBs
        """
        #   step
        self._step = step
        #   time
        self._t = t

        if len(self._t_solution_container) == 1:
            pass
        else:
            self._delta = self.distance_TOL

        self.no_contact()
        # self._status_container = np.append(self._status_container, self.status)

    def evaluate_potential_energy(self, q):
        """

        :return:
        """
        _energy = []
        for contact_point in self._contact_point_obj_list:
            E_i = 0.5 * self.contact_model.K * contact_point._distance_sign ** 2
            _energy.append(E_i)

        if _energy != []:
            return sum(_energy)
        else:
            return 0.

    def create_bounding_box_for_each_body_in_contact(self, bodies):
        """
        Create bounding box object for each body
        :param bodies:  list of bodies of MBD system ojbect
        """
        _body_name_list = ["body_i", "body_j"]

        for body_id, _body, contact_geometry in zip(self.body_id_list, _body_name_list, self.contact_geometry_list):
            #    assign body object from body list to new variable name
            body = bodies[body_id]
            setattr(body, "skin_thickness", 10 * self.distance_TOL)

            #   remove duplicate nodes
            # _nodes, _normals = body_.geom.geom_data.remove_duplicate_2D()

            #   build contact geometry object and append it to lists (in contact and in body)
            # _contact_area = self.properties[_body+".contact_area"]
            if contact_geometry is None:
                contact_geometry_type = self.properties[_body + ".contact_geometry_type"]

                contact_geometry = ContactGeometry(contact_geometry_type, body, parent=self)

                self.contact_geometry_list.append(contact_geometry)

            # add contact geometry object to body list of contact geometries
            body.contact_geometries.append(contact_geometry)

            #    get nodes and normals from contact geometry object
            _nodes = contact_geometry.get_vertices_2D()
            _normals = contact_geometry.get_normals_2D()
            _tangents = contact_geometry.get_tangents_2D()
            _angles = contact_geometry.get_angles_2D()

            #    if value is assigned to position_of_plane_in_z_direction only 2D nodes are put in argument
            if self.z_dim is not None:
                _nodes = _nodes  # [:, 0:2]
                _normals = _normals  # [:, 0:2]

            # add additional parameter to body object - body point that has maximum distance
            #    from center of body mass
            body.uP_i_max = np.linalg.norm(np.array(np.amax(_nodes, axis=0)), ord=2)

            #    head AABBTree object is constructed
            body.AABBtree = AABB2D(nodes_LCS=_nodes, normals_LCS=_normals, tangents_LCS=_tangents, angles_LCS=_angles, parent=self, parent_body=body)
            # print "body_.AABBtree =", body_.AABBtree
            # pprint(vars(body_.AABBtree))
            body.AABBtree._visible = True

            #    if additional properties are defined, than they are assigned to AABB objects
            for key in self.properties:
                if body_id in _body_name_list:
                    if body_id in key and "." in key:
                        setattr(body.AABBtree, key[key.rindex('.') + 1:], self.properties[key])

            #    create AABBTree object with recursion
            body.AABBtree.construct()

            #    append AABB object to list
            self.AABB_list.append(body.AABBtree)

        # after the list of AABBtree for the contact is created (filled with object-data) each AABBtree
        #    is saved as different attribute for later calculations
        # print "self.AABB_list =", self.AABB_list
        if self.body_id_list is not []:
            self.AABB_i, self.AABB_j = self.AABB_list

        self.add_additional_parameters(self.properties)

    def _AABB_AABB_overlap(self, AABB_i, AABB_j):
        """
        Function used to recursively check every pair of AABB object that has no children
        link:
            http://gamemath.com/2011/09/detecting-whether-two-boxes-overlap/
        """
        #    if parameter is defined (not None) overlap of AABB is checked in 2D else in 3D
        #   2D overlap check
        if self.z_dim is not None:
            if AABB_i.x_max_GCS < AABB_j.x_min_GCS:  # 1 is left of 2
                return False
            if AABB_i.x_min_GCS > AABB_j.x_max_GCS:  # 1 is right of 2
                return False
            if AABB_i.y_max_GCS < AABB_j.y_min_GCS:  # 1 is above of 2
                return False
            if AABB_i.y_min_GCS > AABB_j.y_max_GCS:  # 1 is below of 2
                return False
            return True
        # 3D overlap check
        else:
            if AABB_i.x_max_GCS < AABB_j.x_min_GCS:  # 1 is left of 2
                return False
            if AABB_i.x_min_GCS > AABB_j.x_max_GCS:  # 1 is right of 2
                return False
            if AABB_i.y_max_GCS < AABB_j.y_min_GCS:  # 1 is above of 2
                return False
            if AABB_i.y_min_GCS > AABB_j.y_max_GCS:  # 1 is below of 2
                return False
            if AABB_i.z_max_GCS < AABB_j.z_min_GCS:  # 1 is in front of 2
                return False
            if AABB_i.z_min_GCS > AABB_j.z_max_GCS:  # 1 is behind of 2
                return False
            return True

    def check_for_overlap(self, q, AABB_i, AABB_j):
        """
        Function recursively checks if AABB_i from body i and AABB_j from body j overlap
        if they overlap and they have no children or their children do not overlap the two AABBs 
        are stored to new OverlapPair() object
        """
        self.AABB_update_list = [AABB_i, AABB_j]

        #    update frame geometry of main frame of every body
        for AABB, body_id in zip(self.AABB_update_list, self.body_id_list):
            q_i = q2q_body(q, body_id)
            AABB.frame_geometry_GCS(q_i)

        # check for overlap of each AABB frame
        if self._AABB_AABB_overlap(AABB_i, AABB_j):
            #    check between children of AABB_i and AABB_j
            if AABB_i.has_children() and AABB_j.has_children():
                for __AABB_i in AABB_i.children:
                    for __AABB_j in AABB_j.children:  # self.AABB_j.children
                        self.check_for_overlap(q, __AABB_i, __AABB_j)

            # check between AABB_j and children of AABB_i
            elif AABB_i.has_children() and not AABB_j.has_children():
                for __AABB_i in AABB_i.children:
                    self.check_for_overlap(q, __AABB_i, AABB_j)

            # check between AABB_i and children of AABB_j
            elif not AABB_i.has_children() and AABB_j.has_children():
                for __AABB_j in AABB_j.children:
                    self.check_for_overlap(q, AABB_i, __AABB_j)

            # if none of the two AABB has children, the two AABB are used to create contact pair object
            else:
                __AABB_overlap_pair = OverlapPair(_AABB_i=AABB_i, _AABB_j=AABB_j, parent=self)
                #    add created object to list of all AABB overlap pairs
                self.AABB_list_of_overlap_pairs.append(__AABB_overlap_pair)

    def check_for_contact(self, step, t, q):
        """
        Function implemented in subclass
        """

    def _get_contact_geometry_data(self, q):
        """
        Function calculates a vector - point of contact from global coordinates to local coordinates of each body
        edge body, node body
        normals
        tangents

        """
        contact_nodes = []
        for i, overlap_pair in enumerate(self.AABB_list_of_overlap_pairs):
            _contact_nodes = overlap_pair.get_contact_nodes()

            #   append contact nodes of each overlap pair to one list in contact object
            contact_nodes.append(_contact_nodes)

        # define empty list of normals and tangents for contact points
        self._n = []
        self._t = []
        # print "len(contact_nodes) =", len(contact_nodes)

        # for _contact_node in contact_nodes:
        #     print "_contact_node =", _contact_node
        # print "self.AABB_list_of_overlap_pairs =", self.AABB_list_of_overlap_pairs
        # time.sleep(100)
        #    get distance objects that have contacts from overlap pair objects
        min_overlap_pair_of_AABB = min(self.AABB_list_of_overlap_pairs, key=attrgetter('distance_min._distance'))

        #   get contact nodes
        self.min_distance_obj = min_overlap_pair_of_AABB.distance_min

        #   create edge and node body pointers
        self.__edge_body = min_overlap_pair_of_AABB.edge_body
        self.__node_body = min_overlap_pair_of_AABB.node_body

        #    normal and tangent of contact in GCS
        self._n_GCS = self.min_distance_obj.get_normal_2D()
        self._t_GCS = self.min_distance_obj.get_tangent_2D()

        #    list of new names for object attributes
        self.u_P_attribute_name_list = ["u_iP", "u_jP"]
        self.theta_ij_attribute_name_list = ["theta_i", "theta_j"]
        # self.theta_list = []

        #    normal, tangent of each body
        self._n_list = []
        self._t_list = []

        #    new contact properties
        self.u_P_list = []
        #    list of contact points in global coordinates
        self.r_P = self.min_distance_obj.node

        # logging.getLogger("DyS_logger").info("Contact point in GCS: %s", self.r_P)
        #    calculate a vector of contact point for each body in body LCS
        for body, body_id, u_P, theta_str in zip(self.body_list, self.body_id_list, self.u_P_attribute_name_list, self.theta_ij_attribute_name_list):
            #    assign pointer to a body
            # _body = bodies[body_id]
            # self.body_list.append(_body)

            #    get R, theta of body from q
            q_body = q2q_body(q, body_id)

            #    update coordinates and angle in 2D, R, theta
            body.update_coordinates_and_angles_2D(q_body)

            #    distance in local coordinate system of a body
            u_P_val = gcs2cm_lcs(self.r_P, q_body[0:2], q_body[2])

            #    add distance as object attribute
            setattr(self, u_P, u_P_val)
            setattr(self, theta_str, q_body[2])
            self.theta_list.append(q_body[2])

            #    append u_P value to list
            self.u_P_list.append(u_P_val)

        # body i is edge and body j is node
        if self.body_list[0] is self.min_distance_obj.edge_body and self.body_list[1] is self.min_distance_obj.node_body:
            n_i_GCS = self._n_GCS
            t_i_GCS = self._t_GCS

            n_j_GCS = -n_i_GCS
            t_j_GCS = -t_i_GCS

            self.__edge_body_id = self.body_id_list[0]
            self.__node_body_id = self.body_id_list[1]

        # body i is node and body j is edge
        if self.body_list[0] is self.min_distance_obj.node_body and self.body_list[1] is self.min_distance_obj.edge_body:
            n_i_GCS = -self._n_GCS
            t_i_GCS = -self._t_GCS

            n_j_GCS = -n_i_GCS
            t_j_GCS = -t_i_GCS

            self.__edge_body_id = self.body_id_list[1]
            self.__node_body_id = self.body_id_list[0]

        # list of normals and tangents of each body
        self._n_GCS_list = [n_i_GCS, n_j_GCS]
        self._t_GCS_list = [t_i_GCS, t_j_GCS]
        # print "self._n_list =", self._n_list
        # print "self._t_list =", self._t_list

        #   transform contact point data in GCS to LCS of each body (node and edge)
        self._contact_geometry_LCS(q)

        self._contact_geometry_GCS(q)

    def _contact_geometry_LCS(self, q):
        """
        Function calculates the contact geometry from GCS to LCS (only once) 
        This data is then used to calculate the contact data in GCS as a function of R, theta (i, j) at every time step
        of contact during numerical integration.
        :return:
        """
        self.u_P_LCS_list = []
        #    distance in local coordinate system of a body
        _node_GCS = self.min_distance_obj.get_node_2D()
        _edge_GCS = self.min_distance_obj.get_edge_2D()

        #   get body coordinates
        for body in self.body_list:
            #    get R, theta of body from q
            q_body = q2q_body(q, body.body_id)

            #   calculate contact node in GCS
            #   if body is node body - calculate node in GCS
            if body == self.__node_body:
                self.node_LCS = gcs2cm_lcs(_node_GCS, q_body[0:2], q_body[2])
                self.u_P_LCS_list.append(self.node_LCS)

            # if body is edge body - calculate edge nodes in GCS
            elif body == self.__edge_body:
                for i, edge_node_GCS in enumerate(_edge_GCS):
                    __edge_node_LCS = gcs2cm_lcs(edge_node_GCS, q_body[0:2], q_body[2])
                    self.edge_LCS[i] = __edge_node_LCS

                # calculate contact point on edge - normal projection of node to edge
                #   of opposite body
                _node_LCS = gcs2cm_lcs(_node_GCS, q_body[0:2], q_body[2])
                self.u_P_LCS_list.append(_node_LCS)

    def _contact_geometry_GCS(self, q):
        """

        :param q:
        :return:
        """
        for body in self.body_list:
            #    get R, theta of body from q
            q_body = q2q_body(q, body.body_id)

            if body == self.__node_body:
                self.node_GCS = cm_lcs2gcs(self.node_LCS, q_body[0:2], q_body[2])

            elif body == self.__edge_body:
                for i, _edge_node_LCS in enumerate(self.edge_LCS):
                    _edge_node_GCS_val = cm_lcs2gcs(_edge_node_LCS, q_body[0:2], q_body[2])

                    self.edge_GCS[i] = _edge_node_GCS_val

    def contact_velocity(self, q):
        """
        Function evaluates relative contact velocity vectors in normal and tangent direction
        :param q:
        :return:
        """
        dr_P = []
        for body_id, u_P in zip(self.body_id_list, self.u_P_LCS_list):
            #   body velocity, R, theta
            dR = q2dR_i(q, body_id)
            theta = q2theta_i(q, body_id)

            #    dtheta - omega
            dtheta = q2dtheta_i(q, body_id)

            #    point velocity
            dr_P_body = dr_contact_point_uP(dR, theta, dtheta, u_P)

            #    add to list
            dr_P.append(dr_P_body)

        #   relative contact velocity vector
        _dq = dr_P[1] - dr_P[0]
        #   relative contact velocity
        #   normal direction
        _dq_n = np.dot(_dq, self._n_GCS)
        #   tangent direction
        _dq_t = np.dot(_dq, self._t_GCS)

        return _dq_n, _dq_t

    def contact_point_object_properties(self):
        """

        :return:
        """
        for i, contact_point_obj in enumerate(self._contact_point_obj_list):
            print "------------------------------------------------------"
            print i
            pprint(vars(contact_point_obj))

    def write_to_file(self):
        """
        Save contact solution data to file.
        """
        # #   number of steps
        # n = len(self._step_num_solution_container)
        # #    get joint object data
        # data = np.vstack((np.array(self._step_num_solution_container),
        #                   np.array(self._status_container),
        #                   np.array(self._step_size_solution_container),
        #                   np.array(self._t_solution_container),
        #                   np.array(self._distance_solution_container),
        #                   np.array(self._dqn_solution_container),
        #                   np.array(self._dqt_solution_container),
        #                   np.array(self._Fn_solution_container),
        #                   np.array(self._Ft_solution_container),
        #                   np.array(self._u_P_solution_container).reshape(4, n),
        #                   np.array(self._r_P_solution_container).reshape(4, n)))

        #   create solution data object
        # solution_data = SolutionDataContact(name=self._name, MBD_item=self)

        #   add solution data to object
        # solution_data.add_data(data)

        #    write data to file
        self.solution_data.write_to_file()

        self._solution_filename = "solution_data_" + self._name + "_sol" + self._solution_filetype
        self._solution_filename = check_filename(self._solution_filename)

        # self._comments ="contact model: "+self.contact_model._type+"\n"
        #
        # #   save to selected file format
        # if self._solution_filetype == ".dat" or self._solution_filetype == ".sol":
        #     self._write_to_txt_file()
        # elif self._solution_filetype == ".xlsx" or self._solution_filetype == ".xls":
        #     self._write_to_excel_file()
        # elif self._solution_filetype == ".csv":
        #     self._write_to_csv_file()
        #
        # logging.getLogger("DyS_logger").info("Contact data (time and geometry) saved to file:\n%s" % self._solution_filename)

    def get_contact_force(self):
        """
        
        """
        for _force_n, _force_t in zip(self._Fn_list, self._Ft_list):
            _force_n.get_data()
            _force_t.get_data()

    def get_contact_point(self, step=None):
        """
        Method prints contact point in GCS
        :return:
        """
        if step is not None:
            _u_P_solution_container = np.array(self._u_P_solution_container)
            print "uP_GCS =", _u_P_solution_container[int(step), 4:6]

    def plot_Fn_delta(self, color=None, plot_label=None):
        """

        :return:
        """
        _contact_indices = self.__contact()

        label, color = self._plot_info(color)

        x = abs(self._distance_solution_container[_contact_indices])
        y = abs(self._Fn_solution_container[_contact_indices])

        plt.plot(x, y, marker=None, color=color, label=label)

    def plot_Fn_time(self, color=None):
        """

        :return:
        """
        _contact_indices = self.__contact()

        label, color = self._plot_info(color)

        plt.plot(abs(self._t_solution_container[_contact_indices]), abs(self._Fn_solution_container[_contact_indices]), marker=None, color=color, label=label)

    def _contact(self):
        """

        :return:
        """
        #   indices where contact is present, value is 1
        _contact_indices = np.nonzero(self._status_container)  # _Fn_solution_container, _status_container
        return _contact_indices

    def _plot_info(self, color=None):
        """

        :return:
        """
        #   label
        label = self.contact_model._type

        #   plot
        if color is None:
            color = self.color
        else:
            color = color
        return label, color

    def plot(self, q):
        """
        Plot contact geometry defined in subclass
        :param q:
        :return:
        """

    def print_solution_containers(self):
        """

        :return:
        """
        self.solution_data.print_solution_containers()

    def testing(self):
        print "self.contact_models =", self.contact_models
        for contact_model in self.contact_models:
            print contact_model
            pprint(vars(contact_model))


if __name__ == "__main__":
    #   solution_data_contact_12_sol_12.xlsx
    #   solution_data_contact.xlsx
    #   solution_data_contact_12_sol.xlsx
    filename = "solution_data_contact_12_sol.xlsx"

    contact = Contact()

    fig = plt.figure(num=1, figsize=(8, 6), dpi=100, facecolor='w', edgecolor='k')
    ax = plt.subplot(111)
    ax.ticklabel_format(style='sci', scilimits=(-4, 4), axis='both')

    contact.read_file(filename)
    # contact.plot_Fn_delta()

    contact.plot_Fn_time()

    plt.grid()
    plt.show()

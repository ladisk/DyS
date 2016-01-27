"""
Created on 21. feb. 2014

@author: lskrinjar
"""
import os
import time
import logging
from operator import attrgetter
import xlsxwriter
import xlrd
from matplotlib import pyplot as plt
from OpenGL.GL import *
import copy


from MBD_system.MBD_system_items import ContactItem
from MBD_system.force.force import Force
from MBD_system.contact.bounding_box.bounding_box import AABB
from MBD_system.contact.overlap_pair.overlap_pair import OverlapPair
from evaluate_distance import evaluate_distance_2D
from contact_C_matrix import *
from MBD_system.q2dR_i import q2dR_i
from MBD_system.q2theta_i import q2theta_i
from MBD_system.q2dtheta_i import q2dtheta_i
from MBD_system.q2q_body import q2q_body
from MBD_system.transform_cs import cm_lcs2gcs, gcs2cm_lcs, gcs2lcs_z_axis
from MBD_system.dr_contact_point_uP import dr_contact_point_uP
from MBD_system.contact_model.contact_model import ContactModel
from MBD_system.friction_model.friction_model import FrictionModel
from MBD_system.contact.contact_geometry.contact_geometry import ContactGeometry
from MBD_system.check_filename import check_filename
from MBD_system.extract_from_dictionary_by_string_in_key import extract_from_dictionary_by_string_in_key


class Contact(ContactItem):
    """
    classdocs
    """
    __id = itertools.count(0)

    def __init__(self, _type, body_id_i, body_id_j, name=None, contact_model_type=None, friction_model_type=None, properties_dict=[], parent=None):
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

        #    this has to be after attributes contact_id and _name as name is constructed from contact_id and _name
        super(Contact, self).__init__(self._name, parent)

        #   parent
        self._parent = parent

        #    contact type
        self._type = _type

        #   supported types
        self._types = ["general",
                        "revolute clearance joint",
                        "prismatic clearance joint",
                        "contact sphere-sphere",
                        "contact plane-sphere",
                        "pin-slot clearance joint linear",
                        "pin-slot clearance joint radial"]

        #   contact geometry (profile) - list of object pairs
        self.contact_geometry_list = []

        #   position of contact in z direction
        self.z_dim = 0


        #   solution atrributes
        #   save solution - properties
        self.save_to_file = False
        self._solution_filename = None
        #   a solutions list to store solution data of contact object of each simulation
        self.solutions = []
        #   solution options
        self._save_options = "Discard"
        #   solution file type, optional: .dat, .xlsx, .csv
        self._solution_filetype = ".xlsx"


        #    numerical error when calculating distance (point to line), in m
        self.distance_TOL = 1.E-7


        #   visualization properties - open gl properties
        self._color_GL = np.array([0.4, 0.2, 0.2], dtype="float32")
        self._scale_GL = 1E-2

        #    contact body ids
        self.body_id_i = body_id_i
        self.body_id_j = body_id_j
        self.body_id_list = [self.body_id_i, self.body_id_j]
        self.body_list = []
        self.contact_geometry_list = []

        #   contact model
        #   default properties
        self.coef_of_restitution = 1
        self.contact_model_type = contact_model_type
        if self.contact_model_type is None:
            self.contact_model_type = "hertz"
        #   extract only properties that relate to contact model object properties
        self.properties_contact_model = extract_from_dictionary_by_string_in_key(properties_dict, "contact_model.")
        self.contact_model = ContactModel(self.contact_model_type, c_r=self.coef_of_restitution, properties_dict=properties_dict, parent=self)

        #    friction model
        self.friction_model_type = friction_model_type
        if self.friction_model_type is None:
            self.friction_model_type = "ideal"
        self._dq_t_TOL = 1E-3
        self.coef_of_friction_dynamic = 0
        self.coef_of_friction_static = 0
        self.properties_friction_model = extract_from_dictionary_by_string_in_key(properties_dict, "contact_model.")
        self.friction_model = FrictionModel(self.friction_model_type, coef_of_friction_dynamic=self.coef_of_friction_dynamic, coef_of_friction_static=self.coef_of_friction_static, parent=self)

        #    set properties dictionary as object property
        self.properties = properties_dict
        if self.properties is not []:
            self.add_additional_parameters(self.properties)

        #    status of contact
        #    0 - no contact
        #    -1 - contact is detected, two bodies have already collided and are in collision, but the
        #    current initial penetration depth is too large (halve the time step)
        #    1 - contact is detected
        self.status = 0

        #    initialize a list of AABB
        self.AABB_i = None
        self.AABB_j = None
        self.AABB_list = []

        #    initialize a list of contact pairs
        self.AABB_list_of_overlap_pairs = []

        #   contact node-edge points to calculate contact properties
        self.node_LCS = None
        self.node_GCS = None
        self.edge_LCS = [None, None]
        self.edge_GCS = [None, None]

        #   contact point(s) in GCS
        self.u_P_GCS = None
        self.u_iP_GCS = None
        self.u_jP_GCS = None
        self.u_P_list_GCS = [self.u_iP_GCS, self.u_jP_GCS]

        #   contact point in body LCS
        self.u_P_LCS = None
        self.u_iP_LCS = None
        self.u_jP_LCS = None
        self.u_P_LCS_list = [self.u_iP_LCS, self.u_jP_LCS]

        #   list of normals and tangents in LCS of each body in contact
        #   tangents
        self._t_GCS = None
        self._t_i_LCS = None
        self._t_j_LCS = None
        self._t_LCS_list = [self._t_i_LCS, self._t_j_LCS]
        #   normals
        self._n_GCS = None
        self._n_i_LCS = None
        self._n_j_LCS = None
        self._n_LCS_list = [self._n_i_LCS, self._n_j_LCS]
        self._n_list = []

        #    contact distance - penetration depth
        self._delta = 0
        self.contact_detected = False  # True or False
        self.contact_distance_inside_tolerance = False
        self._contact_point_found = False

        #   list of contact forces
        self.contact_bodies_added_to_list = False
        self.list_of_contact_force_objects_constructed = False
        self.Fn = 0
        self.Ft = 0
        self._Fn_list = []
        self._Ft_list = []

        #   define relative contact velocity components as object attributes
        self._dq_n = 0
        self._dq_t = 0
        self.initial_contact_velocity_calculated = False

        if self._parent is not None:
            for body_id in self.body_id_list:
                _Fn = Force(body_id, force_name=self._name + "_Fn_on_body_" + str(body_id))
                # _Fn._visible = True
                _Ft = Force(body_id, force_name=self._name + "_Ft_on_body_" + str(body_id))
                _Ft._scale_GL = self._scale_GL
                #   add pair of contact forces to forces list of MBD system
                self._parent._parent.forces.append(_Fn)
                self._parent._parent.forces.append(_Ft)

                #   add pair of contact forces to forces list of contact
                self._Fn_list.append(_Fn)
                self._Ft_list.append(_Ft)

                self._parent._parent.bodies[body_id].forces.append(_Fn)
                self._parent._parent.bodies[body_id].forces.append(_Ft)

            self.list_of_contact_force_objects_constructed = True
        
        # self._Fn_list[0]._visible = False
        # self._Ft_list[0]._visible = False

        #   create solution containers
        self._solution_containers()

        #   list of markers
        self.markers = self._create_markers()

    def _create_markers(self):
        """
        Function creates markers
        :return:
        """
        markers = []
        return markers

    def _create_Fn_forces(self):
        """
        Function creates list of contact forces as objects
        :return:
        """

    def _create_Ft_forces(self):
        """
        Function creates list of contact forces as objects
        :return:
        """

    def add_additional_parameters(self, dict):
        """
        Add additional attributes from dictionary self.properties
        """
        #    set object properties from dictionary
        for key in dict:
            #    set Contact object additional properties
            if "body" not in key and "contact " not in key:
                if key == "name":
                    val = dict[key]
                    key = "_name"
                    dict[key] = val

                setattr(self, key, dict[key])

                dict.pop(key, None)
            #    remove key from dictionary after added to object attributes
            elif "." not in key:
                dict.pop(key, None)
            else:
                None

    def _solution_containers(self):
        """

        :return:
        """
        #    solution containers over integration times
        self._step_num_solution_container = np.empty(1)
        self._step_num_solution_container.fill(0)

        self._status_container = np.empty(0)
        self._status_container.fill(0)

        self._step_size_solution_container = np.empty(1)
        self._step_size_solution_container.fill(0)

        self._t_solution_container = np.empty(1)
        self._t_solution_container.fill(0)
        
        self._distance_solution_container = np.empty([1])
        self._distance_solution_container.fill(0)

        self._dqn_solution_container = np.empty([1])
        self._dqn_solution_container.fill(0)

        self._dqt_solution_container = np.empty([1])
        self._dqt_solution_container.fill(0)

        self._Fn_solution_container = np.empty([1])
        self._Fn_solution_container.fill(0)

        self._Ft_solution_container = np.empty([1])
        self._Ft_solution_container.fill(0)

        self._u_P_solution_container = []

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

        self.__solution_container_dict = {"delta":self._distance_solution_container,
                                          "Fn":self._Fn_solution_container,
                                          "Ft":self._Ft_solution_container,
                                          "dqn":self._dqn_solution_container,
                                          "dqt":self._dqt_solution_container}

    def _reset_to_initial_state(self):
        """
        Reset main object attributes to initial value
        """
        self.initial_contact_velocity_calculated = False
        self.status = 0
        # self._solution_containers()

        #    reset force object solution data containers
        for force in self._Fn_list:
            force._reset_to_initial_state()
    
    def __set_ACCF_parameters(self):
        """
        Set additional object attributes based if contact type is ACCF
        """
        self._phase = "C"  #    C - compression, E-expansion
    
    def __set_ECF_N_parameters(self):
        """
        Set additional object attributes based if contact type is ECF-N
        """
        self.K = 10
        self.c_r = 0.1
        
    def create_bounding_box_for_each_body_in_contact(self, bodies=[]):
        """
        Create bounding box object for each body
        """
        self.__body_list = ["body_i", "body_j"]

        for body_id, __body in zip(self.body_id_list, self.__body_list):
            #    assign body object from body list to new variable name
            body_ = bodies[body_id]
            setattr(body_, "skin_thickness", 10 * self.distance_TOL)

            #   remove duplicate nodes
            # _nodes, _normals = body_.geom.geom_data.remove_duplicate_2D()

            #   build contact geometry object and append it to lists (in contact and in body)
            _contact_area = self.properties[__body+".contact_area"]
            _contact_geometry_type = self.properties[__body+".contact_geometry_type"]

            if _contact_geometry_type == "area":
                pass
            elif _contact_geometry_type == "edge":
                pass
            elif _contact_geometry_type == "body":
                pass

            _contact_geometry = ContactGeometry(_contact_geometry_type, body_, area=_contact_area, parent=self)

            self.contact_geometry_list.append(_contact_geometry)
            body_.contact_geometries.append(_contact_geometry)

            #    get nodes and normals from contact geometry object
            _nodes, _normals = _contact_geometry.get_data()

            #    if value is assigned to position_of_plane_in_z_direction only 2D nodes are put in argument
            if self.z_dim is not None:
                _nodes = _nodes#[:, 0:2]
                _normals = _normals#[:, 0:2]

            #    add additional parameter to body object - body point that has maximum distance
            #    from center of body mass
            body_.uP_i_max = np.linalg.norm(np.array(np.amax(_nodes, axis=0)), ord=2)
        
            #    head AABBTree object is constructed
            body_.AABBtree = AABB(nodes_LCS=_nodes, normals_LCS=_normals, parent=self, parent_body=body_)
            # print "body_.AABBtree =", body_.AABBtree
            # pprint(vars(body_.AABBtree))
            body_.AABBtree._visible = True
        
            #    is additional properties are defined, than they are assigned to AABB objects
            for key in self.properties:
                if body_id == self.body_id_i and "body_i" in key:
                    setattr(body_.AABBtree, key[key.rindex('.') + 1:], self.properties[key])
        
                if body_id == self.body_id_j and "body_j" in key:
                    setattr(body_.AABBtree, key[key.rindex('.') + 1:], self.properties[key])

            #    create AABBTree object with recursion
            body_.AABBtree.construct()
        
            #    append AABB object to list
            self.AABB_list.append(body_.AABBtree)
        
        #    after the list of AABBtree for the contact is created (filled with object-data) each AABBtree
        #    is saved as different attribute for later calculations
        # print "self.AABB_list =", self.AABB_list
        if self.body_id_list is not []:
            self.AABB_i, self.AABB_j = self.AABB_list

        self.add_additional_parameters(self.properties)

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
            # #    step size
            # self._h = t - self._t_solution_container[-1]
            #
            # if self._step==self._step_num_solution_container[0]:
            #     self._step_num_solution_container[0] = self._step
            # # else:
            # #     pass
            #     # self._step_num_solution_container = np.append(self._step_num_solution_container, self._step)
            #
            # self._t_solution_container = np.append(self._t_solution_container, self.t)
            #
            # self._step_size_solution_container = np.append(self._step_size_solution_container, self._h)
            self._step_solution_accepted = True

        #    replace last value if new value is smaller than last value
        else:
            # #    step size
            # self._h = t - self._t_solution_container[-2]
            #
            # self._step_num_solution_container[-1] = self._step
            # self._t_solution_container[-1] = self.t
            #
            # self._step_size_solution_container[-1] = self._h
            self._step_solution_accepted = False

    def _track_data(self, step, h, t, q):
        """
        Function saves current calculated data at
        :return:
        """
        #   step number
        self._step_num_solution_container = np.append(self._step_num_solution_container, step)

        #   status
        self._status_container = np.append(self._status_container, self.status)

        #   step size
        self._step_size_solution_container = np.append(self._step_size_solution_container, h)

        #   time
        self._t_solution_container = np.append(self._t_solution_container, t)

        #   delta
        self._distance_solution_container = np.append(self._distance_solution_container, self._delta)

        #   dq
        self._dqn_solution_container = np.append(self._dqn_solution_container, self._dq_n)
        self._dqt_solution_container = np.append(self._dqt_solution_container, self._dq_t)

        #   F
        self._Fn_solution_container = np.append(self._Fn_solution_container, self.Fn)
        self._Ft_solution_container = np.append(self._Ft_solution_container, self.Ft)

        #   uPi, uPj
        self._u_P_solution_container.append(self.u_P_list_LCS.flatten())

    def __AABB_AABB_overlap(self, AABB_i, AABB_j):
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
        #   3D overlap check
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
            AABB.update_frame_geometry(q_i)

        #    check for overlap of each AABB frame
        if self.__AABB_AABB_overlap(AABB_i, AABB_j):
            #    check between children of AABB_i and AABB_j
            if AABB_i.has_children() and AABB_j.has_children():
                for __AABB_i in AABB_i.children:
                    for __AABB_j in AABB_j.children:  # self.AABB_j.children
                        self.check_for_overlap(q, __AABB_i, __AABB_j)
            
            #    check between AABB_j and children of AABB_i
            elif AABB_i.has_children() and not AABB_j.has_children():
                for __AABB_i in AABB_i.children:
                    self.check_for_overlap(q, __AABB_i, AABB_j)
            
            #    check between AABB_i and children of AABB_j
            elif not AABB_i.has_children() and AABB_j.has_children():
                for __AABB_j in AABB_j.children:
                    self.check_for_overlap(q, AABB_i, __AABB_j)

            #    if none of the two AABB has children, the two AABB are used to create contact pair object
            else:
                __AABB_overlap_pair = OverlapPair(_AABB_i=AABB_i, _AABB_j=AABB_j, parent=self)
                #    add created object to list of all AABB overlap pairs
                self.AABB_list_of_overlap_pairs.append(__AABB_overlap_pair)
            
    def check_for_contact(self, q):
        """
        Function check for contact between contact pair of AABBs
        returns:
                -1 - first contact indentation is too large, return to previous time step and continue with smaller step size
                0 - no contact
                +1 - contact detected
        """
        # print "check_for_contact() - contact.py"
        #    predefine and reset to empty lists
        _distance_list_values = []
        _inside_list = []
        _distance_list = []

        #    loop through list of overlapped AABB object pairs
        print "number of overlap contact pairs=", len(self.AABB_list_of_overlap_pairs)
        for overlap_pair_obj in self.AABB_list_of_overlap_pairs:
            #    check for contact of overlapped AABB pair (object)
            overlap_pair_obj.check_for_contact_2D(q)
            
            #    get distance list of every AABB overlap pair object of contact and add it to contact attribute
            __distance_list_of_overlap_pair = overlap_pair_obj.distance_list
            
            #    extend distance list from AABB overlap object with distance list of every contact pair object list
            _distance_list.extend(__distance_list_of_overlap_pair)
        # print "_distance_list =", _distance_list
#         pprint(vars(_distance_list[1]))
        # if 0.019 < self.t < 0.022:
        # print "t =", self.t
        # pprint(vars(_distance_list[1]))
        # print "y =", _distance_list[1].edge_node_1[1]
        # print "max_penetration_depth =",  _distance_list[1].max_penetration_depth
        # print "inside =", _distance_list[1]._inside
        # print "normal =", _distance_list[1].normal
        # print "self._distance_vector =", _distance_list[1]._distance_vector
        # print "self._distance =", _distance_list[1]._distance
        # print "self.normal_plus_distance =", _distance_list[1].normal_plus_distance
            # print "_TF1 =", _distance_list[1]._TF1
            # print "_TF2(never TRUE) =", _distance_list[1]._TF2
            # print "_TF3 =", _distance_list[1]._TF3
            # print "_TF4 =", _distance_list[1]._TF4
        # print "------------------------------"
#         for i, _dist in enumerate(_distance_list):
#             print "_dist._distance =", _dist._distance, "_inside =", _dist._inside

        #    distance list values
        _distance_list_values = np.array([_dist._distance for _dist in _distance_list])

        #     print "_dist._inside =", _dist._inside
        #     print "edge body =", _dist.edge_body._name
        #     print "node body =", _dist.node_body._name
        #     print i
        #     pprint(vars(_dist))
#         print "d =", _distance_list[1]._distance
        # time.sleep(100)
        #    predefine min distance attribute
        self._distance = None
        self._distance_sign = None

        _inside_list = [_dist._inside for _dist in _distance_list]
        # print "_inside_list =", _inside_list
        # for dist in _inside_list:
        #     print "------------------------------"
        #     dist
        
#         print "_distance_list_values =", _distance_list_values
        #    position of minimal distance (value)
        min_distance_val_pos = _distance_list_values.argmin(axis=0)

        self._distance_sign_list_values = np.array([_dist._distance_sign for _dist in _distance_list])
        
        self._distance = min(_distance_list_values)
        
        self._distance_sign = min(self._distance_sign_list_values)

        #   add calculated data to container
        # if self._step_solution_accepted:
        #     self._distance_solution_container = np.append(self._distance_solution_container, self._distance_sign)  # _distance_sign, _distance
        # else:
        #     self._distance_solution_container[-1] = self._distance_sign

        #    all calculated distances are greater than minimum tolerance for contact
        #    returns value 1 - to divide time step
        # print "min(_distance_list_values) =", min(_distance_list_values)
        if any(_inside_list) and (_distance_list_values > self.distance_TOL).all():
            self.contact_detected = True
            self.status = -1
            return self.status
        
        #    contact is detected
        if (_distance_list_values <= self.distance_TOL).any():
            self.contact_detected = True
            self.contact_distance_inside_tolerance = True
            self.status = 1
            # self._status_container = np.append(self._status_container, self.status)
            return self.status

        #    all calculated distances are greater than tolerance and bodies are not in contact
        self.status = 0
        # self._status_container = np.append(self._status_container, self.status)
        self.no_contact()

        # print "self.contact_detected =", self.contact_detected
        # print "self._distance =", self._distance
        return self.status

    def _get_contact_geometry_data(self, t, q, bodies):
        """
        Function calculates a vector - point of contact from global coordinates to local coordinates of each body
        edge body, node body
        normals
        tangents

        """
        contact_nodes = []
        for i, overlap_pair in enumerate(self.AABB_list_of_overlap_pairs):
            print "i =", i
            _contact_nodes = overlap_pair.get_contact_nodes()

            #   append contact nodes of each overlap pair to one list in contact object
            contact_nodes.append(_contact_nodes)

        #   define empty list of normals and tangents for contact points
        self._n = []
        self._t = []
        print "len(contact_nodes) =", len(contact_nodes)

        for _contact_node in contact_nodes:
            print "_contact_node =", _contact_node
        print "self.AABB_list_of_overlap_pairs =", self.AABB_list_of_overlap_pairs
        time.sleep(100)
        #    get distance objects that have contacts from overlap pair objects
        min_overlap_pair_of_AABB = min(self.AABB_list_of_overlap_pairs, key=attrgetter('distance_min._distance'))

        #   get contact nodes
        self.min_distance_obj = min_overlap_pair_of_AABB.distance_min

        #   create edge and node body pointers
        self.__edge_body = min_overlap_pair_of_AABB.edge_body
        self.__node_body = min_overlap_pair_of_AABB.node_body

        #    normal and tanget of contact
        self._n = self.min_distance_obj.get_normal_2D()
        self._t = self.min_distance_obj.get_tangent_2D()

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
        for body_id, u_P, theta_str in zip(self.body_id_list, self.u_P_attribute_name_list, self.theta_ij_attribute_name_list):
            #    assign pointer to a body
            _body = bodies[body_id]
            self.body_list.append(_body)

            #    get R, theta of body from q
            q_body = q2q_body(q, body_id)

            #    update coordinates and angle in 2D, R, theta
            _body.update_coordinates_and_angles_2D(q_body)

            #    distance in local coordinate system of a body
            u_P_val = gcs2cm_lcs(self.r_P, q_body[0:2], q_body[2])

            #    add distance as object attribute
            setattr(self, u_P, u_P_val)
            setattr(self, theta_str, q_body[2])
            self.theta_list.append(q_body[2])

            #    append u_P value to list
            self.u_P_list.append(u_P_val)

        #   body i is edge and body j is node
        if self.body_list[0] is self.min_distance_obj.edge_body and self.body_list[1] is self.min_distance_obj.node_body:
            n_i = self._n
            t_i = self._t

            n_j = -n_i
            t_j = -t_i

            self.__edge_body_id = self.body_id_list[0]
            self.__node_body_id = self.body_id_list[1]

        #   body i is node and body j is edge
        if self.body_list[0] is self.min_distance_obj.node_body and self.body_list[1] is self.min_distance_obj.edge_body:
            n_i = -self._n
            t_i = -self._t

            n_j = -n_i
            t_j = -t_i

            self.__edge_body_id = self.body_id_list[1]
            self.__node_body_id = self.body_id_list[0]

        #   list of normals and tangents of each body
        self._n_list = [n_i, n_j]
        self._t_list = [t_i, t_j]
        print "self._n_list =", self._n_list
        print "self._t_list =", self._t_list

        #   transform contact point data in GCS to LCS of each body (node and edge)
        self._contact_geometry_LCS(q)

        self._contact_geometry_GCS(q)

    def _contact_geometry_LCS(self, q):
        """
        Function calculates (only once) the contact geometry from GCS to LCS
        this data is then used to calculate the contact data in GCS as a function of R, theta (i, j) at every time step
        during contact
        :return:
        """
        self.u_P_list_LCS = []
        #    distance in local coordinate system of a body
        _node_GCS = self.min_distance_obj.get_node_2D()
        # u_P_GCS
        print "_node_GCS =", _node_GCS
        # glVertex3f(_node_GCS[0],_node_GCS[1],self.z_dim)

        _edge_GCS = self.min_distance_obj.get_edge_2D()

        #   get body coordinates
        for body in self.body_list:
            #    get R, theta of body from q
            q_body = q2q_body(q, body.body_id)

            #   calculate contact node in GCS
            #   if body is node body - calculate node in GCS
            if body == self.__node_body:
                self.node_LCS = gcs2cm_lcs(_node_GCS, q_body[0:2], q_body[2])
                self.u_P_list_LCS.append(self.node_LCS)

            #   if body is edge body - calculate edge nodes in GCS
            elif body == self.__edge_body:
                for i, edge_node_GCS in enumerate(_edge_GCS):
                    __edge_node_LCS = gcs2cm_lcs(edge_node_GCS, q_body[0:2], q_body[2])
                    self.edge_LCS[i] = __edge_node_LCS

                #   calculate contact point on edge - normal projection of node to edge
                #   of opposite body
                _node_LCS = gcs2cm_lcs(_node_GCS, q_body[0:2], q_body[2])
                self.u_P_list_LCS.append(_node_LCS)

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

    def contact_update(self, step, t, q):
        """
        Function updates status only for (general) contact
        :return:
            status - status of contact 0-no contact, 2-contact
        """
        # print "contact_update()"
        self._step = step

        #   calculate relative contact velocity in normal direction at time t
        self._dq_n, self._dq_t = self._contact_velocity(q)

        #   calculate contact geometry in GCS to calculate contact distance
        self._contact_geometry_GCS(q)

        #   calculate distance: node-edge
        _distance, _inside = evaluate_distance_2D(self.node_GCS, self.edge_GCS[0], self._n, self._t)
        
        #   contact is present
        if _inside and self._distance < self.distance_TOL:
            self._distance = -_distance
            self.status = 1

        #   no contact
        else:
            self._distance = _distance
            self.status = 0
            self.contact_detected = False 
            self.contact_distance_inside_tolerance = False
            self._contact_point_found = False
            self.initial_contact_velocity_calculated = False

            for _Fn in self._Fn_list:
                _Fn._update_force_vector(np.zeros(2))

        # print "self.status - contact_update() =", self.status
        return self.status

    def solve(self, t, q):
        """
        Calculate contact parameters
        returns:
        """
        # print "self._contact_point_found =", self._contact_point_found
        #    calculate coordinates of contact point from global coordinates in local coordinates of each body in contact
        if not self._contact_point_found:
            self._get_contact_geometry_data(q)
            self._contact_point_found = True
        else:
            pass
            # self._distance, self._delta = self._contact_geometry_GCS(q)
            # print "contact update"

        # print "self._contact_point_found =", self._contact_point_found

        #   kinematic properties of contact point
        #   initial contact velocity
        if not self.initial_contact_velocity_calculated:
            self._dq0_n, self._dq0_t = self._contact_velocity(q)
            self.contact_model.set_dq0(self._dq0_n, self._dq0_t)
            self.initial_contact_velocity_calculated = True

        # if self.contact_detected:
        self._distance, self._delta = self._contact_geometry_GCS(q)
        # print "self._dq0_n, self._dq0_t =", self._dq0_n, self._dq0_t

        #   current contact velocity at time t
        self._dq_n, self._dq_t = self._contact_velocity(q)
        # print "self._dq_n, self._dq_t =", self._dq_n, self._dq_t
        # time.sleep(100)
        #   current contact velocity at time t
        # self._dq_n, self._dq_t = self._contact_velocity(q)

        if self._type.lower() in self.__types:#== "general" or self._type.lower() == "revolute clearance joint" or self._type.lower() == "contact sphere-sphere":#ECF-N
            self._solve_ECF_N(t, q, self._delta, self._dq_n, self._dq_t)#self._delta
        else:
            raise ValueError, "Contact type not correct!"

    def _solve_ECF_N(self, t, q, _delta, _dq_n, _dq_t):
        """
        Args:
            _delta - (scalar value) contact depth in normal direction
            _dq_n - contact velocity in normal direction
            _dq_t - contact velocity in tangent direction
        returns:
            F - force vector (Fx, Fy)
        """
        # print "_solve_ECF_N()"

        # print "delta =", _delta
        # print "dqn =", _dq_n
        # print "q =", q
        #    normal and tangent force are calculated with selected contact and friction model
        self.Fn = self.contact_model.contact_force(_delta, _dq_n)# _dq_t, self._n, self._t
        # print "--------------------------------"
        # print "q =", q
        # print "t =%20.15f"%t, "Fn =", self.Fn, "delta =", _delta, "dqn =", _dq_n
        self.Ft = self.friction_model.friction_force(self.Fn, _dq_t)
        # print "Ft =", Ft

        #    contact force in GCS
        #    normal component
        # F_n = Fn * self._n  # add normal to get orientation of force in GCS
        # #    tangent component
        # F_t = Ft * self._t


        #   contact force (normal and tangent values) in GCS of body i and body j
        # F_i = F_n + F_t
        # F_j = F_n - F_t
        #   create list of forces for each body in contact
        # F_list = [F_i, F_j]

        #    check if contact is finished
        # print "TOL =", self.distance_TOL
        # print "_delta =", _delta
        # print "_delta0 =", self._delta0
        # print "abs(_delta) < self.distance_TOL =", abs(_delta) < self.distance_TOL
        #   contact
        if _delta <= self._delta0:#<self._delta0 and self._dq_n <= abs(self._dq0_n):#-self.distance_TOL:# and self._dq_n >=  self._dq0_n:#self.distance_TOL:# and (self._dq_n < 0):, _delta<-1*self.distance_TOL
            # print "exe???"
            # print "self._Fn_list =", self._Fn_list
            # print "self._Ft_list =", self._Ft_list
            # print "self.u_P_list_LCS =", self.u_P_list_LCS
            # print "self._n_list =", self._n_list
            # print "self._t_list =", self._t_list

            for _force_n, _force_t, _u_P, _n, _t in zip(self._Fn_list, self._Ft_list, self.u_P_list_LCS, self._n_list, self._t_list):
                # print "-------------------------"
                # print "body id =", _force_n.body_id
                # print "n =", _n
                # print "self.status =", self.status
                #   contact force values during contact
                if self.status == 1:
                    # print "self.status =", self.status
                    _Fn = self.Fn * _n
                    _Ft = self.Ft * _t
                    # print "_Fn =", _Fn
                    # print "_Ft =", _Ft
                    # print "self.status =", self.status
                else:
                    _Fn = self.Fn = np.zeros(2)
                    _Ft = self.Ft = np.zeros(2)
                    self._contact_point_found = False
                # print "_u_P =", _u_P
                # print "_Fn =", _Fn
                _force_n.update(self._step, F=_Fn, u_P=_u_P)
                _force_t.update(self._step, F=_Ft, u_P=_u_P)
                # pprint(vars(_force_n))
            # time.sleep(100)
            #   store data to solution containers
            # self._distance_solution_container = np.append(self._distance_solution_container, _delta)

            # self._dqn_solution_container = np.append(self._dqn_solution_container, _dq_n)
            # self._dqt_solution_container = np.append(self._dqt_solution_container, _dq_t)

        #   no contact
        else:
            self._contact_point_found = False
            self.initial_contact_velocity_calculated = False
            self.contact_distance_inside_tolerance = False
            self.contact_detected = False
            self.status == 0
            self.no_contact()
            # time.sleep(10)

        # self._status_container = np.append(self._status_container, self.status)

    def _contact_velocity(self, q):
        """
        Function evaluates relative contact velocity vectors in normal and tangent direction
        """
        dr_P = []
        for body_id, u_P in zip(self.body_id_list, self.u_P_list_LCS):
            #   body velocity, R, theta
            dR = q2dR_i(q, body_id)
            theta = q2theta_i(q, body_id)

            #    dtheta - omega
            dtheta = q2dtheta_i(q, body_id)

            #    point velocity
            dr_P_body = dr_contact_point_uP(dR, theta, dtheta, u_P)

            #    add to list
            dr_P.append(dr_P_body)

        #    relative contact velocity vector
        _dq = dr_P[1] - dr_P[0]
        #   relative contact velocity
        #   normal direction
        _dq_n = np.dot(_dq, self._n)
        #   tangent direction
        _dq_t = np.dot(_dq, self._t)

        return _dq_n, _dq_t

    def no_contact(self):
        """
        Method sets values to contact attributes when there is no contact
        """
        self.Fn = 0
        self.Ft = 0
        # self._Fn_solution_container = np.append(self._Fn_solution_container, 0)
        # self._Ft_solution_container = np.append(self._Ft_solution_container, 0)

        self._dq_n = 0
        self._dq_t = 0
        # self._dqn_solution_container = np.append(self._dqn_solution_container, 0)
        # self._dqt_solution_container = np.append(self._dqt_solution_container, 0)

        #   update contact force values with zeros as no contact is present
        for _force_n, _force_t in zip(self._Fn_list, self._Ft_list):
            _force_n.update(self._step)
            _force_t.update(self._step)

    def no_overlap(self):
        """

        """
        if len(self._t_solution_container) == 1:
            pass
        else:
            self._distance_solution_container = np.append(self._distance_solution_container, np.nan)

        self.no_contact()
        # self._status_container = np.append(self._status_container, self.status)

    def mechanical_energy(self):
        """

        :return:
        """
        if self._delta < 0:
            _energy = 0.5 * self.contact_model.K * self._delta**2
            return _energy
        else:
            return 0

    def write_to_file(self):
        """
        Save contact solution data to file.
        """

        self._solution_filename = "solution_data_" + self._name + "_sol"+self._solution_filetype
        self._solution_filename = check_filename(self._solution_filename)

        self._comments ="contact model: "+self.contact_model._type+"\n"

        #   save to selected file format
        if self._solution_filetype == ".dat" or self._solution_filetype == ".sol":
            self._write_to_txt_file()
        elif self._solution_filetype == ".xlsx" or self._solution_filetype == ".xls":
            self._write_to_excel_file()
        elif self._solution_filetype == ".csv":
            self._write_to_csv_file()

        logging.getLogger("DyS_logger").info("Contact data (time and geometry) saved to file:\n%s" % self._solution_filename)

    def _write_to_txt_file(self):
        """

        :return:
        """
        #   column headers
        _header = "i-th step\tstatus\t\tdt\t\t\t\t\t\ttime \t\t\t\t\t delta \t\t\t\t\t dq_n \t\t\t\t\t Fn \t\t\t\t\t Ft"

        #   all data to one matrix
        solution_data = np.hstack((np.array([self._step_num_solution_container]).T, np.array([self._status_container]).T, np.array([self._step_size_solution_container]).T, np.array([self._t_solution_container]).T, np.array([self._distance_solution_container]).T, np.array([self._dqn_solution_container]).T, np.array([self._Fn_solution_container]).T, np.array([self._Ft_solution_container]).T))
        #   save to file
        np.savetxt(self._solution_filename, solution_data, fmt=['%5i', '%5i', '%20.10f', '%20.10f', '%20.10f', '%20.10f', '%20.10f', '%20.10f'], delimiter='\t', header=_header, comments = self._comments)

    def _write_to_excel_file(self):
        """

        :return:
        """
        self.__excel_columns()

        #   column headers
        _header = ["i-th step", "status", "dt", "time", "delta", "dq_n", "dq_t", "Fn", "Ft", "uPi_x", "uPi_y", "uPj_x", "uPj_y"]

        #   create an new Excel file and add a worksheet
        workbook = xlsxwriter.Workbook(self._solution_filename, {'nan_inf_to_errors': True})
        format_1 = workbook.add_format({'num_format': '#0.000000000000000'})
        format_2 = workbook.add_format({'num_format': '#0.00000000'})
        worksheet = workbook.add_worksheet(self.__workbook_name)

        #   write comments
        worksheet.write(0, 0, self._comments)
        #   write header
        worksheet.write_row(1, 0, _header)

        #   write solution data to columns
        #   step number
        worksheet.write_column(self.__column_start_write, self.__col_step_num_solution_container, self._step_num_solution_container)
        #   contact status
        worksheet.write_column(self.__column_start_write, self.__col_status_container, self._status_container)
        #   step size
        worksheet.write_column(self.__column_start_write, 2, self._step_size_solution_container, format_1)
        worksheet.set_column('C:C', 20)
        #   time array/column
        worksheet.write_column(self.__column_start_write, self.__col_t_solution_container, self._t_solution_container, format_1)
        worksheet.set_column('D:D', 20)
        #   distance array/column
        worksheet.write_column(self.__column_start_write, self.__col_distance_solution_container, self._distance_solution_container, format_1)
        worksheet.set_column('E:E', 20)

        #   dqn
        worksheet.write_column(self.__column_start_write, self.__col_dqn_solution_container, self._dqn_solution_container, format_1)
        worksheet.set_column('F:F', 22)

        #   dqt
        worksheet.write_column(self.__column_start_write, self.__col_dqt_solution_container, self._dqt_solution_container, format_1)
        worksheet.set_column('G:G', 22)

        #   Fn
        worksheet.write_column(self.__column_start_write, self.__col_Fn_solution_container, self._Fn_solution_container, format_2)
        worksheet.set_column('H:H', 16)
        #   Ft
        worksheet.write_column(self.__column_start_write, self.__col_Ft_solution_container, self._Ft_solution_container, format_2)
        worksheet.set_column('I:I', 16)

        try:
            for i in range(0, np.shape(np.array(self._u_P_solution_container))[1]):
                worksheet.write_column(2, 9+i, np.array(self._u_P_solution_container)[:, i], format_1)
        except:
            pass
        worksheet.set_column("J:M", 22)

        #   freeze first two rows
        worksheet.freeze_panes(2, 0)

        #   close file
        workbook.close()

        # for force in self._Fn_list:
        #     force.save_solution_data()

    def __excel_columns(self):
        """

        :return:
        """
        self.__workbook_name = "data"

        self.__column_start_write = self.__column_start = 2

        self.__col_step_num_solution_container = 0
        self.__col_status_container = 1
        self.__col_step_size_solution_container = 2
        self.__col_t_solution_container = 3
        self.__col_distance_solution_container = 4
        self.__col_dqn_solution_container = 5
        self.__col_dqt_solution_container = 6
        self.__col_Fn_solution_container = 7
        self.__col_Ft_solution_container = 8

        self.__cols_containers_list = [self.__col_step_num_solution_container,
                                       self.__col_status_container,
                                       self.__col_step_size_solution_container,
                                       self.__col_t_solution_container,
                                       self.__col_distance_solution_container,
                                       self.__col_dqn_solution_container,
                                       self.__col_dqt_solution_container,
                                       self.__col_Fn_solution_container,
                                       self.__col_Ft_solution_container]

    def read_file(self, _file_abs_path):
        """


        :param filename:
        :return:
        """
        #   load data from file if already exists
        _path, _file = os.path.split(_file_abs_path)
        self.filename, self._filetype = os.path.splitext(_file)

        if self._filetype == ".xlsx" or self._solution_filetype == ".xls":
            self._read_excel_file(_file_abs_path)

    def _read_excel_file(self, file_path):
        """

        :return:
        """
        #   set default properties for excel format of data
        self.__excel_columns()
        #   create workbook object
        workbook = xlrd.open_workbook(file_path)
        #   get data sheet from workbook object
        worksheet = workbook.sheet_by_name(self.__workbook_name)

        #   get contact model
        __contact_model_type_info = worksheet.cell(0, 0).value
        __contact_model_type_info = __contact_model_type_info.encode('ascii','ignore')
        self.contact_model.setType(__contact_model_type_info[__contact_model_type_info.index(": ")+2:__contact_model_type_info.index("\n")])

        #   size of first column
        self.__N_col0 = len(worksheet.col(0))-2

        for sol_container_name, col_sol_container in zip(self.__solution_container_names_list, self.__cols_containers_list):
            #   read column and returns values as a list
            _list = worksheet.col_values(col_sol_container,
                                         start_rowx=self.__column_start,
                                         end_rowx=self.__N_col0)

            # print _list
            #   convert list to numpy array
            __solution_container = np.array(_list)

            #   set values to attribute of object
            setattr(self, sol_container_name, __solution_container)

    def get_contact_force(self):
        """
        
        """
        for _force_n, _force_t in zip(self._Fn_list, self._Ft_list):
            _force_n.get_data()
            _force_t.get_data()

    def _paint_GL(self, step=None):
        """
        Display opengl data of contact. Contact point in LCS of a body and
        :return:
        """
        if step is None:
            if self.contact_detected:
                for u_CP in self.u_P_list_LCS:
                    glVertex3f(u_CP[0], u_CP[1], 0)
        else:
            try:
                _u_P_solution_container = np.array(self._u_P_solution_container)[step, :]
                for i, u_CP in enumerate(self.u_P_list):
                    _u_CP = _u_P_solution_container[2*i:2*i+2]
                    glVertex3f(_u_CP[0], _u_CP[1], 0)
            except:
                pass

    def _paing_GL_GCS(self, step=None):
        """
        Paint contact point in GCS
        :return:
        """
        if step is None:
            if self.u_P_GCS is not None and self.contact_detected:
                glVertex3f(self.u_P_GCS[0], self.u_P_GCS[1], 0)
        else:
            None

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

    def plot(self, x, y, color=None, label=None, str_id=None):
        """

        :param x: options = (delta)
        :param y: options = (Fn, dqn)
        :param color:
        :return:
        """
        _contact_indices = self.__contact()

        label, color = self._plot_info(color)

        if x == "delta":
            x_data = abs(self._distance_solution_container[_contact_indices])
        elif x == "t":
            x_data = self._t_solution_container[_contact_indices] - self._t_solution_container[_contact_indices][0]
        else:
            raise ValueError, "x not corrent!"

        if y == "Fn":
            y_data = abs(self._Fn_solution_container[_contact_indices]) / 0.015
        elif y == "dqn":
            y_data = self._dqn_solution_container[_contact_indices]
        else:
            raise ValueError, "y not corrent!"

        plt.plot(x_data, y_data, marker=None, color=color, label=r"%s"%label)
        # plt.loglog(x_data, y_data, basex=10, color=color, label=label)
        #   add label to graph
        # plt.text(abs(min(self._distance_solution_container[_contact_indices])), max(y_data), str_id+" "+label)

    def __contact(self):
        """

        :return:
        """
        #   indices where contact is present, value is 1
        _contact_indices = np.nonzero(self._status_container)#_Fn_solution_container, _status_container
        return _contact_indices

    def _plot_info(self, color=None):
        """

        :return:
        """
        #   label
        label = self.contact_model._type

        #   plot
        if color is None:
            color = self._color_GL
        else:
            color = color
        return label, color

if __name__ == "__main__":
    #   solution_data_contact_12_sol_12.xlsx
    #   solution_data_contact.xlsx
    #   solution_data_contact_12_sol.xlsx
    filename = "solution_data_contact_12_sol.xlsx"

    contact = Contact()

    fig = plt.figure(num=1, figsize=(8, 6), dpi=100, facecolor='w', edgecolor='k')
    ax = plt.subplot(111)
    ax.ticklabel_format(style='sci',scilimits=(-4,4), axis='both')

    contact.read_file(filename)
    # contact.plot_Fn_delta()

    contact.plot_Fn_time()

    plt.grid()
    plt.show()


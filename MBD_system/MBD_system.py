"""
Created on 19. feb. 2014

@author: luka.skrinjar
"""
import copy
import time
import logging
import logging.handlers
from pprint import pprint


import numpy as np
import dill


import dprj_file_
import body.list_of_bodies as list_of_bodies
import contact.list_of_contacts as list_of_contacts
import force.list_of_forces as list_of_forces
import joint.list_of_joints as list_of_joints
import spring.list_of_springs as list_of_springs
import measure.list_of_measures as list_of_measures
import variable.list_of_variables as list_of_variables
import motion.list_of_motions as list_of_motions
from MBD_system_items import *
from list_flatten import list_flatten
from body.body import Body
from body.body_rigid import RigidBody
from body.body_ground import GroundBody
from body.body_flexible import FlexibleBody
from body.point_mass import PointMass
from body.read_body_data_file.read_body_data_file import read_body_data_file
from solution_data.solution_data import SolutionData
from global_variables import GlobalVariables
from extract_from_dictionary_by_string_in_key import extract_from_dictionary_by_string_in_key
from q2q_body import q2q_body


class MBDsystem(MBDsystemItem):
    """
    classdocs
    """

    def __init__(self, MBD_file_abs_path=[], MBD_folder_abs_path=[], MBD_filename="Model_1", dys=None, parent=None):
        """
        Constructor creates a world (MBD system) a object has object parameters:
        - list of bodies
        - list of forces
        - list of joints
        - list of contacts
        - gen_body_list_func - function creates body system from 1 body
        - gen_joint_list_func - function creates joints of body system between bodies
        """
        super(MBDsystem, self).__init__(MBD_filename, parent)
        self.list_of_object_groups = ["Bodies", "Forces", "Joints", "Contacts", "Springs", "Motions", "Measures", "Variables", "Markers"]

        #   paths to files
        self.abs_path_to_bodies = None
        self.abs_path_to_joints = None
        self.abs_path_to_forces = None
        self.abs_path_to_motions = None
        self.abs_path_to_springs = None
        self.abs_path_to_contacts = None
        self.abs_path_to_measures = None
        self.abs_path_to_variables = None

        #   groups
        self.Bodies = None
        self.Joints = None
        self.Forces = None
        self.Motions = None
        self.Springs = None
        self.Contacts = None
        self.Measures = None
        self.Variables = None

        #   list of all paths
        self.abs_paths = []

        self.get_folder_and_filename(MBD_file_abs_path)

        #   dys application properties
        self.dys = dys
        self.simulation_control_widget = None

        #   simulation settings
        #   use baumgarte stabilization method - BSM
        self.use_BSM = False
        self.integrationMethod = "RKF45"
        self.integrationMethods = ["Euler", "RKF45" , "HHT-I3"]

        #   type of analysis, options:
        #   kinematic
        #   dynamic
        self.analysis_type = None

        #   monte carlo simulation settings
        self.monte_carlo = False
        self.monte_carlo_number_of_simulations = 4

        #   filetype settings
        #   project filetype
        self._project_filetype = ".dprj"

        #   solution properties
        #   predefined list
        self.solutions = []
        #   a pointer to solution data object that has been loaded
        self.loaded_solution = None
        #   filetype, options:
        #   .dat, (.xlsx, .csv - not implemented yet)
        self._solution_filetype = ".sol" #.sol,.xlsx
        self._solution_filetypes = [".sol",
                                    ".xlsx",
                                    ".dat"]
        self._soltuion_filetype_software = ["Solution file", "Excel file", "Text file"]
        #   options:
        #   discard
        #   overwrite
        #   save to new
        self._solution_save_options = "discard"

        #   extension of data files
        self._data_filetype = ".dat"
        
        #   parameters - used for monte carlo simulations
        self.parameters = AnalysisTreeItem(self._name)
        self.parameters_selected = [AnalysisTreeItem(self._name)]
        self.__parameters_by_type()

        #   angle units, deg or rad
        self._angle_units = "deg"  #    deg or rad
        self._input_angle_units = "deg"

        #    simulation settings
        self.q = []
        self.analysis_type = "dynamic"
        self.time = 0.
        self.step_num = 0
        self.t_n = None
        self.h = 0.
        self.stepsNumber = None
        self.Hmax = 1E-3
        self.Hmin = 1E-4
        self.Hcontact = 1E-8
        self.Houtput = 2 * self.Hmax
        self.absTol = 1E-8
        self.relTol = 1E-8
        #   error control
        self.errorControl = True
        #   epsilon tolerance for norm of newton differences vector
        self.TOL_dq_i = 1E-9
        #   epsilon tolerance for norm of constraint equations of MBS system - C(q, t)
        self.TOL_C_i = 1E-9
        #   tolerance for constraint equations - this is checked before numerical analysis starts
        self.TOL_C = 1E-9
        #   solver
        self.solver = None

        #   energy
        self._mechanical_energy = 0.
        self._kinetic_energy = 0.
        self._potential_energy = 0.
        self._elastic_strain_energy = 0.
        
        #    set and create logging file
        LOG_FILENAME = MBD_filename + '.log'
        print "LOG_FILENAME =", LOG_FILENAME
        self._log_file = LOG_FILENAME
        self.loadSolutionFileWhenFinished = False
        self.restoreInitialConditionsWhenFinished = False

        self._logger = logging.getLogger('DyS_logger')
        self._logger.addHandler(logging.StreamHandler())
        # logging.getLogger().addHandler(logging.StreamHandler())
        self._log_file_abs_path = os.path.normpath(os.path.join(MBD_folder_abs_path, self._log_file))
        
        #    visualization properties
        self._update_display_type = "dt"
        if self.t_n is not None:
            self._dt = self.t_n / 100
        else:
            self._dt = 1.
        self.updateEveryIthStep = 10

        #    check if log file exists and if not create it
        # if not os.path.exists(_log_file_abs_path):
        #     pass#open(_log_file_abs_path, 'w+')
        # else:
        self._log_file_abs_path = os.path.join(MBD_folder_abs_path, self._log_file)
        # print "_log_file_abs_path =", _log_file_abs_path
        open(self._log_file_abs_path, 'w+')

        #   set up logging
        hdlr = logging.FileHandler(self._log_file_abs_path, mode='w')
        formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
        hdlr.setFormatter(formatter)
        self._logger.addHandler(hdlr) 
        self._logger.setLevel(logging.INFO)

#        logging.basicConfig(filename=self._log_file)
#        self._logging = logging.getLogger("test")
#        self._logging.basicConfig(filename=self._log_file)
#        self._logging.setLevel(logging.DEBUG)
        
#        self._logging.info("testing")
#        self._logging = logging.basicConfig(filename=self._log_file, level=logging.INFO,)        
#         self._logger.critical("critical log")
#         self._logger.error('error log')
#         self._logger.warning('warning log')
#         self._logger.info('info log')
#         self._logger.debug('debug log')

        self.MBD_system_constructed = False
        self._name = MBD_filename
        self.saved = False
        
        if MBD_folder_abs_path != []:
            if os.path.isabs(MBD_folder_abs_path):
                self.MBD_folder_abs_path = os.path.abspath(MBD_folder_abs_path)
            else:
                self.MBD_folder_abs_path = os.path.abspath(MBD_folder_abs_path)

        GlobalVariables.MBDsystem_folder = self.MBD_folder_abs_path
        os.chdir(self.MBD_folder_abs_path)
        #    create groups as items (only for tree view hierarchy)
        self.__create_groups()

        if MBD_file_abs_path == []:
            for file in os.listdir(self.MBD_folder_abs_path):
                if file.endswith(".dprj"):
                    self._name = file
                    self.MBD_file_abs_path = os.path.join(os.path.abspath(self.MBD_folder_abs_path), file)
        else:
            self.MBD_file_abs_path = MBD_file_abs_path

        #   save screen shots options
        self.save_screenshots = False
        self.saved_screenshots_folder_name = []

        #    simulation-solver properties
        self.joint_list_counted = False
        self.q0_created = False
        self.q0 = None
        self.q0_filename = None
        #    check if absolute path to MBS folder is defined and if not, the current working directory is set to absolute path
#         if MBS_folder_abs_path == []:
#             folder_path = os.getcwd()
#             #    absolute path to folder with folder name
#             self.folder_abs_path = os.path.join(folder_path, self.MBS_folder_name)
        
        #    add ground name to the list of body names
        self.__list_of_body_names = []
        self.ground = None
        self.bodies = []
        self.number_of_bodies = len(self.bodies)
        self.forces = []
        self.springs = []
        self.joints = []
        self.contacts = []
        self.motions = []
        self.measures = []
        self.markers = []
        
        #    gravity properties
        self.gravity_magnitude = 9.83
        self.gravity_vector = np.array([0, -1, 0])
        self.gravity = self.gravity_magnitude * self.gravity_vector

        #   simulation time info
        self.start_time_simulation_info_in_sec_UTC = None
        self.end_time_simulation_info_in_sec_UTC = None

        #    properties dictionary attribute
        self.MBD_file_properties_dict = {}

        #   visualization settings
        self.scale = 1
        self.GCS_visible = True
        
        if self.MBD_file_abs_path is None:
            #   create ground body object
            self.create_ground()

            self.C_q_number_of_rows = 0
        else:
            self.construct_MBD_system(MBD_file_abs_path=self.MBD_file_abs_path, MBD_folder_abs_path=self.MBD_folder_abs_path, _name=self._name)

    def __getstate__(self):
        """

        :return:
        """
        if "dys" in self.__dict__:
            self.__dict__["dys"] = None

        if "_log_file_abs_path" in self.__dict__:
            self.__dict__["_log_file_abs_path"] = None

        if "_log_file" in self.__dict__:
            self.__dict__["_log_file"] = None

        return self.__dict__

    def read_ascii(self, filename):
        """
        Read .dprj project file ascii coded
        :param filename:
        :return:
        """

    def write(self, filename=None):
        """
        Save to file
        :return:
        """
        if filename is not None:
            self._name = filename

    def setSolver(self, solver):
        """

        :param solver:
        :return:
        """
        self.solver = solver

    def __parameters_by_type(self):
        """
        Function defines which parameters can be accesed via GUI by user for monte-carlo analysis
        :return:
        """
#         self.__body_parameters = ["max_penetration_depth",
#                                   "J_zz",
#                                   "_max_penetration_depth",
#                                   "skin_thickness",
#                                   "mass"]
#         self.__force_parameters = []
#         self.__spring_parameters = []
#         self.__joint_parameters = ["u_iP",
#                                    "u_jP"]
#         self.__contact_parameters = ["coef_of_friction_static",
#                                      "coef_of_friction_dynamic",
#                                      "contact_model",
#                                      "distance_TOL",
#                                      "_dq_t_TOL",
#                                      "coef_of_restitution",
#                                      "z_dim",
#                                      "contact_model_type",
#                                      "friction_model_type",
#                                      "friction_model"]

        self.__body_parameters = ["_name",
                                  "body_id",
                                  "color_GL",
                                  "mass"]
        
        self.__force_parameters = []
        self.__spring_parameters = []
        self.__joint_parameters = ["joint_id",
                                   "joint_type"]

        self.__contact_parameters = ["contact_id",
                                     "_type"]

        self.__contact_model_parameters = []
        self.__friction_model_parameters = []

    def parameters_construct_tree_view(self, item=None):
        """
        Function lists all parameters (attributes) of all objects in MBD system recursively
        """
        #   first run, the item uses self object attributes, and goes through its children, then it recursively runs through all children items
        #    and items that are inputed as arguments
        if item is None:
            item = self

        for child in item._children:
            if child._typeInfo == "group":
                _parent = self.parameters
                AnalysisTreeItem(child._name, parent=_parent)

                #   this is recursion to call this function but with different argument
                self.parameters_construct_tree_view(child)

            elif child._typeInfo == "body":
                if child._name != "Ground":
                    _parent = [child_ for child_ in self.parameters._children if child_._name == "Bodies"][0]
                    child_tree_item = AnalysisTreeItem(child._name, parent=_parent)
                    for attr in  child.__dict__:
                        if attr in self.__body_parameters:
                            AnalysisTreeItem(attr, item=child, mean_value=getattr(child, attr), parent=child_tree_item)
                        # child._children.append(attr_obj)

            elif child._typeInfo == "joint":
                _parent = [child_ for child_ in self.parameters._children if child_._name == "Joints"][0]
                child_tree_item = AnalysisTreeItem(child._name, parent=_parent)
                for attr in  child.__dict__:
                    if attr in self.__joint_parameters:
                        AnalysisTreeItem(attr, item=child, mean_value=getattr(child, attr), parent=child_tree_item)

            elif child._typeInfo == "contact":
                _parent = [child_ for child_ in self.parameters._children if child_._name == "Contacts"][0]
                child_tree_item = AnalysisTreeItem(child._name, parent=_parent)
                for attr in  child.__dict__:
                    if attr in self.__contact_parameters:
                        AnalysisTreeItem(attr, item=child, mean_value=getattr(child, attr), parent=child_tree_item)

#             else:
#                 _parent = self.parameters
#                 
#                 AnalysisTreeItem(child._name, parent=_parent)
#             
#                 self.parameters_construct_tree_view(child)
            
            
#             if child._typeInfo == "group":
#                 AnalysisTreeItem("group", item=child, parent=self)
#                 self.list_parameters(item=child)
#                 
#             if child._typeInfo == "body":
#                 AnalysisTreeItem(child._name, item=child, parent = child._parent)
#                 pass
                
#             if child._typeInfo == "body" or child._typeInfo == "joint":
#                 AnalysisTreeItem(child._name, item=child, parent = child._parent)
#             # if child._typeInfo == "group" or child._typeInfo == "MBDsystem" or child._typeInfo == "body":
#                 #    loop through all oboject's attributes if _typeInfo is right
#                 for attr in child.__dict__:
# #                     if attr[0] is not "_":
#                     print "%s = %s" % (attr, getattr(child, attr))
#                     #    create analysis tree item object for display in analysis tree view of parameters
#                     AnalysisTreeItem(attr, item=child, mean_value=getattr(child, attr), parent = child)# child._typeInfo, THIS IS OK -> (but not exccelent)self.parameters
#                      
#                 self.list_parameters(item=child)
#             else:
# #                 AnalysisTreeItem(attr, item=child, mean_value=getattr(child, attr), parent = self)
#                  
#                 self.list_parameters(item=self.parameters)
    
    def parameters_construct_table_view(self):
        """
        Function is used to display only selected tree items in table view
        """

    def __create_groups(self):
        for group_name in self.list_of_object_groups:
            group_obj = GroupItem(group_name, parent=self)
            setattr(self, group_name, group_obj)

            #    create abs path to group items filenames
            self._create_abs_path_to_group_file(group_name)

        group_obj = SettingsGroupItem("Settings", self)
        setattr(self, group_obj._name, group_obj)

        group_obj = SolutionGroupItem("Solution", self)
        setattr(self, group_obj._name, group_obj)

    def _create_abs_path_to_group_file(self, group_name):
        abs_path_by_group = "abs_path_to_" + group_name.lower()
        abs_path = os.path.join(self.MBD_folder_abs_path, group_name.lower() + self._data_filetype)
        setattr(self, abs_path_by_group, abs_path)

        self.abs_paths.append(abs_path)

    def get_properties(self):
        """
        Function returns MBD file properties to display
        """

    def get_folder_and_filename(self, MBD_file_abs_path):
        """
        
        """
        self.MBD_file_abs_path = MBD_file_abs_path
        if MBD_file_abs_path is not None:
            self.MBD_folder_abs_path, self.filename_ = os.path.split(MBD_file_abs_path)
        else:
            self.MBD_folder_abs_path = "c:\Temp"
            self.filename_ = "Model_1"
            
    def construct_MBD_system(self, MBD_file_abs_path="", MBD_folder_name="", MBD_folder_abs_path="", _name=""):
        """

        :param MBD_file_abs_path:
        :param MBD_folder_name:
        :param MBD_folder_abs_path:
        :param _name:
        :return:
        """
        self._name = _name
        self.MBD_folder_name = MBD_folder_name
        self.MBD_folder_abs_path = MBD_folder_abs_path
        self.MBD_file_abs_path = MBD_file_abs_path

        #   dictionary of MBD properties
        self.MBD_file_properties_dict = dprj_file_.dprj_file_read(self.MBD_file_abs_path)
        #    create ground body
        self.create_ground()
        #    create bodies
        self.create_bodies()
        #    create joints
        self.create_joints()
        #    create forces
        self.create_forces()
        #   create motions
        self.create_motions()
        #    create contacts
        self.create_contacts()
        #    create springs
        self.create_springs()
        #    count joints by type
        self.count_joints_by_type()
        #   create measures
        self.create_measures()
        #   create variables
        self.create_variables()

        #   add properties from dictionary to object
        self.dict_items_2_object_propeties(self.MBD_file_properties_dict)

        if self.save_screenshots and self.MBD_folder_abs_path != []:
            self.saved_screenshots_folder_name = "screenshots"
            self.saved_screenshots_folder_abs_path = os.path.join(self.MBD_folder_abs_path, self.saved_screenshots_folder_name)
            #    create folder
            try:
                os.stat(self.saved_screenshots_folder_abs_path)
            except:
                os.mkdir(self.saved_screenshots_folder_abs_path) 

        #   if file is defined overwrite initial conditions of MBD system from file
        if self.q0_filename is not None:
            if os.path.isfile(self.q0_filename):
                #   create solution object to read initial conditions file
                _sol = SolutionData(_file=self.q0_filename)
                _sol.read_file()
                #   assign initial conditions to MBD system object attribute
                self.q0 = _sol._q_sol_container[-1, :]

                #   delete object
                del _sol

                #   set initial conditions (q, dq) to each body object
                self.set_q(self.q0)

        else:
            self.evaluate_q()

        #   status of MBD system
        self.MBD_system_constructed = True

        #   gravity
        self.gravity = self.gravity_magnitude * self.gravity_vector

    def dict_items_2_object_propeties(self, _dict):
        for item in _dict:
            _value = _dict[item]
            self.add(item, _value)
    
    def add(self, key, value):
        setattr(self, key, value)
    
    def delete_MBD_system(self):
        """
        Delete MBD_system object additional properties.
        """
        self._name = "Model_1"
        self.MBD_folder_abs_path = []
        
        self.create_ground()

        for body in self.bodies:
            body.delete_VBO()

        self.bodies_file_list = []
        self.list_of_body_names = []
 
        self.bodies = []
        
        self.number_of_bodies = 0
        self.forces = []
        self.joints = []
        self.spring = []
        self.contacts = []
        
        self.MBD_system_constructed = False
        
        self.count_joints_by_type()

    def create_ground(self):
        """
        Create ground object in object MBD_system.
        """
        self.ground = GroundBody(parent=self)

        #   dictionary of properties for grid visualization object
        grid_properties_dict = extract_from_dictionary_by_string_in_key(self.MBD_file_properties_dict, "grid.")

        for key in grid_properties_dict:
            setattr(self.ground.grid, key, grid_properties_dict[key])

    def create_bodies(self, bodies_function=None):
        self.abs_path_to_bodies = os.path.join(self.MBD_folder_abs_path, "bodies" + self._data_filetype)
        if bodies_function is None:
            #    reads bodies list - file bodies.txt

            self.list_of_body_names = list_of_bodies.create_list(self.abs_path_to_bodies)

            #    create list of bodies
            for body_name in self.list_of_body_names:
                #     #    add body
                self.addBody(body_name)
        else:
            pass
            #TODO

        #    number of all bodies
        self.number_of_bodies = len(self.bodies)
 
        if self.number_of_bodies == 0:
            print "No bodies (objects) created when finished reading folder:", self.MBD_folder_name

        #   count body coordinates
        self.evaluate_q_i_size()

    def addBody(self, body_name):
        """
        Create a body object and add it to the list of bodies
        """
        file_path = os.path.join(self.MBD_folder_abs_path, body_name + self._data_filetype)

        #   read body data file as dict
        body_data_dict = read_body_data_file(file_path)

        #   create body
        if body_data_dict:
            if body_data_dict["body_type"] == "rigid body":
                body = RigidBody(name=body_name, file_path=file_path, _dict=body_data_dict, parent=self.Bodies)
                # self.bodies.append(body)

            elif body_data_dict["body_type"] == "flexible body":
                body = FlexibleBody(name=body_name, file_path=file_path, _dict=body_data_dict, parent=self.Bodies)
                # self.bodies.append(body)
                # body.check_slope_discontinuity(parent=self.Joints)
                for slope_discontinuity in body.mesh.slope_discontinuities:
                    self.joints.append(slope_discontinuity)

            elif body_data_dict["body_type"] == "point mass":
                body = PointMass(name=body_name, file_path=file_path, _dict=body_data_dict, parent=self.Bodies)

            else:
                raise "Body data filename %s for body %s not found"%(file_path, body_name)

            if body is not None:
                self.bodies.append(body)

        else:
            print "Body properties are not defined!"

        del body_data_dict

    def create_joints(self):
        """
        Create list of joint objects
        Supported types of joint:
        fixed
        revolute
        prismatic
        hinged support
        fixed support
        revolute joint rigid-flexible
        revolute joint flexible-flexible
        rigid joint flexible-flexible
        rigid joint rigid-flexible
        """
        self.joints = list_of_joints.create_list(self.joints, self.abs_path_to_joints, parent=self.Joints)

        if self.joints == [] and os.path.isfile(self.abs_path_to_joints):
            print "No joints (objects) created when finished reading filename: joints.txt. Check file: %s" % self.abs_path_to_joints

    def create_forces(self):
        """
        Create list of force objects
        """
        self.forces = list_of_forces.create_list(self.abs_path_to_forces, parent=self.Forces)

    def create_motions(self):
        """
        Create list of motion objects
        """
        self.motions = list_of_motions.create_list(self.abs_path_to_motions, parent=self.Motions)

    def create_springs(self):
        """
        Create list of spring objects
        """
        self.springs += list_of_springs.create_list(self.abs_path_to_springs, parent=self.Springs)

    def create_contacts(self):
        """
        Create list of contact objects
        """
        #    list of contact objects
        self.contacts = list_of_contacts.create_list(self.abs_path_to_contacts, parent=self.Contacts)

        #   if list is not empty prints lower info message
        if self.contacts:
            #    for each contact object in a list do
            logging.getLogger("DyS_logger").info("Creating AABB trees for every contact pair.")

        for contact in self.contacts:
            if contact._contact_type.lower() in ["roughness profile", "general"]:
                #    create bounding box
                contact.create_bounding_box_for_each_body_in_contact(bodies=self.bodies)

        if not self.contacts:
            logging.getLogger("DyS_logger").info("No AABB trees at contact pairs created!")
        else:
            logging.getLogger("DyS_logger").info("AABB trees for every contact pair created!")

    def create_measures(self):
        """
        Create list of measures objects to measure and display data during simulation
        :return:
        """
        #    list of measure objects
        self.measures = list_of_measures.create_list(self.abs_path_to_measures, parent=self.Measures)

    def create_variables(self):
        """

        :return:
        """
        self.variables = list_of_variables.create_list(self.abs_path_to_variables, parent=self.Variables)

    def evaluate_q0(self):
        """
        Function returns q vector (positions and velocities) of MBD system from body object data
        """
        _q = []
        #   q
        for body in self.bodies:
            q_body = body.evaluate_q0()
            _q.extend(q_body.tolist())

        #   dq
        for body in self.bodies:
            dq_body = body.evaluate_dq0()
            _q.extend(dq_body.tolist())

        self.q = np.array(_q).flatten()

        #   create initial conditions when this function if first called
        if self.q0 is None:
            self.q0 = copy.copy(self.q)

        return self.q

    def evaluate_q(self):
        """
        Function returns q vector (positions and velocities) of MBD system from body object data
        """
        _q = []
        #   q
        for i, body in enumerate(self.bodies):
            q_body = body.evaluate_q()
            _q.extend(q_body.tolist())

        #   dq
        for i, body in enumerate(self.bodies):
            dq_body = body.evaluate_dq()
            _q.extend(dq_body.tolist())

        self.q = np.array(_q).flatten()

        #   create initial conditions when this function if first called
        if self.q0 is None:
            self.q0 = copy.copy(self.q)

        return self.q

    def set_q(self, q):
        """

        :param q:
        :return:
        """
        self.update_coordinates_and_angles_of_all_bodies(t=0, q=q, step=0)

    def _restore_initial_conditions(self):
        """
        Restore initial conditions
        """
        # if not hasattr(self, "q0"):
        #     self.q0 = self.evaluate_q()
        self.update_coordinates_and_angles_of_all_bodies(t=0., q=self.q0)

        for i, contact in enumerate(self.contacts):
            contact.reset()

        for measure in self.measures:
            measure.reset()

        for i, force in enumerate(self.forces):
            force.reset(self.q0)

        for spring in self.springs:
            spring.reset(self.q0)

        for joint in self.joints:
            joint.reset(self.q0)

        self.update_vtk_data(0., self.q0)

    def update_vtk_data(self, t, q):
        """

        :return:
        """
        for body in self.bodies:
            body.step = self.step_num
            body.update_vtk_data(t, q)

            if body.AABBtree is not None:
                if hasattr(body.AABBtree, "vtk_actor"):
                    if body.AABBtree.vtk_actor is not None:
                        body.AABBtree.update_vtk_data(q)

        for i, force in enumerate(copy.copy(self.forces)):
            if not force.active and force.remove:
                self.forces.remove(force)

        for i, force in enumerate(self.forces):
            force.update_vtk_data(t, q)

        for spring in self.springs:
            spring.update_vtk_data(q)

        for joint in self.joints:
            joint.update_vtk_data(q)

        for contact in self.contacts:
            contact.update_vtk_data(q)

    def update_coordinates_and_angles_of_all_bodies(self, t, q, step=None):
        """

        :param q:
        :param step:
        :return:
        """
        #   set step value in global variables
        GlobalVariables.step = step

        q_, dq_ = self.q2positions_velocities(q)
        for body in self.bodies:
            # q_b = q_[3 * body.body_id:3 * body.body_id + 3]
            q_b = q2q_body(q, body.body_id)
            body.update_coordinates_and_angles_2D(q_b)

            body.update_vtk_data(t, q)

        for i, spring in enumerate(self.springs):
            spring.update_vtk_data(q)

    def update_velocities_of_all_bodies(self, q):
        """

        :param q:
        :return:
        """
        q_, dq_ = self.q2positions_velocities(q)
        for body in self.bodies:
            dq_b = dq_[3 * body.body_id:3 * body.body_id + 3]
            body.update_velocities_2D(dq_b)        

    def update_simulation_properties(self, time, step_num, h):
        """

        :param time:
        :param step_num:
        :return:
        """
        self.time = time
        self.step_num = step_num
        self.h = h

    def update_positions_and_velocities_of_all_bodies(self, q):
        """

        :param q:
        :return:
        """
        q_, dq_ = self.q2positions_velocities(q)
        
        for body in self.bodies:
            q_b = q_[3 * body.body_id:3 * body.body_id + 3]
            dq_b = dq_[3 * body.body_id:3 * body.body_id + 3]
            body.update_coordinates_and_angles_2D(q_b)
            body.update_velocities_2D(dq_b)  

    def q2positions_velocities(self, q):
        """
        
        """
        q = q[0:3 * self.number_of_bodies]
        dq = q[3 * self.number_of_bodies::]
        return q, dq

    def count_joints_by_type(self):
        """
        Count number of joints by type.
        """
        #   create vector of constrained coordinates per joint
        C_q_i_dim = []
        self.C_q_number_of_rows = 0
        for i, joint in enumerate(self.joints):
            C_q_i_dim.append(joint.C_q_dim)

            self.C_q_number_of_rows += joint.C_q_dim[0]

        GlobalVariables.C_q_i_dim = C_q_i_dim

        #   for rigid body
        #   number of fixed joints
        self.number_of_fixed_joints = sum(1 for joint in self.joints if joint.joint_type == "fixed")

        #   nuimber of revolute joints
        self.number_of_revolute_joints = sum(1 for joint in self.joints if joint.joint_type == "revolute")

        #   nuimber of prismatic joints
        self.number_of_prismatic_joints = sum(1 for joint in self.joints if joint.joint_type == "prismatic")

        #   for flexible body (ANCF mesh)
        #   number of fixed support joins
        self.number_of_fixed_supports = sum(1 for joint in self.joints if joint.joint_type == "fixed support")

        #   number of hinged support joints
        self.number_of_hinged_supports = sum(1 for joint in self.joints if joint.joint_type == "hinged support")

        #   number of roller support joints
        self.number_of_roller_supports = sum(1 for joint in self.joints if joint.joint_type == "roller support")

        #   number of rows of C_q matrix
        # self.C_q_number_of_rows = 2 * self.number_of_revolute_joints + 3 * self.number_of_fixed_joints + 2 * self.number_of_prismatic_joints

        self.joint_list_counted = True

    def evaluate_q_i_size(self):
        """
        Function creates a global vector that
        :return:
        """

        self.q_i_dim = []
        for body in self.bodies:
            _q_i_dim = body.evaluate_q_i_size()
            self.q_i_dim.append(_q_i_dim)

        self.q_i_dim = np.array(self.q_i_dim)

        #   set as global variable
        self.q_i_dim = np.array(self.q_i_dim)
        GlobalVariables.q_i_dim = self.q_i_dim

    def evaluate_M_size(self):
        """
        Matrix evaluates size of mass matrix of MBD system based on type of bodies
        :return:
        """
        self.M_dim = 0
        self.M_dim_rigid = 0
        self.M_dim_flexible = 0

        for i, body in enumerate(self.bodies):
            M_i = body.evaluate_M_size()
            if body.body_type == "rigid body":
                self.M_dim_rigid += M_i
            else:
                self.M_dim_flexible += M_i

            self.M_dim += M_i

        return self.M_dim

    def evaluate_C_number_of_rows(self):
        """
        Function evaluates number of rows of constraint vector equations, based on number of joints
        and number of motions, if kinematic analysis is done (for kinematically driven mechanical systems)
        :return:
        """
        self.C_number_of_rows_motions = 0
        _q = 0
        if self.analysis_type == "kinematic":
            for motion in self.motions:
                _q += motion.C_size

            self.C_number_of_rows_motions = self.C_number_of_rows_motions + _q

        C_q_number_of_rows = self.C_q_number_of_rows + self.C_number_of_rows_motions

        return C_q_number_of_rows

    def evaluate_Hmin(self):
        """
        Function evaluates Hmin based on mass properties and stiffness properties of the system
        :return:
        """
        if self.bodies != []:
            #   get min mass properties
            mass_properties = []
            for body in self.bodies:
                mass_properties.append(body.mass)
                mass_properties.append(body.J_zz)

            #   get minimum value
            m = min(mass_properties)

            #   get max stiffness properties
            stiffness_properties = []
            for spring in self.springs:
                stiffness_properties.append(spring.k)
            
            for contact in self.contacts:
                if hasattr(contact.contact_model, "K"):
                    if not isinstance(contact.contact_model, list):
                        stiffness_properties.append(contact.contact_model.K)

            stiffness_properties = list_flatten(stiffness_properties)
            
            #   get maximum value
            if stiffness_properties != [] and m > 0:
                
                k = max(stiffness_properties) 

                if k > 0:
                    t_min = 2*np.pi*(np.sqrt(k/m)**(-1))
                else:
                    self.Hmin = 1E-4
                    t_min = self.Hmin
                Hmin = (10**np.floor(np.log10(t_min)))*1E-1
                print "Suggested value of Hmin is %4.3e"%Hmin
                logging.getLogger("DyS_logger").info("Suggested value of Hmin is %4.3e" %Hmin)
                # QtGui.QMessageBox.information(self._3, "Information!", "Suggested value of Hmin is %s"%t_min,QtGui.QMessageBox.Ok, QtGui.QMessageBox.NoButton,QtGui.QMessageBox.NoButton)

    def evaluate_Hmax(self):
        """

        :return:
        """
        if self.t_n/10 < self.Hmax:
            self.Hmax = 0.01*self.t_n
            print "Hmax modified in %s"%self.__class__.__name__

        if self.Hcontact > self.Hmax:
            self.Hcontact = self.Hmax

    def evaluate_mechanical_energy(self, gravity=None, q=None):
        """

        :param q:
        :return:
        """
        if gravity is None:
            gravity = self.gravity

        #   kinetic energy
        self._kinetic_energy = self.evaluate_kinetic_energy(q, gravity)

        #   potential energy
        self._potential_energy = self.evaluate_potential_energy(q, gravity)

        #   elastic strain energy
        self._elastic_strain_energy = self.evaluate_elastic_strain_energy(q, gravity)

        #   total mechanical energy
        self._mechanical_energy = self._kinetic_energy + self._potential_energy + self._elastic_strain_energy

        return self._mechanical_energy, self._kinetic_energy, self._potential_energy, self._elastic_strain_energy

    def evaluate_kinetic_energy(self, q, gravity):
        """

        :return:
        """
        #    predefine zero array
        E_k_b = np.zeros(self.number_of_bodies)
        for i, body in enumerate(self.bodies):
            E_k_b[i] = body.evaluate_kinetic_energy(q=q)

        E_k = np.sum(E_k_b)

        return E_k

    def evaluate_potential_energy(self, q, gravity):
        """

        :return:
        """
        #   predefine zero array
        #   bodies
        E_p_b = np.zeros(self.number_of_bodies)
        for i, body in enumerate(self.bodies):
            E_p_b[i] = body.evaluate_potential_energy(q=q, gravity=gravity)

        #   springs
        E_p_s = np.zeros(len(self.springs))
        for i, spring in enumerate(self.springs):
            E_p_s[i] = spring.evaluate_potential_energy(q=q)

        #   contacts
        # E_p_c = np.zeros(len(self.contacts))
        # for i, contact in enumerate(self.contacts):
        #     E_p_c[i] = contact.evaluate_potential_energy(q=q)

        E_p = np.sum(E_p_b) + np.sum(E_p_s)# + np.sum(E_p_c)

        return E_p

    def evaluate_elastic_strain_energy(self, q, gravity):
        """

        :return:
        """
        #   predefine zero array
        #   bodies
        E_es_b = np.zeros(self.number_of_bodies)
        for i, body in enumerate(self.bodies):
            if hasattr(body, "evaluate_elastic_strain_energy") and hasattr(body, "_elastic_strain_energy"):
                E_es_b[i] = body.evaluate_elastic_strain_energy(q=q)

        E_p = np.sum(E_es_b)

        return E_p

    def testing(self):
        """

        :return:
        """
        print "testing()@",__name__
        print "N =", len(self.forces)

        # for contact_point in self.contacts[0]._contact_point_obj_list:
        #     contact_point.active = False
        #
        #
        # self.contacts[0].no_contact()
        # print "self.contacts[0]._Fn_list =", self.contacts[0]._Fn_list
        # print "self.contacts[0]._Ft_list =", self.contacts[0]._Ft_list
        #
        # for i, force in enumerate(self.forces):
        #     print "force =", force
        #     if not force.active and force.remove:
        #         self.forces.remove(force)
        #         print "force removed!", force, force._name
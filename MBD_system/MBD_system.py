"""
Created on 19. feb. 2014

@author: luka.skrinjar
"""
import os
import sys
import logging
import logging.handlers
from pprint import pprint


import numpy as np
from PyQt4 import QtCore, QtGui


from body.body import Body
import body.list_of_bodies as list_of_bodies
import contact.list_of_contacts as list_of_contacts
import dprj_file_
import force.list_of_forces as list_of_forces
import joint.list_of_joints as list_of_joints
import spring.list_of_springs as list_of_springs
from MBD_system_items import *
from motion import list_of_motions

class MBDsystem(MBDsystemItem):
    """
    classdocs
    """
    def __init__(self, MBD_file_abs_path=[], MBD_folder_name=[], MBD_folder_abs_path=[], MBD_filename="Model_1", parent=None):
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
        self.list_of_object_groups = ["Bodies", "Forces", "Joints", "Contacts", "Springs", "Motions", "Measures"]
        
        self.get_folder_and_filename(MBD_file_abs_path)

        #   simulation settings
        #   use baumgarte stabilization method - BSM
        self.use_BSM = True
        self.t_n = 1
        self.integrationMethod = None

        #   type of analysis, options:
        #   kinematic
        #   dynamic
        self.analysis_type = None

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
        #   options:
        #   discard
        #   overwrite
        #   save to new
        self._save_options = "discard"

        #   extension of data files
        self._data_filetype = ".dat"
        
        #   parameters - used for monte carlo simulations
        self.parameters = AnalysisTreeItem(self._name)
        self.parameters_selected = [AnalysisTreeItem(self._name)]
        self.__parameters_by_type()

        #   angle units, deg or rad
        self._angle_units = "deg"  #    deg or rad
        
        #    set and create logging file
        LOG_FILENAME = MBD_filename + '.log'
        self._log_file = LOG_FILENAME
        self.loadSolutionFileWhenFinished = False
        self.updateEveryIthStep = 10

        self._logger = logging.getLogger('DyS_logger')

        _log_file_abs_path = os.path.join(MBD_folder_abs_path, self._log_file)
#         print "_log_file_abs_path =", _log_file_abs_path
#         
#         if not os.path.exists(_log_file_abs_path):
#             pass
#         else:
#             _log_file_abs_path = os.path.join(os.getcwd(), self._log_file)
            
#         print "_log_file_abs_path =", _log_file_abs_path
        hdlr = logging.FileHandler(_log_file_abs_path, mode='w')
        formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
        hdlr.setFormatter(formatter)
        self._logger.addHandler(hdlr) 
        self._logger.setLevel(logging.INFO)
        
        
#        logging.basicConfig(filename=self._log_file)
#        self._logging = logging.getLogger("test")
#        self._logging.basicConfig(filename=self._log_file)
#        self._logging.setLevel(logging.DEBUG)
        
#        
        
#        self._logging.info("testing")
#        self._logging = logging.basicConfig(filename=self._log_file, level=logging.INFO,)        
#         self._logger.critical("critical log")
#         self._logger.error('error log')
#         self._logger.warning('warning log')
#         self._logger.info('info log')
#         self._logger.debug('debug log')

        self.MBD_system_constructed = False
        self._name = MBD_filename
        
        
        if MBD_folder_abs_path != []:
            if os.path.isabs(MBD_folder_abs_path):
                self.MBD_folder_abs_path_ = os.path.abspath(MBD_folder_abs_path)
            else:
                self.MBD_folder_abs_path_ = os.path.abspath(MBD_folder_abs_path)
        
        #    create groups as items (only for tree view hierarchy)
        self.__create_groups()
        
        
        if MBD_file_abs_path == []:
            for file in os.listdir(self.MBD_folder_abs_path_):
                if file.endswith(".dprj"):
                    self._name = file
                    self.MBD_file_abs_path = os.path.join(os.path.abspath(self.MBD_folder_abs_path_), file)
        else:
            self.MBD_file_abs_path = MBD_file_abs_path

        #   save screen shots options
        self.save_screenshots = False
        self.saved_screenshots_folder_name = []

        #    simulation-solver properties
        self.joint_list_counted = False
        self.q0_created = False
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
        
        #    gravity properties
        self.gravity = 9.81
        self.gravity_vector = np.array([0, -1, 0])

        #    properties dictionary attribute
        self.MBD_file_properties_dict = {}
        
        #    simulation settings
        self.time = 0
        self.step_num = 0
        self.Hmax = 0
        self.Hmin = 0
        self.absTol = 0
        self.relTol = 0
        #   epsilon tolerance for norm of newton differences vector
        self.TOL_dq_i = 1E-9
        #   epsilon tolerance for norm of constraint equations of MBS system - C(q, t)
        self.TOL_C_i = 1E-9
        #   tolerance for constraint equations - this is checked before numerical analysis starts
        self.TOL_C = 1E-9

        #   visualization settings
        self.scale_GL = 1
        self.GCS_visible = True
        
        if self.MBD_file_abs_path == []:
            #   create ground body object
            self.create_ground()

            self.C_q_number_of_rows = 0
        else:
            self.construct_MBD_system(MBD_file_abs_path=self.MBD_file_abs_path, MBD_folder_abs_path=self.MBD_folder_abs_path_, _name=self._name)

    def read_ascii(self, filename):
        """
        Read .dprj project file ascii coded
        :param filename:
        :return:
        """


    def write(self):
        """

        :return:
        """

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
#         print "item =", item
        #    loop through items 
#         print "children =", item._children
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
        setattr(self, "abs_path_to_" + group_name.lower(), os.path.join(self.MBD_folder_abs_path_, group_name.lower() + self._data_filetype))

    def get_properties(self):
        """
        Function returns MBD file properties to display
        """

    def get_folder_and_filename(self, MBD_file_abs_path):
        """
        
        """
        self.MBD_file_abs_path = MBD_file_abs_path
        if MBD_file_abs_path != []:
            self.MBD_folder_abs_path_, self.filename_ = os.path.split(MBD_file_abs_path)
        else:
            self.MBD_folder_abs_path_ = "c:\Temp"
            self.filename_ = "Model_1"
            
    def construct_MBD_system(self, MBD_file_abs_path=[], MBD_folder_name=[], MBD_folder_abs_path=[], _name=[]):

#        if self.MBD_system_constructed:
#            self.delete_MBD_system()
            
            
        self._name = _name
        self.MBD_folder_name = MBD_folder_name
        self.MBD_folder_abs_path_ = MBD_folder_abs_path
        self.MBD_file_abs_path = MBD_file_abs_path

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
        
        self.MBD_file_properties_dict = dprj_file_.dprj_file_read(self.MBD_file_abs_path)
#         print "self.MBD_file_properties_dict =", self.MBD_file_properties_dict

        self.dict_items_2_object_propeties(self.MBD_file_properties_dict)
    
    

        if self.save_screenshots and self.MBD_folder_abs_path_ != []:
            self.saved_screenshots_folder_name = "screenshots"
            self.saved_screenshots_folder_abs_path = os.path.join(self.MBD_folder_abs_path_, self.saved_screenshots_folder_name)
            #    create folder
            try:
                os.stat(self.saved_screenshots_folder_abs_path)
            except:
                os.mkdir(self.saved_screenshots_folder_abs_path) 


        self.MBD_system_constructed = True
        
    
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
        self.MBD_folder_abs_path_ = []
        
        self.ground = None
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
        self.ground = Body(parent=self, name="Ground", MBD_folder_abs_path=self.MBD_folder_abs_path_)

    def create_bodies(self, bodies_function=[]):
        if bodies_function == []:
            #    reads bodies list - file bodies.txt
            self.list_of_body_names = list_of_bodies.create_list(filename="bodies" + self._data_filetype, MBD_folder_abs_path=self.MBD_folder_abs_path_)
            #    reads the content of the folder and returns the list of files that have .stl file of geometry
#             list_of_body_names = self.__list_of_body_names + list_of_bodies.create_list(MBD_folder_abs_path=self.MBD_folder_abs_path_)
 
            #    create list of bodies
            for body_name_ in self.list_of_body_names:
                if body_name_.lower() == "ground":
                    None
                else:
                    #    add body
                    self.addBody(body_name_)
        else:
            None

        #    number of all bodies
        self.number_of_bodies = len(self.bodies)
 
        if self.number_of_bodies == 0:
            print "No bodies (objects) created when finished reading folder:", self.MBD_folder_name

    def create_VBOs(self):
        for body in self.bodies:
            body.create_VBO()

    def addBody(self, body_name_):
        """
        Create a body object and add it to the list of bodies
        """
        body_ = Body(name=body_name_, MBD_folder_abs_path=self.MBD_folder_abs_path_, parent=self.Bodies)
        self.bodies.append(body_)

    def create_joints(self):
        """
        Create list of joint objects
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
        self.springs = list_of_springs.create_list(self.abs_path_to_springs, parent=self.Springs)

    def create_contacts(self):
        """
        Create list of contact objects
        """
        #    list of contact objects
        self.contacts = list_of_contacts.create_list(self.abs_path_to_contacts, parent=self.Contacts)

        #    for each contact object in a list do
        for contact in self.contacts:
#            contact.build_2D_contact_geometry(bodies=self.bodies)

            if contact._type.lower() == "general":
                #    create bounding box
                contact.create_bounding_box_for_each_body_in_contact(bodies=self.bodies)

    def get_q(self):
        """
        Function returns q vector of MBD system from body object data
        """
        _q = []
        #   q
        for body in self.bodies:
            q_body = np.append(body.R[0:2], body.theta[-1])
            _q.append(q_body)

        #   dq
        for body in self.bodies:
            dq_body = np.append(body.dR[0:2], body.dtheta[-1])
            _q.append(dq_body)

        __q = np.array(_q)

        return __q.flatten()

    def evaluate_q0(self):
        self.q0 = np.zeros(6 * self.number_of_bodies)

        self.bodies.sort(key=lambda body: body.transparent_GL, reverse=False)

        for body in self.bodies:
            #    insert positions and rotations
            self.q0[3 * body.body_id:3 * body.body_id + 3] = np.append(body.R[[0, 1]], body.theta[-1])
            #    insert velocities (translational and rotational)
            self.q0[3 * self.number_of_bodies + 3 * body.body_id:3 * self.number_of_bodies + 3 * body.body_id + 3] = np.append(body.dR[[0, 1]], body.dtheta[-1])
        
        self.q0_created = True
        
        return self.q0

    def _restore_initial_conditions(self):
        """
        Restore initial conditions
        """
        if not hasattr(self, "q0"):
            self.q0 = self.create_q0()

        self.update_coordinates_and_angles_of_all_bodies(self.q0)
        for contact in self.contacts:
            contact._reset_to_initial_state()

    def update_coordinates_and_angles_of_all_bodies(self, q, step=None):
        """

        :param q:
        :param step:
        :return:
        """
        q_, dq_ = self.q2positions_velocities(q)
        for body in self.bodies:
            q_b = q_[3 * body.body_id:3 * body.body_id + 3]
            body.update_coordinates_and_angles_2D(q_b)

    def update_velocities_of_all_bodies(self, q):
        q_, dq_ = self.q2positions_velocities(q)
        for body in self.bodies:
            dq_b = dq_[3 * body.body_id:3 * body.body_id + 3]
            body.update_velocities_2D(dq_b)        
    
    
    def update_simulation_properties(self, time, step_num):
        self.time = time
        self.step_num = step_num
        
        
    def update_positions_and_velocities_of_all_bodies(self, q):
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
        #   number of fixed joints
        self.number_of_fixed_joints = sum(1 for joint in self.joints if joint.joint_type == "fixed")
        #   nuimber of revolute joints
        self.number_of_revolute_joints = sum(1 for joint in self.joints if joint.joint_type == "revolute")
        #   nuimber of prismatic joints
        self.number_of_prismatic_joints = sum(1 for joint in self.joints if joint.joint_type == "prismatic")
        
        #   number of rows of
        self.C_q_number_of_rows = 2 * self.number_of_revolute_joints + 3 * self.number_of_fixed_joints + 2 * self.number_of_prismatic_joints

        self.joint_list_counted = True

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
                    stiffness_properties.append(contact.contact_model.K)

            #   get maximum value
            if stiffness_properties != []:
                k = max(stiffness_properties)

                if k > 0:
                    t_min = 2*np.pi*(np.sqrt(k/m)**(-1))
                else:
                    self.Hmin = 1E-4
                    t_min = self.Hmin
                Hmin = (10**np.floor(np.log10(t_min)))*1E-1
                if t_min > self.Hmin:
                    print "Suggested value of Hmin is %4.3e"%Hmin
                    # QtGui.QMessageBox.information(self._3, "Information!", "Suggested value of Hmin is %s"%t_min,QtGui.QMessageBox.Ok, QtGui.QMessageBox.NoButton,QtGui.QMessageBox.NoButton)

                else:
                    self.Hmin = Hmin

    def evaluate_Hmax(self):
        """

        :return:
        """
        if self.MBD_system.t_n/10 < self.MBD_system.Hmax:
            self.MBD_system.Hmax = 0.01*self.MBD_system.t_n
            print "Hmax modified in %s"%self.__class__.__name__
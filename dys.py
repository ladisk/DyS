"""
Created on 30. sep. 2013

@author: lskrinjar
"""
from pprint import pprint
import ctypes
import logging
import logging.handlers
import os
import sys
import warnings
import ConfigParser
import time
import copy
import sip
import socket

import dill
import numpy as np
from PyQt4 import QtCore, QtGui

try:
    from binaryornot.check import is_binary
except:
    pass


import qrc_resources
from MBD_system import MBD_system_items
from MBD_system import convert_bytes_to_
from MBD_system import read_and_write
from MBD_system.MBD_system import MBDsystem
from job_list_widget.job_list_widget import JobListWidget
from preferences_widget.preferences_widget import PreferencesWidget
from simulation_control_widget.simulation_control_widget import SimulationControlWidget
from tree_view_widget.tree_view_widget import TreeViewWidget
# from python_console.python_console import PythonConsole


np.set_printoptions(precision=12, threshold=None, edgeitems=100, linewidth=1000, suppress=False, nanstr=None, infstr=None)
# sys.setrecursionlimit(3000)
# print sys.getrecursionlimit()


try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s


#   authorship information
__authors__ = ['Luka Skrinjar', 'Ales Turel', 'Janko Slavic']
__author__ = ','.join(__authors__)
__credits__ = []
__copyright__ = 'Copyright (c) 2015'
__license__  = 'GPL'


# maintenance information
__maintainer__ = 'Luka Skrinjar'
__email__  = 'skrinjar.luka@gmail.com'
__version__ = "0.0.6"


filename = os.path.relpath(__file__)
current_working_directory = os.getcwd()
job_data_folder = ""


warnings.filterwarnings("ignore")


class DySMainWindow(QtGui.QMainWindow):
    """
    classdocs
    """
    
    def __init__(self):
        """
        Class constructor

        >>> 1 + 1
        2
        >>>
        """
        super(DySMainWindow, self).__init__()
        # self.setAttribute(QtCore.Qt.WA_TranslucentBackground, 90)
        #    minimum size
        self.setMinimumSize(500, 500)

        #   computer information
        self._hostname = socket.gethostname()
        print "self._hostname =", self._hostname
        #    style
        # self.setStyle("windows")
        #     plastique
        #     cde
        #     motif
        #     sgi
        #     windows
        #     cleanlooks
        #     mac
        
        #    icon
        self.setWindowIcon(QtGui.QIcon(":/application.png"))
        QtGui.QSystemTrayIcon(QtGui.QIcon(":/application.png"))
        
        #    location of window on self.screen_geometry
        self.window_offset_x = 250
        self.window_offset_y = 50

        #    language settings
        self.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.UnitedStates))

        #    window name
        self.setWindowTitle("DyS")

        #   font
        self.font = QtGui.QFont("Consolas", 10)
        
        #    mouse tracking
        self.setMouseTracking(True)
        
        #    context menu
        self.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        
        #    main window flags
        self.flags = QtCore.Qt.Window

        __main_dir = current_working_directory
        self._working_directory = QtCore.QString(os.path.abspath(__main_dir).replace("/", "\\"))

        #   get screen geometry size
        self.screen_geometry = QtGui.QDesktopWidget().screenGeometry()
        #    move
        # self.move(self.window_offset_x, self.window_offset_y)
        dy = .05*self.screen_geometry.height()

        #   size of application
        #   position main widget on self.screen_geometry
        self.setGeometry(.25*self.screen_geometry.width(), dy, .5*self.screen_geometry.width(), .5*self.screen_geometry.width())

        #    MBD system
        MBD_folder_name_ = None
        self.MBD_file_abs_path = None

        # self.MBD_filename = "0_0_0_double_pendulum.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_0_0_double_pendulum"

        # self.MBD_filename = "dys_0_.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_"

        # self.MBD_filename = "dys_0_2_1.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_2_1"

        # self.MBD_filename = "dys_0_3_prismatic_joint.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_3_prismatic_joint"

        # self.MBD_filename = "dys_0_4.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_4"

        # self.MBD_filename = "dys_0_5.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_5"

        # self.MBD_filename = "dys_0_6_1.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_6_1"

        # self.MBD_filename = "dys_0_6_2.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_6_2"

        # self.MBD_filename = "dys_0_7_0_1.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_0_1"

        # self.MBD_filename = "dys_0_7_0_2.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_0_2"

        # self.MBD_filename = "0_7_3_3_kinematic_analysis_1arm.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_3_3_kinematic_analysis_1arm"

        # self.MBD_filename = "0_7_3_3_kinematic_analysis_2arm.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_3_3_kinematic_analysis_2arm"

        # self.MBD_filename = "0_7_3_3_kinematic_analysis_slider_crank.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_3_3_kinematic_analysis_slider_crank"

        # self.MBD_filename = "0_7_3_3_kinematic_analysis_4bar.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_3_3_kinematic_analysis_4bar"

        # self.MBD_filename = "0_7_3_3_kinematic_analysis_efi.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_3_3_kinematic_analysis_efi"

        # self.MBD_filename = "0_7_3_3_dynamic_analysis_efi.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_3_3_dynamic_analysis_efi"
    
        # self.MBD_filename = "0_7_3_3_dynamic_analysis_efi_V2.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_3_3_dynamic_analysis_efi_V2"

        # self.MBD_filename = "dys_0_7_3_0.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_3_0"
        
        # self.MBD_filename = "0_7_3_0_plane_sphere.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_3_0_plane_sphere"

        # self.MBD_filename = "dys_0_7_2_revolute_clearence_joint.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_2_revolute_clearence_joint"

        # self.MBD_filename = "0_7_2_revolute_clearance_joint_1.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_2_revolute_clearance_joint_1"

        # self.MBD_filename = "0_7_3_0_contact_models.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_3_0_contact_models"

        self.MBD_filename = "0_7_3_0_contact_models_cylinder.dprj"
        MBD_folder_name_ = "dynamic_systems\\0_7_3_0_contact_models_cylinder"

        # self.MBD_filename = "0_7_3_2_pin_slot_clearance_joint.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_3_2_pin_slot_clearance_joint"

        # self.MBD_filename = "0_7_3_2_pin_slot_clearance_joint_2.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_3_2_pin_slot_clearance_joint_2"

        # self.MBD_filename = "0_7_3_2_pin_slot_clearance_joint_paper.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_3_2_pin_slot_clearance_joint_paper"

        # self.MBD_filename = "0_7_3_2_pin_slot_clearance_joint_paper_1.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_3_2_pin_slot_clearance_joint_paper_1"

        # self.MBD_filename = "0_7_3_2_pin_slot_clearance_joint_paper_2.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_3_2_pin_slot_clearance_joint_paper_2"

        # self.MBD_filename = "0_7_1_0_spring_translational.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_1_0_spring_translational"

        # self.MBD_filename = "0_8_0_surface_roughness.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_8_0_surface_roughness"

        # self.MBD_filename = "0_9_0_brush_slip-ring.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_9_0_brush_slip-ring"

        # self.MBD_filename = "0_9_1_brush_slip-ring.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_9_1_brush_slip-ring"

        # self.MBD_filename = "0_9_2_multiple_contacts.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_9_2_multiple_contacts"

        # self.MBD_filename = "0_9_3_brush_slip-ring.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_9_3_brush_slip-ring"

        # self.MBD_filename = "0_9_4_brush_slip-ring.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_9_4_brush_slip-ring"

        # self.MBD_filename = "0_9_5_brush_energy.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_9_5_brush_energy"

        # self.MBD_filename = "0_9_6_brush_slip-ring.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_9_6_brush_slip-ring"

        # self.MBD_filename = "1_0_0_ancf-beam.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_0_ancf-beam"

        # self.MBD_filename = "1_0_0_ancf-cantilever_1element.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_0_ancf-cantilever_1element"

        # self.MBD_filename = "1_0_0_ancf-cantilever_2elements.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_0_ancf-cantilever_2elements"

        # self.MBD_filename = "1_0_0_ancf-cantilever_2elements.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_0_ancf-cantilever_2elements"
        #
        # self.MBD_filename = "1_0_0_ancf-cantilever_4elements.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_0_ancf-cantilever_4elements"
        #
        # self.MBD_filename = "1_0_0_ancf-cantilever_8elements.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_0_ancf-cantilever_8elements"
        #
        # self.MBD_filename = "1_0_0_ancf-cantilever_16elements.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_0_ancf-cantilever_16elements"
        #
        # self.MBD_filename = "1_0_0_ancf-cantilever_32elements.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_0_ancf-cantilever_32elements"

        # self.MBD_filename = "1_0_0_ancf-free_falling_flexible_pendulum.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_0_ancf-free_falling_flexible_pendulum"

        # self.MBD_filename = "1_0_0_ancf-cantilever_M.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_0_ancf-cantilever_M"

        # self.MBD_filename = "1_0_0_ancf-beam_contact.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_0_ancf-beam_contact"

        # self.MBD_filename = "1_0_1_ancf_frame.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_1_ancf_frame"

        # self.MBD_filename = "1_0_1_geometry_test.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_1_geometry_test"

        # self.MBD_filename = "1_0_1_external_force.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_1_external_force"

        # self.MBD_filename = "1_0_1_paper_2_eigenfrequencies.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_1_paper_2_eigenfrequencies"

        # self.MBD_filename = "1_0_1_paper_2_dynamic_system_point_mass.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_1_paper_2_dynamic_system_point_mass"

        # self.MBD_filename = "1_0_1_paper_2_dynamic_system_2.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_1_paper_2_dynamic_system_2"

        # self.MBD_filename = "1_0_1_paper_2_dynamic_system_sub.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_1_paper_2_dynamic_system_sub"

        # self.MBD_filename = "1_0_1_paper_2_dynamic_system_2.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_1_paper_2_dynamic_system_2"

        # self.MBD_filename = "1_0_1_paper_2_dynamic_system_2_y_spring=-1mm.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_1_paper_2_dynamic_system_2_y_spring=-1mm"

        # self.MBD_filename = "1_0_1_paper_2_dynamic_system_2_y_spring=+1mm.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_1_paper_2_dynamic_system_2_y_spring=+1mm"

        # self.MBD_filename = "1_0_1_paper_2_dynamic_system_2_prestress_2x.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_1_paper_2_dynamic_system_2_prestress_2x"

        #self.MBD_filename = "1_0_1_paper_2_dynamic_system_2_prestress_params_2x.dprj"
        #MBD_folder_name_ = "dynamic_systems\\1_0_1_paper_2_dynamic_system_2_prestress_params_2x"

        # self.MBD_filename = "1_0_1_paper_2_dynamic_system_2_prestress_0.5x.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_1_paper_2_dynamic_system_2_prestress_0.5x"

        # self.MBD_filename = "1_0_1_paper_2_dynamic_system_rigid.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_1_paper_2_dynamic_system_rigid"

        # self.MBD_filename = "1_0_1_paper_2_dynamic_system_rigid_2.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_1_paper_2_dynamic_system_rigid_2"

        # self.MBD_filename = "0_7_3_2_pin_slot_clearance_joint_testing.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_3_2_pin_slot_clearance_joint_testing"

        # self.MBD_filename = "1_0_1_paper_2_rigid_joint_rigid_flexible.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_1_paper_2_rigid_joint_rigid_flexible"

        # self.MBD_filename = "1_0_1_paper_2_revolute_joint_flexible_flexible.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_1_paper_2_revolute_joint_flexible_flexible"

        # self.MBD_filename = "1_0_1_paper_2_revolute_joint_rigid_flexible.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_1_paper_2_revolute_joint_rigid_flexible"

        # self.MBD_filename = "1_0_1_paper_2_rigid_joint_point_mass_flexible.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_1_paper_2_rigid_joint_point_mass_flexible"

        # self.MBD_filename = "1_0_1_paper_2_rigid_joint_point_mass_flexible_offset.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_1_paper_2_rigid_joint_point_mass_flexible_offset"

        # self.MBD_filename = "1_0_1_paper_2_rigid_joint_rigid_flexible.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_1_paper_2_rigid_joint_rigid_flexible"

        # self.MBD_filename = "1_0_1_paper_2_rigid_joint_rigid_flexible_2.dprj"
        # MBD_folder_name_ = "dynamic_systems\\1_0_1_paper_2_rigid_joint_rigid_flexible_2"

        #   predefine attributes
        self.simulation_control_widget = None
        self._visualization_widget = "vtk"

        #   file is not defined
        if MBD_folder_name_ is None:
            self.MBD_folder_abs_path = os.getcwd()
            self.MBD_filename = project_filename = "Model_1"
            print "File is not defined"

        #   file is defined
        else:
            self._working_directory = os.path.join(os.getcwd(), "..")
            self.MBD_file_abs_path = os.path.abspath(os.path.join(self._working_directory, MBD_folder_name_, self.MBD_filename))

            self.MBD_folder_abs_path = os.path.join(self._working_directory, MBD_folder_name_)
            project_filename = os.path.basename(self.MBD_folder_abs_path)

        self._file_types = "Dynamic System Project File (*.dprj);;Multibody System File (*.mbd)"
        self.project = MBD_system_items.MBDsystemItem("MBDProject", parent=None)
        self.project.dys = self

        # print "self.MBD_folder_abs_path =", self.MBD_folder_abs_path
        self.MBD_system = MBDsystem(MBD_file_abs_path=self.MBD_file_abs_path,
                                    MBD_folder_abs_path=self.MBD_folder_abs_path,
                                    MBD_filename=project_filename,
                                    dys=self,
                                    parent=self.project)

        #   tree view widget
        self.tree_view_widget = TreeViewWidget(self.project, self.MBD_system, parent=self)
        self.tree_view_widget.setWindowFlags(self.flags)
        # self.tree_view_widget.show()
        #    move
        self.tree_view_widget.move(.18*self.screen_geometry.width() - self.tree_view_widget.frameGeometry().width(), dy)#self.frameGeometry().width(), self.frameGeometry().height())
        self.tree_view_widget.resize(.20*self.screen_geometry.width(), 0.8 * self.screen_geometry.height())

        #    simulation control widget
        self.simulation_control_widget = SimulationControlWidget(MBD_system=self.MBD_system, parent=self)
        self.simulation_control_widget.setWindowFlags(self.flags)
        # self.simulation_control_widget.show()

        #    move
        self.simulation_control_widget.move(.25*self.screen_geometry.width() + self.geometry().width(), dy)

        #   preferences widget
        self.preferences_widget = PreferencesWidget(self.simulation_control_widget, parent=None)

        #   graph widget - test
        for measure in self.MBD_system.measures:
            measure.graph_widget.setWindowFlags(self.flags)
            measure.graph_widget.show()

        #   set central widget
        self.setCentralWidget(self.simulation_control_widget.vtkWidget)

            # self.simulation_control_widget.vtkWidget.displayText()
            # print "2 =", self.simulation_control_widget.vtkWidget.vtkWidget.GetSize()
        # if self.simulation_control_widget.opengl_widget is not None:
        #     self.setCentralWidget(self.simulation_control_widget.opengl_widget)
        #    init graph widget - test
#         graph_widget = GraphWidget(parent=self)
#         graph_widget.setWindowFlags(self.flags)
#         graph_widget.show()
        #    output widget
        #    maybe put in new thread? if output will be large
#         self.OutputWidget = OutputWidget(parent=self)
#         self.OutputWidget.setWindowFlags(self.flags)
#         self.OutputWidget.show()

        #    job list widget
        self.JobListWidget = JobListWidget(job_list_=["job1", "job2", "job3"], parent=self)
        self.JobListWidget.setWindowFlags(self.flags)
        
        #    python console
        self.python_console = None#PythonConsole()

        #    resize
        # self.resize(self.simulation_control_widget.opengl_widget.initial_window_width, self.simulation_control_widget.opengl_widget.initial_window_height)
        
        #    create actions
        self.create_actions()
        
        #    create menus
        self.create_menus()
        
        #    create status bar
        self.create_status_bar()

        #    signals and connections
        self.connect_signals()

    def load_settings(self, configfile='settings.ini'):
        """
        Load settings from .ini file.
        ;param configfile: Path to .ini configuration file.
        """
        self.settings = ConfigParser.ConfigParser()
        self.settings.read(configfile)

        #   example
        # self.conv_tol = float(Config.get('DIC', 'convergence_tolerance'))
        # self.max_iter = int(Config.get('DIC', 'max_iterations'))
        # self.int_order = int(Config.get('DIC', 'interpolation_order'))
        # self.sequence_increment = int(Config.get('DIC', 'sequence_increment'))
        # self.min_roi_size = tuple(int((vel)) for vel in Config.get('GUI', 'min_roi_size').split(', '))
        # self.initial_roi_size = tuple(int((vel)) for vel in Config.get('GUI', 'initial_roi_size').split(', '))

    def connect_signals(self):
        """

        :return:
        """
        if self.simulation_control_widget.solver.analysis is not None:
            self.simulation_control_widget.solver.analysis.step_signal.signal_step.connect(self.update_statusbar_text2)

            #   simulation status
            self.simulation_control_widget.solver.analysis.signal_simulation_status.signal_simulation_status.connect(self.update_statusbar_simulation_status)

            self.simulation_control_widget.solver.analysis.energy_signal.signal_energy.connect(self.update_statusbar_text1)

        self.simulation_control_widget.signal_simulation_status.signal_simulation_status.connect(self.update_statusbar_simulation_status)
        self.simulation_control_widget.step_num_signal.signal_step.connect(self.update_statusbar_text2)

        self.tree_view_widget.filename_signal.signal_filename.connect(self.simulation_control_widget.load_solution_file)

        self.simulation_control_widget.energy_signal.signal_energy.connect(self.update_statusbar_text1)

        #   save file
        # self.connect(QtGui.QShortcut(QtGui.QKeySequence.Save, self), QtCore.SIGNAL('activated()'), self.saveFile)

    def create_actions(self):
        """
        Actions
        :return:
        """
        #    file - menu
        #    new
        self.newAction = QtGui.QAction('New', self, shortcut='Ctrl+N', statusTip='Create new')
        self.newAction.triggered.connect(self.createNewProject)
        #    open
        self.openAction = QtGui.QAction('Open', self, shortcut='Ctrl+O', statusTip='Open')
        self.openAction.triggered.connect(self.showOpenFileDialog)
        #    close
        self.closeAction = QtGui.QAction('Close Project', self, statusTip='Close Project')
        self.closeAction.triggered.connect(self.closeProjectFile)
        #   save
        self.saveAction = QtGui.QAction('Save', self)
        self.saveAction.setStatusTip("Save file")
        self.saveAction.setShortcut("Ctrl+S")
        self.saveAction.triggered.connect(self.saveFile)

        #    exit
        self.exitAction = QtGui.QAction('Exit', self, shortcut='Alt+F4', statusTip='Exit Application')
        self.exitAction.triggered.connect(self.close)

        #    view - menu
        self.take_snapshot_action = QtGui.QAction('Take snap shot',
                                                    self,
                                                    statusTip='Take snap shot')
        self.take_snapshot_action.triggered.connect(self.take_snapshot)
        
        #    window - menu
        #    window
        self.show_control_panel_action = QtGui.QAction(self.tr("&Control Panel"),
                                                       self,
                                                       statusTip='Control Panel')
        self.show_control_panel_action.triggered.connect(self.show_control_panel)
        
        self.show_tree_view_action = QtGui.QAction(self.tr("&Tree View"),
                                                       self,
                                                       statusTip='Tree View')
        self.show_tree_view_action.triggered.connect(self.show_tree_view)
        
        self.show_output_widget_action = QtGui.QAction(self.tr("&Output Window"),
                                                       self,
                                                       statusTip='Output Window')
        self.show_output_widget_action.triggered.connect(self.show_output_widget)
        
        self.show_job_list_action = QtGui.QAction(self.tr("&Job List"),
                                                       self,
                                                       statusTip='Job List')
        self.show_job_list_action.triggered.connect(self.show_job_list)
        
        self.show_python_console_action =  QtGui.QAction(self.tr("&Python console"),
                                                       self,
                                                       statusTip='Python console')
        self.show_python_console_action.triggered.connect(self.show_python_console_widget)

        #    settings - menu
        self.show_preferences_widget_action = QtGui.QAction('Preferences',
                                                                   self,
                                                                   shortcut="None",
                                                                   statusTip='Preferences')
        self.show_preferences_widget_action.triggered.connect(self.preferences_widget._show)
        
        #    help - menu
        #    about
        self.aboutAction = QtGui.QAction(self.tr("&About"), self, statusTip='Information about the application for simulation of multibody dynamics.')
        self.aboutAction.triggered.connect(self.about)

    def keyPressEvent(self, event):
        """
        User key pressed events
        """
        if type(event) == QtGui.QKeyEvent:
            if event.key() == QtCore.Qt.Key_Right and self.simulation_control_widget.ui.forwardButton.isEnabled():
                self.simulation_control_widget.animation_forward()
            
            if event.key() == QtCore.Qt.Key_Left and self.simulation_control_widget.ui.backwardButton.isEnabled():
                self.simulation_control_widget.animation_backward()

    def create_menus(self):
        """
        Create menus
        :return None:
        """
        #    file menu
        self.fileMenu = QtGui.QMenu(self.tr("&File"), self)
        
        #    create new
        createNew = QtGui.QAction('New', self)
        createNew.setShortcut('Ctrl+N')
        createNew.setStatusTip('Create New Project')
        self.fileMenu.addAction(createNew)
        self.fileMenu.addSeparator()

        #    open
        showOpenFileDialog = QtGui.QAction('Open', self)
        showOpenFileDialog.setShortcut('Ctrl+O')
        showOpenFileDialog.setStatusTip('Open Project')
        showOpenFileDialog.triggered.connect(self.showOpenFileDialog)
        self.fileMenu.addAction(showOpenFileDialog)
        self.fileMenu.addSeparator()

        #    save file
        saveFile = QtGui.QAction('Save', self)
        saveFile.setShortcut('Ctrl+S')
        saveFile.setStatusTip('Save Project')
        saveFile.triggered.connect(self.saveFile)
        self.fileMenu.addAction(saveFile)

        #    save as file
        saveAsFile = QtGui.QAction('Save As', self)
        saveAsFile.setShortcut('Ctrl+Shift+S')
        saveAsFile.setStatusTip('Save As')
        saveAsFile.triggered.connect(self.showSaveAsFileDialog)
        self.fileMenu.addAction(saveAsFile)
        self.fileMenu.addSeparator()

        #    close
        closeFile = QtGui.QAction('Close Project', self)
        closeFile.setStatusTip('Close Project')
        closeFile.triggered.connect(self.closeProjectFile)
        self.fileMenu.addAction(closeFile)

        #    exit
        self.fileMenu.addSeparator()
        self.fileMenu.addAction(self.exitAction)

        #    view menu
        self.viewMenu = QtGui.QMenu(self.tr("&View"), self)
        if self.simulation_control_widget.vtkWidget is not None:
            viewActionFront = self.viewMenu.addAction("Front")
            viewActionFront.triggered.connect(self.simulation_control_widget.vtkWidget._viewFront)

            viewActionBack = self.viewMenu.addAction("Back")
            viewActionBack.triggered.connect(self.simulation_control_widget.vtkWidget._viewBack)

            viewActionBottom = self.viewMenu.addAction("Bottom")
            viewActionBottom.triggered.connect(self.simulation_control_widget.vtkWidget._viewBottom)

            viewActionTop = self.viewMenu.addAction("Top")
            viewActionTop.triggered.connect(self.simulation_control_widget.vtkWidget._viewTop)

            viewActionLeft = self.viewMenu.addAction("Left")
            viewActionLeft.triggered.connect(self.simulation_control_widget.vtkWidget._viewLeft)

            viewActionRight = self.viewMenu.addAction("Right")
            viewActionRight.triggered.connect(self.simulation_control_widget.vtkWidget._viewRight)

            self.viewMenu.addSeparator()
            viewActionIsometric = self.viewMenu.addAction("Isometric")
            viewActionIsometric.triggered.connect(self.simulation_control_widget.vtkWidget._viewIsometric)

            self.viewMenu.addSeparator()
            viewActionIsometric = self.viewMenu.addAction("Save snapshot")
            viewActionIsometric.triggered.connect(self.simulation_control_widget.vtkWidget._saveSnapShot)
        
        #    window
        self.windowMenu = QtGui.QMenu(self.tr("&Window"), self)
        self.windowMenu.addAction(self.show_control_panel_action)
        self.windowMenu.addAction(self.show_tree_view_action)
        self.windowMenu.addAction(self.show_output_widget_action)
        self.windowMenu.addAction(self.show_job_list_action)
        self.windowMenu.addAction(self.show_python_console_action)
        
        #    settings
        self.settingsMenu = QtGui.QMenu(self.tr("&Settings"), self)
        self.settingsMenu.addAction(self.show_preferences_widget_action)
        
        #    about
        self.helpMenu = QtGui.QMenu(self.tr("&Help"), self)
        self.helpMenu.addAction(self.aboutAction)
        
        #    create menus
        self.menuBar().addMenu(self.fileMenu)
        self.menuBar().addMenu(self.viewMenu)
        self.menuBar().addMenu(self.settingsMenu)
        self.menuBar().addMenu(self.windowMenu)
        self.menuBar().addMenu(self.helpMenu)

    def viewActionFront(self):
        """

        :return:
        """
    def viewActionBack(self):
        """

        :return:
        """
    def viewActionBottom(self):
        """

        :return:
        """
    def viewActionTop(self):
        """

        :return:
        """
    def viewActionLeft(self):
        """

        :return:
        """
    def viewActionRight(self):
        """

        :return:
        """
    def viewActionIsometric(self):
        """

        :return:
        """

    def newFile(self):
        a = 1

    def update_data(self):
        """
        
        """
        self.simulation_control_widget.opengl_widget.update_data(self.MBD_system.bodies)
        self.simulation_control_widget.opengl_widget.repaintGL()

    def createNewProject(self):
        """
        
        """

    def closeProjectFile(self):
        """
        
        """
        self.MBD_system.delete_MBD_system()
        self.tree_view_widget.repaint()

    def showSaveAsFileDialog(self):
        """
        Save project under new filename - opens a save file dialog
        """
        saveAs_dialog = QtGui.QFileDialog()
        filename, file_type = saveAs_dialog.getSaveFileNameAndFilter(self, 'Save File', self.MBD_system.MBD_folder_abs_path, ("DyS Project Files (*.dprj)"))

        read_and_write.write_pickle(self.MBD_system, filename)
        # read_and_write.write_dill(self.MBD_system, filename)

    def showOpenFileDialog(self):
        """
        Open project file
        :return:
        """
        file_dialog = QtGui.QFileDialog()

        filename, file_type = file_dialog.getOpenFileNameAndFilter(self,
                                                                   caption='Open file',
                                                                   directory=QtCore.QString(self.MBD_folder_abs_path),
                                                                   filter=self._file_types)

        filename = str(filename)
        if filename:
            #   use dill to open binary coded object file
            if is_binary(filename):
                # self.MBD_system = None
                global MBD_system_new
                MBD_system_new = None

                MBD_system_new = read_and_write.read_dill(filename)

            #   use pickle to open ascii coded MBD object file
            else:
                pass

            # self.MBD_system.construct_MBD_system(MBD_file_abs_path=filename, _name=filename)
            self.project.removeChild(0)
            del self.MBD_system
            self.MBD_system = MBD_system_new
            self.MBD_system.dys = self
            self.MBD_system._parent = self.project

            self.project.addChild(self.MBD_system)
            # self.simulation_control_widget.MBD_system = MBD_system_new
            # self.simulation_control_widget.opengl_widget.MBD_system = MBD_system_new
            # self.project.addChild(self.MBD_system)

            self.tree_view_widget.setProjectItem(self.project, self.MBD_system)

        # print "self.MBD_system._name =", self.MBD_system._name
        # print "self.MBD_system =", self.MBD_system
        # self.simulation_control_widget.opengl_widget.context().reset()
        # self.simulation_control_widget.opengl_widget.close()
        # self.simulation_control_widget.opengl_widget = None
        # openql_widget = OpenGLWidget(MBD_system=copy.copy(self.MBD_system), parent=self)
        # self.setCentralWidget(openql_widget)
        # self.simulation_control_widget.opengl_widget =
        # print self.MBD_system
        # self.simulation_control_widget.setMBDsystem(self.MBD_system)
        # self.setCentralWidget(self.simulation_control_widget.opengl_widget)
        # sip.delete(self.simulation_control_widget.opengl_widget)
        # self.simulation_control_widget.opengl_widget = None
        self.simulation_control_widget.setMBDsystem(self.MBD_system)
        self.simulation_control_widget.opengl_widget.setMBDsystem(self.MBD_system)
        # self.simulation_control_widget.opengl_widget = OpenGLWidget(MBD_system=self.MBD_system, parent=self)
        # self.simulation_control_widget.opengl_widget.setMBDsystem(self.MBD_system)
        # self.simulation_control_widget.opengl_widget.geometry()
        # self.setCentralWidget(self.simulation_control_widget.opengl_widget)

        # except:
        #     pass
        # self.simulation_control_widget.opengl_widget.initializeGL()

        # self.simulation_control_widget.opengl_widget.geometry()

        # self.MBD_system.construct_MBD_system()
        # pprint(vars(self.MBD_system))
        # self.update_data()
#            pprint(vars(body))
#            pprint(vars(body.geom))
#        #    absolute file path of opened file
#        try:
#            self.MBD_file_abs_path = _open_file[0]
#            print "MBD_file_name =", self.MBD_file_abs_path
#        except:
#            None
#        if os.path.exists(self.MBD_file_abs_path):
#            project_folder_abs_path, MBD_filename = os.path.split(str(self.MBD_file_abs_path))
#
#            project_folder_name = os.path.basename(project_folder_abs_path)
#
        #

    def saveFile(self):
        """
        Save MBD system object to file
        """
        if not self.MBD_system.saved:
            #   file dialog
            file_dialog = QtGui.QFileDialog()
            file_dialog.setFileMode(QtGui.QFileDialog.ExistingFiles)

            self.MBD_file_abs_path = file_dialog.getSaveFileName(self,
                                                    caption=QtCore.QString("Save file"),
                                                    directory=QtCore.QString(self.MBD_folder_abs_path),
                                                    filter=QtCore.QString(self._file_types))

            self.MBD_file_abs_path = str(self.MBD_file_abs_path)

        if self.MBD_file_abs_path:
            with open(self.MBD_file_abs_path, "wb") as _file:
                save_data = self.MBD_system
                # pickle.dump(MBD_system, save, -1)
                dill.dump(save_data, _file)

                _file.close()
                self.MBD_system.saved = True

                logging.getLogger("DyS_logger").info("Project saved to file: %s. Size is %s", self.MBD_file_abs_path, convert_bytes_to_.convert_size(os.path.getsize(self.MBD_file_abs_path)))

    def show_control_panel(self):
        """
        Function shows control panel when hidden and selecte by user to be displayed
        :return:
        """
        self.simulation_control_widget.move(self.window_offset_x + self.simulation_control_widget.opengl_widget.initial_window_width + 10, 6)
        self.simulation_control_widget.show()
        
    def show_tree_view(self):
        self.tree_view_widget.move(self.window_offset_x - 228, 6)
        self.tree_view_widget.show()
    
    def show_output_widget(self):
        self.OutputWidget.show()

    def show_job_list(self):
        self.JobListWidget.move(self.window_offset_x + self.simulation_control_widget.opengl_widget.initial_window_width + 10 + 220, 6)
        self.JobListWidget.show()
        
    def show_python_console_widget(self):
        # self.python_console.move(0, 0)
        self.python_console.show()

    def take_snapshot(self, filename=None):
        """
        Function takes snapshot and opens save as dialog to save it as .png picture file.
        """
        captured_figure = self.simulation_control_widget.opengl_widget.takeSnapShot()

        _save_file = QtGui.QFileDialog()
        _save_file.setDirectory(self.MBD_system.MBD_folder_abs_path)
        abs_save_file_path, filetype = _save_file.getSaveFileNameAndFilter(self, "Save file", self.MBD_system.MBD_folder_abs_path, ("png (*.png)"))
        abs_save_file_path = str(abs_save_file_path).replace("/", "\\")
        captured_figure.save(abs_save_file_path)

    def create_status_bar(self):
        """
        Status bar widget       
        """
        self.infoText1 = QtGui.QLabel(self.statusBar())
        self.statusBar().addWidget(self.infoText1, stretch=2)
        self.infoText2 = QtGui.QLabel(self.statusBar())
        self.statusBar().addWidget(self.infoText2, stretch=1)
        self.simulation_status_text = QtGui.QLabel(self.statusBar())
        self.statusBar().addWidget(self.simulation_status_text, stretch=1)
  
        self.infoText1.setText("Energy (delta): ")
        self.infoText2.setText("Step No.: 0")
        self.simulation_status_text.setText("Status: Ready")
        self.repaint()

    def update_statusbar_text1(self, energy, energy_delta):
        self.infoText1.setText("Energy (delta): %4.3E"%energy+" (%4.3E"%energy_delta+")")

    def update_statusbar_text2(self, value):
        self.infoText2.setText("Step No.: " + str(value))

    @QtCore.pyqtSlot(str)
    def update_statusbar_simulation_status(self, simulation_status):
        """

        :param simulation_status:
        :return:
        """
        self.simulation_status_text.setText("Status: " + QtCore.QString(simulation_status))

    def about(self):
        QtGui.QMessageBox.about(self, self.tr("About Visualization Engine"), self.tr(
              "<b>Visualization Engine</b> version: %s <br><br> " 
              "<b>Motion functions and buttons</b>:<br>"
              "<b>Zoom</b>: Scroll with middle mouse button.<br>"
              "<b>Rotate</b>: Middle mouse button and mouse motion.<br>"
              "<b>Translate</b>: Middle mouse button + CTRL and mouse motion.<br><br>"
              "Visualization engine was designed for simulate dynamic motion of bodies from STL files." % (__version__)))

    def view(self):
        QtGui.QWidget()

    def _show(self):
        self.show()
        self.simulation_control_widget.show()
        self.tree_view_widget.show()

    def close(self):
        QtGui.qApp.quit()
    

if __name__ == '__main__':
    simulation_started = False
    app = QtGui.QApplication(sys.argv)  # 
    
    #    path to icon
    myappid = ":/application.png"
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

    dys = DySMainWindow()
    dys._show()

    # if win.simulation_control_widget.vtkWidget is not None:
    #     win.simulation_control_widget.vtkWidget.displayText()

    #***********************************************************
#     #    run simulation
#     if not simulation_started:
#         win.SimulationControlWidget.simulationStart()
#         win.SimulationControlWidget.solver.start_solver()
#         simulation_started = True
    #***********************************************************
    #   testing
    # print "testing"
    # q = dys.MBD_system.evaluate_q()
    # print dys.MBD_system.joints[1].spring.evaluate_Q_e(q)
    # print dys.MBD_system.springs
    # print "C(q, t) =", dys.simulation_control_widget.solver.analysis.DAE_fun.evaluate_C(dys.MBD_system.q0, 0.)
    # print "e =", dys.MBD_system.bodies[1].q
    # print "e_mesh =", dys.MBD_system.bodies[1].mesh.e
    # print "nodes =", dys.MBD_system.bodies[1].mesh.nodes
    # print "e_mesh =", dys.MBD_system.bodies[1].mesh.e
    # print "e ="

    # print dys.MBD_system.contacts[0]._contact_geometry_GCS(q)
    # print dys.MBD_system.contacts[0].pin_in_section_iP
    # print dys.MBD_system.contacts[0].pin_in_section_iR
    # print dys.MBD_system.contacts[0].pin_in_section_iPiR
    # print dys.MBD_system.contacts[0].pin_in_section
    # print dys.MBD_system.joints[1]
    # print dys.MBD_system.joints[1].evaluate_C(q)
    # dys.simulation_control_widget.solver.analysis.create_DAE()
    # print dys.simulation_control_widget.solver.analysis.DAE_fun.check_C(q, 0.)
    # pprint(vars(dys.MBD_system.joints[1]))
    # print dys.MBD_system.joints[1]
    # print "C() ="
    # print dys.MBD_system.joints[1].evaluate_C(q)
    # print "C_q(q) ="
    # print dys.MBD_system.joints[1].evaluate_C_q(q)
    # print "Cqi ="
    # print Cqi
    # print "Cqj ="
    # print Cqj
    # print dys.MBD_system.contacts[0]._contact_geometry_GCS(dys.MBD_system.q0)
    #
    # print "check_for_contact()"
    # print dys.MBD_system.contacts[0].check_for_contact(0, 0., dys.MBD_system.q0)
    # for force in dys.MBD_system.forces:
    #     print force._name, force._parent

    sys.exit(app.exec_())
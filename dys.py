"""
Created on 30. sep. 2013

@author: lskrinjar
"""
import ctypes
import logging
import logging.handlers
import os
import sys
import warnings

import dill
import numpy as np
from PyQt4 import QtCore, QtGui
from binaryornot.check import is_binary

from MBD_system import MBD_system_items
from MBD_system import convert_bytes_to_
from MBD_system import read_and_write
from MBD_system.MBD_system import MBDsystem
from job_list_widget.job_list_widget import JobListWidget
from preferences_widget.preferences_widget import PreferencesWidget
from simulation_control_widget.simulation_control_widget import SimulationControlWidget
from tree_view_widget.tree_view_widget import TreeViewWidget

np.set_printoptions(precision=6, threshold=None, edgeitems=100, linewidth=1000, suppress=False, nanstr=None, infstr=None)
# sys.setrecursionlimit(3000)
# print sys.getrecursionlimit()


try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s


# define authorship information
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


class MainWindow(QtGui.QMainWindow):
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
        super(MainWindow, self).__init__()
        
        #    minimum size
        self.setMinimumSize(500, 500)
        
        #    location of window on screen
        self.window_offset_x = 250
        self.window_offset_y = 50

        #    language settings
        self.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.UnitedStates))

        #    window name
        self.setWindowTitle("DyS")
        
        #    mouse tracking
        self.setMouseTracking(True)
        
        #    main window flags
        self.flags = QtCore.Qt.Window

        __main_dir = current_working_directory
        self._working_directory = QtCore.QString(os.path.abspath(__main_dir).replace("/", "\\"))

        #   get screen size
        screen = QtGui.QDesktopWidget().screenGeometry()
        #    move
        # self.move(self.window_offset_x, self.window_offset_y)
        dy = .05*screen.height()

        #   size of application
        #   position main widget on screen
        self.setGeometry(.25*screen.width(), dy, .5*screen.width(), .5*screen.height())

        #    MBD system
        MBD_folder_name_ = []
        self.MBD_file_abs_path = []

        # self.MBD_filename = "dys_0.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0"
        # self.MBD_filename = "dys_0_2_2.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_2_2"
#         self.MBD_filename = "dys_0_3_prismatic_joint.dprj"
#         MBD_folder_name_ = "dynamic_systems\\0_3_prismatic_joint"
#         self.MBD_filename = "dys_0_5.dprj"
#         MBD_folder_name_ = "dynamic_systems\\0_5"
#         self.MBD_filename = "dys_0_6_1.dprj"
#         MBD_folder_name_ = "dynamic_systems\\0_6_1"
#         self.MBD_filename = "dys_0_6_2.dprj"
#         MBD_folder_name_ = "dynamic_systems\\0_6_2"
#         self.MBD_filename = "dys_0_7_0_1.dprj"
#         MBD_folder_name_ = "dynamic_systems\\0_7_0_1"
#         self.MBD_filename = "dys_0_7_0_2.dprj"
#         MBD_folder_name_ = "dynamic_systems\\0_7_0_2"

#         self.MBD_filename = "0_7_3_0_contact_models.dprj"
#         MBD_folder_name_ = "dynamic_systems\\0_7_3_0_contact_models"

        # self.MBD_filename = "0_7_3_3_kinematic_analysis_1arm.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_3_3_kinematic_analysis_1arm"

        # self.MBD_filename = "0_7_3_3_kinematic_analysis_2arm.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_3_3_kinematic_analysis_2arm"

        # self.MBD_filename = "0_7_3_3_kinematic_analysis_slider_crank.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_3_3_kinematic_analysis_slider_crank"

#         self.MBD_filename = "0_7_3_3_kinematic_analysis_4bar.dprj"
#         MBD_folder_name_ = "dynamic_systems\\0_7_3_3_kinematic_analysis_4bar"

        # self.MBD_filename = "0_7_3_3_kinematic_analysis_efi.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_3_3_kinematic_analysis_efi"

        # self.MBD_filename = "0_7_3_3_dynamic_analysis_efi.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_3_3_dynamic_analysis_efi"
    
#         self.MBD_filename = "0_7_3_3_dynamic_analysis_efi_V2.dprj"
#         MBD_folder_name_ = "dynamic_systems\\0_7_3_3_dynamic_analysis_efi_V2"

#         self.MBD_filename = "dys_0_7_3_0.dprj"
#         MBD_folder_name_ = "dynamic_systems\\0_7_3_0"
        
        # self.MBD_filename = "0_7_3_0_plane_sphere.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_3_0_plane_sphere"

        # self.MBD_filename = "dys_0_7_2_revolute_clearence_joint.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_2_revolute_clearence_joint"

        # self.MBD_filename = "0_7_3_0_contact_models_cylinder.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_3_0_contact_models_cylinder"

        self.MBD_filename = "0_7_3_2_pin_slot_clearance_joint.dprj"
        MBD_folder_name_ = "dynamic_systems\\0_7_3_2_pin_slot_clearance_joint"

        # self.MBD_filename = "0_7_3_2_pin_slot_clearance_joint_paper.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_3_2_pin_slot_clearance_joint_paper"

        if MBD_folder_name_ == []:
            self.MBD_folder_abs_path = os.getcwd()
            project_filename = "Model_1"
        else:
            self._working_directory = os.path.join(os.getcwd(), "..")
            self.MBD_file_abs_path = os.path.abspath(os.path.join(self._working_directory, MBD_folder_name_, self.MBD_filename))

            self.MBD_folder_abs_path = os.path.join(self._working_directory, MBD_folder_name_)
            project_filename = os.path.basename(self.MBD_folder_abs_path)

        self._file_types = "Dynamic System Project File (*.dprj);;Multibody System File (*.mbd)"
        self.project = MBD_system_items.MBDsystemItem("self.project")
        
        
        self.MBD_system = MBDsystem(MBD_file_abs_path=self.MBD_file_abs_path, MBD_folder_name=MBD_folder_name_, MBD_folder_abs_path=self.MBD_folder_abs_path, MBD_filename=project_filename, parent=self.project)

        #   tree view widget
        self.TreeViewWidget = TreeViewWidget(self.project, self.MBD_system, parent=self)
        self.TreeViewWidget.setWindowFlags(self.flags)
        self.TreeViewWidget.show()
        #    move
        self.TreeViewWidget.move(.25*screen.width() - self.TreeViewWidget.frameGeometry().width(), dy)#self.frameGeometry().width(), self.frameGeometry().height())
        # self.TreeViewWidget.move(self.window_offset_x - self.TreeViewWidget.frameGeometry().width(), self.window_offset_y)

        #    simulation control widget
        self.simulation_control_widget = SimulationControlWidget(MBD_system=self.MBD_system, parent=self)
        self.simulation_control_widget.setWindowFlags(self.flags)
        self.setCentralWidget(self.simulation_control_widget.opengl_widget)
        self.simulation_control_widget.show()
        #    move
        self.simulation_control_widget.move(.25*screen.width() + self.geometry().width(), dy)

        #   preferences widget
        self.preferences_widget = PreferencesWidget(self.simulation_control_widget, parent=None)

        #    output widget
        #    maybe put in new thread? if output will be large
#         self.OutputWidget = OutputWidget(parent=self)
#         self.OutputWidget.setWindowFlags(self.flags)
#         self.OutputWidget.show()

        #    job list widget
        self.JobListWidget = JobListWidget(job_list_=["job1", "job2", "job3"], parent=self)
        self.JobListWidget.setWindowFlags(self.flags)

        #    resize
        self.resize(self.simulation_control_widget.opengl_widget.initial_window_width, self.simulation_control_widget.opengl_widget.initial_window_height)
        
        #    create actions
        self.create_actions()
        
        #    create menus
        self.create_menus()
        
        #    create status bar
        self.create_status_bar()

        #    signals and connections
        self.connect_signals()

    def connect_signals(self):
        """

        :return:
        """
        self.simulation_control_widget.solver.analysis.step_signal.signal_step.connect(self.update_statusbar_text2)

        self.simulation_control_widget.solver.running_signal.signal_running.connect(self.update_statusbar_text3)
        self.simulation_control_widget.solver.stopped_signal.signal_stopped.connect(self.update_statusbar_text3)
        self.simulation_control_widget.solver.analysis.finished_signal.signal_finished.connect(self.update_statusbar_text3)


        # self.simulation_control_widget.solver.solveODE.finished_signal.signal_finished.connect(self.update_statusbar_text3)

        # self.simulation_control_widget.solver.analysis.filename_signal.signal_filename.connect(self.TreeViewWidget.add_solution_data)
        self.simulation_control_widget.step_num_signal.signal_step.connect(self.update_statusbar_text2)
        self.simulation_control_widget.status_signal.signal_status.connect(self.update_statusbar_text3)
        
        
        self.TreeViewWidget.filename_signal.signal_filename.connect(self.simulation_control_widget.load_solution_file)
        # self.TreeViewWidget._widget.repaintGL_signal.signal_repaintGL.connect(self.simulation_control_widget.opengl_widget._repaintGL)
        
        
        self.simulation_control_widget.solver.analysis.energy_signal.signal_energy.connect(self.update_statusbar_text1)
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
        self.openAction.triggered.connect(self.openFile)
        #    close
        self.closeAction = QtGui.QAction('Close Project', self, statusTip='Close Project')
        self.closeAction.triggered.connect(self.closeProjectFile)
        #   save
        self.saveAction = QtGui.QAction('Save', self)
        self.saveAction.setStatusTip("Save file")
        self.saveAction.setShortcut("Ctrl+S")
        self.saveAction.triggered.connect(self.saveFile)

        #    exit
        self.exitAction = QtGui.QAction('Exit', self, shortcut='Ctrl+Q', statusTip='Exit Application')
        self.exitAction.triggered.connect(self.close)

        #    view - menu
        #    view
        self.viewActionFront = QtGui.QAction('Front', self, shortcut='Ctrl+F', statusTip='Front view')
        self.viewActionFront.triggered.connect(self.simulation_control_widget.opengl_widget.viewFront)
        
        self.viewActionBack = QtGui.QAction('Back', self, shortcut='Ctrl+B', statusTip='Back view')
        self.viewActionBack.triggered.connect(self.simulation_control_widget.opengl_widget.viewBack)
        
        self.viewActionBottom = QtGui.QAction('Bottom', self, statusTip='Bottom view')
        self.viewActionBottom.triggered.connect(self.simulation_control_widget.opengl_widget.viewBottom)
        
        self.viewActionTop = QtGui.QAction('Top', self, statusTip='Top view')
        self.viewActionTop.triggered.connect(self.simulation_control_widget.opengl_widget.viewTop)
        
        self.viewActionLeft = QtGui.QAction('Left', self, shortcut='Ctrl+L', statusTip='Left view')
        self.viewActionLeft.triggered.connect(self.simulation_control_widget.opengl_widget.viewLeft)
        
        self.viewActionRight = QtGui.QAction('Right', self, shortcut='Ctrl+R', statusTip='Right view')
        self.viewActionRight.triggered.connect(self.simulation_control_widget.opengl_widget.viewRight)
        
        self.viewActionIsometric = QtGui.QAction('Isometric', self, shortcut='Ctrl+I', statusTip='Isometric view')
        self.viewActionIsometric.triggered.connect(self.simulation_control_widget.opengl_widget.viewIsometric)
        
        
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
            if event.key() == QtCore.Qt.Key_Right:
                self.simulation_control_widget.animation_forward()
            
            if event.key() == QtCore.Qt.Key_Left:
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
        openFile = QtGui.QAction('Open', self)
        openFile.setShortcut('Ctrl+O')
        openFile.setStatusTip('Open Project')
        openFile.triggered.connect(self.openFile)
        self.fileMenu.addAction(openFile)
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
        self.viewMenu.addAction(self.viewActionFront)
        self.viewMenu.addAction(self.viewActionBack)
        self.viewMenu.addAction(self.viewActionBottom)
        self.viewMenu.addAction(self.viewActionTop)
        self.viewMenu.addAction(self.viewActionLeft)
        self.viewMenu.addAction(self.viewActionRight)
        self.viewMenu.addSeparator()
        self.viewMenu.addAction(self.viewActionIsometric)
        self.viewMenu.addSeparator()
        self.viewMenu.addAction(self.take_snapshot_action)
        
        #    window
        self.windowMenu = QtGui.QMenu(self.tr("&Window"), self)
        self.windowMenu.addAction(self.show_control_panel_action)
        self.windowMenu.addAction(self.show_tree_view_action)
        self.windowMenu.addAction(self.show_output_widget_action)
        self.windowMenu.addAction(self.show_job_list_action)
        
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

    def newFile(self):
        a = 1

    def update_data(self):
        """
        
        """
        self.simulation_control_widget.opengl_widget.update_data(self.MBD_system.bodies)
        self.simulation_control_widget.opengl_widget.repaintGL()

    def showSaveAsFileDialog(self):
        """
        
        """
        saveAs_dialog = QtGui.QFileDialog
        filename, file_type = saveAs_dialog.getSaveFileNameAndFilter(self, 'Save File', self.current_working_directory, ("DyS Project Files (*.dprj)"))
        
        read_and_write.write(self.MBD_system, filename)
#        f = open(filename, 'w')
#        filedata = self.text.toPlainText()
#        f.write(filedata)
#        f.close() 

    def createNewProject(self):
        """
        
        """

    def closeProjectFile(self):
        """
        
        """
        self.MBD_system.delete_MBD_system()
        self.TreeViewWidget.repaint()

    def openFile(self):
        """
        Open file
        :return:
        """
        file_dialog = QtGui.QFileDialog()

        filename, file_type = file_dialog.getOpenFileNameAndFilter(self,
                                                                   caption='Open file',
                                                                   directory=QtCore.QString(self.MBD_folder_abs_path),
                                                                   filter=self._file_types)

        filename = str(filename)
        if filename:
            if is_binary(filename):
                with open(filename, 'rb') as _file:
                    print "_file =", _file
                    self.MBD_system = dill.load(_file)
            else:
                pass
        # self.MBD_system = read_and_write.read(filename)
        #
        # self.MBD_system.construct_MBD_system()
        #
        # self.update_data()
#            pprint(vars(body))
#            pprint(vars(body.geom))
#        #    absolute file path of opened file
#        try:
#            self.MBD_file_abs_path = _open_file[0]
#            print "MBD_file_name =", self.MBD_file_abs_path
#        except:
#            None
#
#
#        if os.path.exists(self.MBD_file_abs_path):
#            project_folder_abs_path, MBD_filename = os.path.split(str(self.MBD_file_abs_path))
#
#            project_folder_name = os.path.basename(project_folder_abs_path)
#
        # self.MBD_system.construct_MBD_system(MBD_file_abs_path=filename, _name=filename)

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
        self.TreeViewWidget.move(self.window_offset_x - 228, 6)
        self.TreeViewWidget.show()
    
    def show_output_widget(self):
        self.OutputWidget.show()

    def show_job_list(self):
        self.JobListWidget.move(self.window_offset_x + self.simulation_control_widget.opengl_widget.initial_window_width + 10 + 220, 6)
        self.JobListWidget.show()

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
        self.infoText3 = QtGui.QLabel(self.statusBar())
        self.statusBar().addWidget(self.infoText3, stretch=1)
  
        self.infoText1.setText("Energy (delta): ")
        self.infoText2.setText("Step No.: " + str(self.simulation_control_widget.solver.analysis.step))
        self.infoText3.setText("Status: Ready")
        self.repaint()

    def update_statusbar_text1(self, energy, energy_delta):
        self.infoText1.setText("Energy (delta): %4.3E"%energy+" (%4.3E"%energy_delta+")")

    def update_statusbar_text2(self, value):
        self.infoText2.setText("Step No.: " + str(value))

    def update_statusbar_text3(self, simulation_status_string):
        # print "simulation_status_string =", simulation_status_string
        if self.simulation_control_widget.solver.analysis.running == True:
            simulation_status_string = "Running"

        elif self.simulation_control_widget.solver.analysis.stopped == True:
            simulation_status_string = "Stopped"

        elif self.simulation_control_widget.solver.analysis.finished == True:
            simulation_status_string = "Finished"

        else:
            simulation_status_string = simulation_status_string
        # print "animation"
        self.simulation_control_widget.solver.analysis.finished_signal
        simulation_status_string = simulation_status_string

        self.infoText3.setText("Status: " + QtCore.QString(simulation_status_string))

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

    def close(self):
        QtGui.qApp.quit()
    

if __name__ == '__main__':
    simulation_started = False
#     pyqtRemoveInputHook()
    app = QtGui.QApplication(sys.argv)  # 
    app.setStyle("windows")
#     plastique
#     cde
#     motif
#     sgi
#     windows
#     cleanlooks
#     mac
    app.setWindowIcon(QtGui.QIcon(":/application.png"))
    QtGui.QSystemTrayIcon(QtGui.QIcon(":/application.png"))
    
    #    path to icon
    myappid = ":/application.png"
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

    win = MainWindow()
    win.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
    win.show()

    #***********************************************************
#     #    run simulation
#     if not simulation_started:
#         win.SimulationControlWidget.simulationStart()
#         win.SimulationControlWidget.solver.start_solver()
#         simulation_started = True
    #***********************************************************
    sys.exit(app.exec_())

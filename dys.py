"""
Created on 30. sep. 2013

@author: lskrinjar
"""
import warnings
import os
import sys
import inspect
import time
import numpy as np
import ctypes
from pprint import pprint
from PyQt4 import QtCore, QtGui
import subprocess


from MBD_system.MBD_system import MBDsystem
from MBD_system import MBD_system_items
from job_list_widget.job_list_widget import JobListWidget
from output_widget.output_widget import OutputWidget
from simulation_control_widget.simulation_control_widget import SimulationControlWidget
from tree_view_widget.tree_view_widget import TreeViewWidget
from MBD_system import read_and_write


import qrc_resources


np.set_printoptions(precision=4, threshold=None, edgeitems=100, linewidth=1000, suppress=False, nanstr=None, infstr=None)
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
        dy = .1*screen.height()

        #   size of application
        #   position main widget on screen
        self.setGeometry(.25*screen.width(), dy, .5*screen.width(), .5*screen.height())

        #    MBD system
        MBD_folder_name_ = []
        self.MBD_file_abs_path_ = []


        # self.MBD_filename = "dys_0.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0"
        # self.MBD_filename = "dys_0_2_2.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_2_2"
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
#         self.MBD_filename = "dys_0_7_2_revolute_clearence_joint.dprj"
#         MBD_folder_name_ = "dynamic_systems\\0_7_2_revolute_clearence_joint"
#         self.MBD_filename = "0_7_3_0_contact_models.dprj"
#         MBD_folder_name_ = "dynamic_systems\\0_7_3_0_contact_models"
        self.MBD_filename = "0_7_3_0_contact_models_cylinder.dprj"
        MBD_folder_name_ = "dynamic_systems\\0_7_3_0_contact_models_cylinder"


#         self.MBD_filename = "dys_0_7_3_0.dprj"
#         MBD_folder_name_ = "dynamic_systems\\0_7_3_0"
        
        # self.MBD_filename = "0_7_3_0_plane_sphere.dprj"
        # MBD_folder_name_ = "dynamic_systems\\0_7_3_0_plane_sphere"

        if MBD_folder_name_ == []:
            MBD_folder_abs_path_ = os.getcwd()
            project_filename = "Model_1"
        else:
            self._working_directory = os.path.join(os.getcwd(), "..")
            self.MBD_file_abs_path_ = os.path.abspath(os.path.join(self._working_directory, MBD_folder_name_, self.MBD_filename))

            MBD_folder_abs_path_ = os.path.join(self._working_directory, MBD_folder_name_)
            project_filename = os.path.basename(MBD_folder_abs_path_)

        projectNode = MBD_system_items.MBDsystemItem("projectNode")


        self.MBD_system = MBDsystem(MBD_file_abs_path=self.MBD_file_abs_path_, MBD_folder_name=MBD_folder_name_, MBD_folder_abs_path=MBD_folder_abs_path_, MBD_filename=project_filename, parent=projectNode)

        #   tree view widget
        self.TreeViewWidget = TreeViewWidget(projectNode, parent=self)
        self.TreeViewWidget.setWindowFlags(self.flags)
        self.TreeViewWidget.show()
        #    move
        self.TreeViewWidget.move(.25*screen.width() - self.TreeViewWidget.frameGeometry().width(), dy)#self.frameGeometry().width(), self.frameGeometry().height())
        # self.TreeViewWidget.move(self.window_offset_x - self.TreeViewWidget.frameGeometry().width(), self.window_offset_y)

        #    simulation control widget
        self.SimulationControlWidget = SimulationControlWidget(MBD_system=self.MBD_system, parent=self)
        self.SimulationControlWidget.setWindowFlags(self.flags)
        self.setCentralWidget(self.SimulationControlWidget.OpenGLWidget)
        self.SimulationControlWidget.show()
        #    move
        self.SimulationControlWidget.move(.25*screen.width() + self.geometry().width(), dy)

        #    output widget
        #    maybe put in new thread? if output will be large
#         self.OutputWidget = OutputWidget(parent=self)
#         self.OutputWidget.setWindowFlags(self.flags)
#         self.OutputWidget.show()

        #    job list widget
        self.JobListWidget = JobListWidget(job_list_=["job1", "job2", "job3"], parent=self)
        self.JobListWidget.setWindowFlags(self.flags)

        #    resize
        self.resize(self.SimulationControlWidget.OpenGLWidget.initial_window_width, self.SimulationControlWidget.OpenGLWidget.initial_window_height)
        
        #    create actions
        self.create_actions()
        
        #    create menus
        self.create_menus()
        
        #    create status bar
        self.create_status_bar()
        
        #    signals and connections
        self.SimulationControlWidget.solver.solveODE.step_signal.signal_step.connect(self.update_statusbar_text2)

        self.SimulationControlWidget.solver.running_signal.signal_running.connect(self.update_statusbar_text3)
        self.SimulationControlWidget.solver.stopped_signal.signal_stopped.connect(self.update_statusbar_text3)
        self.SimulationControlWidget.solver.solveODE.finished_signal.signal_finished.connect(self.update_statusbar_text3)


        # self.SimulationControlWidget.solver.solveODE.finished_signal.signal_finished.connect(self.update_statusbar_text3)

        self.SimulationControlWidget.solver.solveODE.filename_signal.signal_filename.connect(self.TreeViewWidget.add_solution_data)
        self.SimulationControlWidget.step_num_signal.signal_step.connect(self.update_statusbar_text2)
        self.SimulationControlWidget.status_signal.signal_status.connect(self.update_statusbar_text3)
        
        
        self.TreeViewWidget.filename_signal.signal_filename.connect(self.SimulationControlWidget.load_solution_file)
        # self.TreeViewWidget._widget.repaintGL_signal.signal_repaintGL.connect(self.SimulationControlWidget.OpenGLWidget._repaintGL)
        
        
        self.SimulationControlWidget.solver.solveODE.energy_signal.signal_energy.connect(self.update_statusbar_text1)
        self.SimulationControlWidget.energy_signal.signal_energy.connect(self.update_statusbar_text1)

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
        
        #    exit
        self.exitAction = QtGui.QAction('Exit', self, shortcut='Ctrl+Q', statusTip='Exit Application')
        self.exitAction.triggered.connect(self.close)
        
        
        #    view - menu
        #    view
        self.viewActionFront = QtGui.QAction('Front', self, shortcut='Ctrl+F', statusTip='Front view')
        self.viewActionFront.triggered.connect(self.SimulationControlWidget.OpenGLWidget.viewFront)
        
        self.viewActionBack = QtGui.QAction('Back', self, shortcut='Ctrl+B', statusTip='Back view')
        self.viewActionBack.triggered.connect(self.SimulationControlWidget.OpenGLWidget.viewBack)
        
        self.viewActionBottom = QtGui.QAction('Bottom', self, statusTip='Bottom view')
        self.viewActionBottom.triggered.connect(self.SimulationControlWidget.OpenGLWidget.viewBottom)
        
        self.viewActionTop = QtGui.QAction('Top', self, statusTip='Top view')
        self.viewActionTop.triggered.connect(self.SimulationControlWidget.OpenGLWidget.viewTop)
        
        self.viewActionLeft = QtGui.QAction('Left', self, shortcut='Ctrl+L', statusTip='Left view')
        self.viewActionLeft.triggered.connect(self.SimulationControlWidget.OpenGLWidget.viewLeft)
        
        self.viewActionRight = QtGui.QAction('Right', self, shortcut='Ctrl+R', statusTip='Right view')
        self.viewActionRight.triggered.connect(self.SimulationControlWidget.OpenGLWidget.viewRight)
        
        self.viewActionIsometric = QtGui.QAction('Isometric', self, shortcut='Ctrl+I', statusTip='Isometric view')
        self.viewActionIsometric.triggered.connect(self.SimulationControlWidget.OpenGLWidget.viewIsometric)
        
        
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
        self.show_change_background_color_dialog_action = QtGui.QAction('Background color',
                                                                   self,
                                                                   shortcut="None",
                                                                   statusTip='Status tip')
        self.show_change_background_color_dialog_action.triggered.connect(self.show_change_background_color_dialog)
        
        
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
                self.SimulationControlWidget.animation_forward()
            
            if event.key() == QtCore.Qt.Key_Left:
                self.SimulationControlWidget.animation_backward()

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
        openFile.triggered.connect(self.showOpenFileDialog)
        self.fileMenu.addAction(openFile)
        self.fileMenu.addSeparator()

        #    save file
        saveFile = QtGui.QAction('Save', self)
        saveFile.setShortcut('Ctrl+S')
        saveFile.setStatusTip('Save Project')
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
        self.settingsMenu.addAction(self.show_change_background_color_dialog_action)
        
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
        self.SimulationControlWidget.OpenGLWidget.update_data(self.MBD_system.bodies)
        self.SimulationControlWidget.OpenGLWidget.repaintGL()

    def showOpenFileDialog(self):

        _open_file_ = QtGui.QFileDialog()
        _open_file_.setDirectory(self._working_directory)
        
        filename, file_type = _open_file_.getOpenFileNameAndFilter(self, 'Open file', self.current_working_directory, ("DyS Project Files (*.dprj)"))
        
        self.MBD_system = read_and_write.read(filename)
        
        self.MBD_system.construct_MBD_system()
        
        self.update_data()
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
        self.MBD_system.construct_MBD_system(MBD_file_abs_path=filename, _name=filename)

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

    def saveFile(self):
        """
        
        """
        filename = QtGui.QFileDialog.getSaveFileNameAndFilter(self)

    def show_change_background_color_dialog(self):
        color = QtGui.QColorDialog.getColor()
        if color.isValid():
            self.SimulationControlWidget.OpenGLWidget.qglClearColor(color)
            self.SimulationControlWidget.OpenGLWidget.updateGL()

    def show_control_panel(self):
        """
        Function shows control panel when hidden and selecte by user to be displayed
        :return:
        """
        self.SimulationControlWidget.move(self.window_offset_x + self.SimulationControlWidget.OpenGLWidget.initial_window_width + 10, 6)
        self.SimulationControlWidget.show()
        
    def show_tree_view(self):
        self.TreeViewWidget.move(self.window_offset_x - 228, 6)
        self.TreeViewWidget.show()
    
    def show_output_widget(self):
        self.OutputWidget.show()

    def show_job_list(self):
        self.JobListWidget.move(self.window_offset_x + self.SimulationControlWidget.OpenGLWidget.initial_window_width + 10 + 220, 6)
        self.JobListWidget.show()

    def take_snapshot(self, filename=[]):
        """
        Function takes snapshot and opens save as dialog to save it as .png picture file.
        """
        captured_figure = self.SimulationControlWidget.OpenGLWidget.takeSnapShot()

        _save_file = QtGui.QFileDialog()
        _save_file.setDirectory(self.MBD_system.MBD_folder_abs_path_)
        abs_save_file_path = _save_file.getSaveFileNameAndFilter(self, "Save file", self.MBD_system.MBD_folder_abs_path_, ("png (*.png)"))
        abs_save_file_path.replace("/", "\\")
        captured_figure.save(abs_save_file_path, 'png')

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
        self.infoText2.setText("Step No.: " + str(self.SimulationControlWidget.solver.solveODE.step))
        self.infoText3.setText("Status: Ready")
        self.repaint()

    def update_statusbar_text1(self, energy, energy_delta):
        self.infoText1.setText("Energy (delta): %4.3E"%energy+" (%4.3E"%energy_delta+")")

    def update_statusbar_text2(self, value):
        self.infoText2.setText("Step No.: " + str(value))

    def update_statusbar_text3(self, simulation_status_string):
        # print "simulation_status_string =", simulation_status_string
        if self.SimulationControlWidget.solver.solveODE.running == True:
            simulation_status_string = "Running"

        elif self.SimulationControlWidget.solver.solveODE.stopped == True:
            simulation_status_string = "Stopped"

        elif self.SimulationControlWidget.solver.solveODE.finished == True:
            simulation_status_string = "Finished"

        else:
            simulation_status_string = simulation_status_string
        # print "animation"
        self.SimulationControlWidget.solver.solveODE.finished_signal
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

"""
Created on 3. mar. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
"""
import os
import datetime
from pprint import pprint
import sys
import time
import gc
import shutil
import cProfile

import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from moviepy.editor import ImageSequenceClip



from opengl_widget.opengl_widget import OpenGLWidget
from simulation_control_widget_ui import Ui_Form
from solver.solver import Solver
from signals import StatusSignal
from MBD_system.check_filename import check_filename


try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s


class stepSignal(QtCore.QObject):
    signal_step = QtCore.pyqtSignal(int, name='signal_step')
    
    
class EnergySignal(QtCore.QObject):
    signal_energy = QtCore.pyqtSignal(float, float, name='energy')


class SimulationControlWidget(QtGui.QWidget):
    """
    Control panel GUI
    """
    def __init__(self, MBD_system, parent=None, flags=0):
        """
        Constructor
        """
        super(SimulationControlWidget, self).__init__(parent)
        self._parent = parent
        self.ui = Ui_Form()
        self.ui.setupUi(self)
        self.ui.simulationStopButton.setEnabled(False)
        self.ui.simulationResetButton.setEnabled(False)
        
        self.ui.forwardButton.setEnabled(False)
        self.ui.backwardButton.setEnabled(False)
        self.ui.playButton.setEnabled(False)
        
        self.MBD_system = MBD_system

        if self.MBD_system.loadSolutionFileWhenFinished:
            self.ui.loadSolutionFileStatus.setChecked(True)
        else:
            self.ui.loadSolutionFileStatus.setChecked(False)

        #    sets opengl widget in central widget position
        self.OpenGLWidget = OpenGLWidget(MBD_system=MBD_system, parent=self._parent)
        
        self.__status = "simulation"  # simulation or animation

        #   set integration method to display
        try:
            _index = self.ui.integrationMethodComboBox.findText(QtCore.QString(self.MBD_system.integrationMethod.title()))
        except:
            _index = 0

        if _index != -1:
            self.ui.integrationMethodComboBox.setCurrentIndex(_index)
        else:
            self.ui.integrationMethodComboBox.setCurrentIndex(_index)
        
        #    signals
        self.step_num_signal = stepSignal()
        self.energy_signal = EnergySignal()
        self.status_signal = StatusSignal()
        
#         self.graphWidget = None
#         self.graphWidget = GraphWidget(MBD_system=MBD_system_, parent=parent)
#         self.graphWidget.setWindowFlags(parent.windowFlags())
#         self.graphWidget.show()

        #    set validators to limit user input
        __validator_dbl = QtGui.QDoubleValidator()
        __validator_int = QtGui.QIntValidator()
        
        #    predefined values
        #    end time
        self.ui.endTime.setText(str(self.MBD_system.t_n))
        self.ui.endTime.setValidator(__validator_dbl)

        # self.ui.simulation_progressBar.setMinimum(0)
        # self.ui.simulation_progressBar.setMaximum(int(1/self.MBD_system.t_n))
        # self.ui.simulation_progressBar.setValue(0)
        
        #    Hmax
        if self.MBD_system.t_n/100 < self.MBD_system.Hmax:
            self.MBD_system.Hmax = 0.01*self.MBD_system.t_n

        #   Hmin
        self.MBD_system.evaluate_Hmin()

        self.ui.Hmax.setText(str(self.MBD_system.Hmax))
        self.ui.Hmax.setValidator(__validator_dbl)
        
        #    Hmin
        self.ui.Hmin.setText(str(self.MBD_system.Hmin))
        self.ui.Hmin.setValidator(__validator_dbl)
        #    abs tol
        self.ui.absTol.setText(str(self.MBD_system.absTol))
        self.ui.absTol.setValidator(__validator_dbl)
        #    rel tol
        self.ui.relTol.setText(str(self.MBD_system.relTol))
        self.ui.relTol.setValidator(__validator_dbl)

        
        #    create solver thread
        self.solver = Solver(MBD_system=MBD_system, parent=self)#parent
        self._solver_thread = QThread()
        self._solver_thread.start()
        self.solver.moveToThread(self._solver_thread)


        #    update display on every i-th simulation step
        self.ui.updateDisplayStep.setText(str(int(self.solver.solveODE.update_opengl_widget_every_Nth_step)))
        self.ui.updateDisplayStep.setValidator(__validator_int)
        self._delta_step = self.solver.solveODE.update_opengl_widget_every_Nth_step
        
        
        self.ui.currentStep.setEnabled(False)
        self.ui.currentStep.setValidator(__validator_int)

        #    connections and signals
        self.ui.simulationStartButton.clicked.connect(self.simulationStart)
        self.ui.simulationStartButton.clicked.connect(self.solver.start_solver)
        
        self.ui.simulationStopButton.clicked.connect(self.simulationStop)
        self.ui.simulationStopButton.clicked.connect(self.solver.stop_solver)

        
        self.solver.solveODE.finished_signal.signal_finished.connect(self.simulationFinished)
        self.solver.solveODE.filename_signal.signal_filename.connect(self.__automaticaly_load_solution_file)

        self.solver.solveODE.solution_signal.solution_data.connect(self.__automaticaly_load_solution_file)
        self.solver.solveODE.solution_signal.solution_data.connect(self._parent.TreeViewWidget.add_solution_data)


        self.ui.simulationResetButton.clicked.connect(self.simulationReset)
        self.ui.endTime.textChanged.connect(self.__update_endTime)
        self.ui.Hmax.textChanged.connect(self.__update_Hmax)
        self.ui.Hmin.textChanged.connect(self.__update_Hmin)
        
        self.ui.loadSolutionFileStatus.stateChanged.connect(self.__update_loadSolutionFileWhenFinished)
        
        self.ui.currentStep.textChanged.connect(self.__update_currentStep)
        self.ui.updateDisplayStep.textChanged.connect(self.__update_updateDisplayStep)


        self.ui.backwardButton.clicked.connect(self.animation_backward) #clicked
        self.ui.forwardButton.clicked.connect(self.animation_forward)
        self.ui.playButton.clicked.connect(self.animationPlay)
        

        #    signal repaintGL.signal_repaintGL from self.solver triggers self.OpenGLWidget.repaintGL
        self.solver.solveODE.repaintGL_signal.signal_repaintGL.connect(self.OpenGLWidget.repaintGL)

        #   signal for take a snapshot
        self.solver.solveODE.save_screenshot_signal.signal_saveScreenshot.connect(self.take_snapshot)

        #   signal time integration error
        self.solver.solveODE.error_time_integration_signal.signal_time_integration_error.connect(self._time_integration_error)

        #   change integration method
        self.ui.integrationMethodComboBox.currentIndexChanged.connect(self.selectedIntegrationMethod)

        self._parent.TreeViewWidget.create_animation_file.signal_createAnimationFile.connect(self._create_animation_file)

        #   profiler
        self.profile = cProfile.Profile()

    @pyqtSlot()
    def profile_functions(self):
        """
        Function profiles main functions that are used in numerical integration
        :return:
        """
        print "profile here"
        self.profile.enable()
        self.solver.solveODE.ode_fun.create_M()
        self.profile.disable()
        self.profile.print_stats()

    def _time_integration_error(self):
        """

        :return:
        """
        QtGui.QMessageBox.critical(self._parent, "Error!", "Hmin exceeded! Procedure failed!",QtGui.QMessageBox.Ok, QtGui.QMessageBox.NoButton,QtGui.QMessageBox.NoButton)

    def __update_loadSolutionFileWhenFinished(self):
        """
        
        """

        if self.ui.loadSolutionFileStatus.isChecked():
            self.MBD_system.loadSolutionFileWhenFinished = True
        else:
            self.MBD_system.loadSolutionFileWhenFinished = False

    def __update_currentStep(self):
        """
        
        """
        try:
            self.__step = int(self.ui.currentStep.text())
            self.MBD_system.update_coordinates_and_angles_of_all_bodies(self.q[self.__step, :])
            self.MBD_system.update_simulation_properties(time=self.t[self.__step], step_num=self.__step)
            self.step_num_signal.signal_step.emit(self.__step)
            self.__update_GL()
        except:
            pass

    def __automaticaly_load_solution_file(self, solution_data_object_ID=None, filename=None):
        """
        Function loads solution data from input:
        if input is filename solution data is read from specified file
        if input is solution data object, solution is read from object attribute solution_data
        :param filename:                a filename of the solution data file
        :param solution_data_object:    a solution data object with attribute solution_data
        :return: None
        """
        if self.MBD_system.loadSolutionFileWhenFinished:
            self.load_solution_file(solution_data_object_ID, filename)

    def load_solution_file(self, solution_object_id=None, filename=None):
        """
        Function loads solution data from object (first) or from file (second)
        """
        if solution_object_id is not None:
            #   assign solution data from solution data object to ne variable
            solution_data = self.solver.solveODE._solution_data.load_solution_data()
        elif filename is not None:
            pass
        else:
            print "Solution data not loaded!"

        # if filename is not None and solution_object_id is None:
        #     solution_data = self.solver.solveODE.load_simulation_solution_from_file(filename)

        #   assign a solution data object to pointer of object attribute
        for sol in self.MBD_system.solutions:
            if id(sol) == solution_object_id:
                self.MBD_system.loaded_solution = sol

        self.step = solution_data[:, 0]
        self.energy = solution_data[:, 1]
        self.error = solution_data[:, 2]
        self.dt = solution_data[:, 3]
        self.t = solution_data[:, 4]
        self.q = solution_data[:, 5:]

        self.__step = 0
        
        self.__status = "animation"
        
        self.ui.forwardButton.setEnabled(True)
        self.ui.backwardButton.setEnabled(False)
        self.ui.playButton.setEnabled(True)
        self.ui.currentStep.setEnabled(True)
        
        self.ui.currentStep.setText(str(int(self.__step)))
        
        self.ui.solutionFileLoaded_display.setText(filename)
        self.ui.numberOfSteps.setText(str(int(len(self.step)-1)))
        self.__update_GL()
        
        # self.status_signal.signal_status.emit("Animation")

        self.simulationReset()

    def __update_GL(self):
        """
        
        """
        self.ui.currentStep.setText(str(int(self.__step)))
        self.MBD_system.update_coordinates_and_angles_of_all_bodies(self.q[int(self.__step), :])
        self.OpenGLWidget.repaintGL(int(self.__step))
        
        #    energy data signal
        _energy = self.energy[self.__step]
        _energy_delta = _energy - self.energy[int(self.__step)-1]
        self.energy_signal.signal_energy.emit(_energy, _energy_delta)

    def animation_forward(self):
        """
        Go to next solution time step
        """
        if (self.__step + self._delta_step) <= self.step[-1]:
            self.__step += int(self._delta_step)
            self.__update_GL()
            self.ui.backwardButton.setEnabled(True)
        else:
            self.ui.forwardButton.setEnabled(False)

    def animation_backward(self):
        """
        Return to prevouos solution time step
        """
        if (self.__step - self._delta_step) >= self.step[0]:
            self.__step -= int(self._delta_step)
            self.__update_GL()
            self.ui.forwardButton.setEnabled(True)
        else:
            self.ui.backwardButton.setEnabled(False)

    def take_snapshot(self):
        """
        Create a snapshot
        """
        captured_figure = self.OpenGLWidget.takeSnapShot()
        captured_figure.save(self.solver.solveODE.screenshot_filename_abs_path + '.png', 'png')

    def selectedIntegrationMethod(self, int):
        """
        Assign a selected integration method to object attribute
        """
        self.MBD_system.integrationMethod = str(self.ui.integrationMethodComboBox.currentText())
    
    def __update_updateDisplayStep(self):
        """

        :return:
        """
        self.solver.solveODE.update_opengl_widget_every_Nth_step = int(self.ui.updateDisplayStep.text())
        self._delta_step = int(self.ui.updateDisplayStep.text())
        self.MBD_system.updateEveryIthStep = self._delta_step

    def __update_Hmax(self):
        try:
            self.Hmax = float(self.ui.Hmax.text()) 
            self.MBD_system.Hmax = self.Hmax
#             self.solver.update_simulation_control_parameters(dt_=self.stepSize)
        except:
            None

    def __update_Hmin(self):
        try:
            self.Hmin = float(self.ui.Hmin.text())
            self.MBD_system.Hmin = self.Hmin
            self.solver.update_simulation_control_parameters(dt_=self.stepSize)
        except:
            None

    def __update_endTime(self):
        try:
            self.endTime = float(self.ui.endTime.text())
            self.MBD_system.t_n = self.endTime
        except:
#             QtGui.QMessageBox.about(self, 'Error', 'Input can only be a number')
            pass

    def setWindowFlags(self, flags):
        super(SimulationControlWidget, self).setWindowFlags(flags)

    def simulationStart(self):
        if not self.solver.solveODE.running:
            self.job_info_time_started = time.clock()
            self.solver.solveODE.running = True
            self.solver.solveODE.stopped = False
            self.solver.solveODE.finished = False

            self.ui.simulationStartButton.setEnabled(False)
            self.ui.simulationResetButton.setEnabled(False)
            self.ui.simulationStopButton.setEnabled(True)
            
    def simulationStop(self):
        if self.solver.solveODE.running:
            self.solver.solveODE.stopped = True
            
            self.ui.simulationStartButton.setEnabled(True)
            self.ui.simulationStopButton.setEnabled(False)
            self.ui.simulationResetButton.setEnabled(True)

    def simulationFinished(self):
        self.ui.simulationStopButton.setEnabled(False)
        self.ui.simulationResetButton.setEnabled(True)
        self.job_info_time_finished = time.clock()

    def simulationReset(self):
        self.solver.solveODE.restore_initial_condition()
        # self.ui.simulationResetButton.setEnabled(False)
        self.ui.simulationStartButton.setEnabled(True)

    def animationPlay(self):
        """
        Function plays animation when solution data is loaded
        :return:
        """
        for _step in xrange(0, len(self.step), int(self._delta_step)):
            self.__step = int(_step)
            self.__update_GL()

            time.sleep(self.ui.playBackSpeed_doubleSpinBox.value()*1E-2)

    def _create_animation_file(self):
        """
        Function creates video file of animation
        :return:
        """
        #   check if folder exists and if not create it
        folder = "_animation"
        if not os.path.exists(folder):
            os.makedirs(folder)

        #   absolute path to store animation images to use them to create video file
        animation_folder_abs_path = os.path.join(os.getcwd(), folder)

        #   predefine empty list to store image objects
        # images = []
        for _step in xrange(0, len(self.step), int(self._delta_step)):
            filename = "V1_step=%03d"%_step+".png"

            #   assign step and repaint GL widget
            self.__step = int(_step)
            self.__update_GL()

            #   get image of current opengl widget scene
            image = self.OpenGLWidget.grabFrameBuffer()

            #   abs path to image object of current simulation time step
            file_abs_path = os.path.join(animation_folder_abs_path, filename)
            #   save image to file
            image.save(file_abs_path)

        #   create video object
        video = ImageSequenceClip(animation_folder_abs_path, fps=24)#24
        #   write data to video object
        __animation_filename = "animation.avi"

        #   check filename
        __animation_filename = check_filename(__animation_filename)

        video.write_videofile(__animation_filename,codec='mpeg4')

        #   delete animation folder with animation figures
        shutil.rmtree(animation_folder_abs_path)
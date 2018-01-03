"""
Created on 3. mar. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
"""
import os
import psutil
import time
import cProfile
from pprint import pprint
import time
import datetime


import numpy as np
from PyQt4 import QtCore, QtGui


from vtk_widget.vtk_widget import VTKWidget
from simulation_control_widget_ui import Ui_Form
from solver.solver import Solver
from solver.solver_thread_manager import SolverThreadManager
from solver.solver_process_manager import SolverProcessManager
from signals import StatusSignal
from movie_maker.movie_maker import MovieMaker
from signals import SignalSimulationStatus
from MBD_system.solution_data.solution_data import SolutionData


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

        #   set size independent of screen resolution to be equal for all resolutions
        self.setFixedWidth(.2 * self._parent.screen_geometry.width())

        self.ui.simulationStopButton.setEnabled(False)
        self.ui.simulationResetButton.setEnabled(False)
        
        self.ui.forwardButton.setEnabled(False)
        self.ui.backwardButton.setEnabled(False)
        self.ui.playButton.setEnabled(False)

        #   movie maker
        self.movie_maker = None

        #   pool of tasks
        self.pool = QtCore.QThreadPool.globalInstance()
        self.pool.setMaxThreadCount(1)
        
        self.MBD_system = MBD_system
        # self.simulation_control_widget = self

        #   when simulation is finishes
        #   automatically load solution file
        self.ui.loadSolutionFileStatus.setChecked(self.MBD_system.loadSolutionFileWhenFinished)
        #   restore initial conditions
        self.ui.restoreInitialConditionsStatus.setChecked(self.MBD_system.restoreInitialConditionsWhenFinished)

        #    set visualization widget in central widget position
        self.vtkWidget = VTKWidget(MBD_system=self.MBD_system, parent=self._parent)
        self._status = "simulation"  # simulation or animation

        #   set integration method to display
        self.ui.integrationMethodComboBox.addItems(self.MBD_system.integrationMethods)
        _index = self.ui.integrationMethodComboBox.findText(self.MBD_system.integrationMethod)
        self.ui.integrationMethodComboBox.setCurrentIndex(_index)
        
        #    signals
        self.step_num_signal = stepSignal()
        self.energy_signal = EnergySignal()
        self.status_signal = StatusSignal()
        self.signal_simulation_status = SignalSimulationStatus()

        #    use BSM - baumgarte stabilization method
        self.ui.useBSM_checkBox.setChecked(self.MBD_system.use_BSM)

        #    set validators to limit user input
        __validator_dbl = QtGui.QDoubleValidator()
        __validator_int = QtGui.QIntValidator()
        
        #    predefined values
        #    end time
        self.ui.endTime.setValidator(__validator_dbl)
        if self.MBD_system.t_n is not None:
            self.ui.endTime.setText(str(self.MBD_system.t_n))

        #   number of steps
        self.ui.stepsNumber.setValidator(__validator_int)
        if self.MBD_system.stepsNumber is not None:
            self.MBD_system.stepsNumber = int(self.MBD_system.stepsNumber)
            self.ui.stepsNumber.setText(str(self.MBD_system.stepsNumber))
        
        #   simulation time step
        self.step = 0
        self._step = 0

        #   animation properties
        self._delta_step = None
            
        #   H min
        self.MBD_system.evaluate_Hmin()

        #   H max
        self.ui.Hmax.setText(str(self.MBD_system.Hmax))
        self.ui.Hmax.setValidator(__validator_dbl)
        
        #   H min
        self.ui.Hmin.setText(str(self.MBD_system.Hmin))
        self.ui.Hmin.setValidator(__validator_dbl)
        
        #   H contact
        self.ui.Hcontact.setText(str(self.MBD_system.Hcontact))
        
        #    abs tol
        self.ui.absTol.setText(str(self.MBD_system.absTol))
        self.ui.absTol.setValidator(__validator_dbl)
        
        #    rel tol
        self.ui.relTol.setText(str(self.MBD_system.relTol))
        self.ui.relTol.setValidator(__validator_dbl)
        
        #   tolerance for newton differences
        self.ui.TOL_dq_i.setText(str(self.MBD_system.TOL_dq_i))
        
        #   tolerance for constraint equations C
        self.ui.TOL_C.setText(str(self.MBD_system.TOL_C))

        #    create solver thread
        #   solver thread
        self._solver_thread = QtCore.QThread()
        self._solver_thread.start()
        if not MBD_system.monte_carlo:
            self.solver = Solver(MBD_system, parent=self)

        else:
            #   jobs of threads
            # self.solver = SolverThreadManager(MBD_system=self.MBD_system, parent=self)
            #   jobs as processes
            self.solver = SolverProcessManager(MBD_system=self.MBD_system, parent=self)

        self.solver.moveToThread(self._solver_thread)

        # self.solver_process = psutil.Process(id(self._solver_thread))

        #   timer to update cpu, memory data in simulation control widget
        self._data_timer = QtCore.QTimer(self)
        self._data_timer.timeout.connect(self._update_cpu_memory_data)

        #   display update
        self._update_display_type = self.MBD_system._update_display_type
        _index = self.ui.updateDisplay_comboBox.findText(self._update_display_type)
        self.ui.updateDisplay_comboBox.setCurrentIndex(_index)
        #   available options
        self._update_display_types = ["dt",
                                      "step"]
        #   update display on every i-th simulation step
        if hasattr(self.solver.analysis, "update_opengl_widget_every_Nth_step"):
            self.ui.updateStep_lineEdit.setText(str(int(self.solver.analysis.update_opengl_widget_every_Nth_step)))
            self._delta_step = self.solver.analysis.update_opengl_widget_every_Nth_step

        self.ui.updateStep_lineEdit.setValidator(__validator_int)
        self.ui.updateStep_lineEdit.setText(str(int(self.MBD_system.updateEveryIthStep)))

        #   update display on dt of simulation time
        #   default value
        if self.MBD_system._dt == 0:
            self._dt = self.MBD_system.t_n / 100.
        else:
            self._dt = self.MBD_system._dt

        self.ui.updateDt_lineEdit.setValidator(__validator_dbl)
        self.ui.updateDt_lineEdit.setText(str(self._dt))

        self.ui.currentStep_lineEdit.setEnabled(False)
        self.ui.currentStep_lineEdit.setValidator(__validator_int)

        #   profiler
        self.profile = cProfile.Profile()

        #   analysis type
        self.ui.analysisTypeComboBox.currentIndexChanged.connect(self._analysis_type_changed)

        #   set analysis type to display it in combobox
        if self.MBD_system.analysis_type is not None:
            _index = self.ui.analysisTypeComboBox.findText(QtCore.QString(self.MBD_system.analysis_type.title()))
            self.ui.analysisTypeComboBox.setCurrentIndex(_index)

        #   tab widget of simulation control widget
        self.ui.tabWidget.setCurrentIndex(0)

        #   video maker attribute object
        self.video_maker = None

        #   initial start time and date
        self.ui.simulationStartTime_dateTimeEdit.setDate(QtCore.QDate().currentDate())

        #   progress bar style
        css_file = open(os.path.join(os.path.dirname(os.path.realpath(__file__)), "progress_bar_style_no_error.css"), 'r')
        self.progress_bar_no_error_css_stype = css_file.read()
        self.ui.simulation_progressBar.setStyleSheet(self.progress_bar_no_error_css_stype)

        css_file = open(os.path.join(os.path.dirname(os.path.realpath(__file__)), "progress_bar_style_error.css"), 'r')
        self.progress_bar_error_css_stype = css_file.read()

        #   signals
        self.connect_signals()

        #   buttons
        self.connect_buttons()

        #   connections
        self.connect_UIWidgetItems2variables()

        #   check if solutionfile is defined and load it
        if self.MBD_system.solutionFilename is not None:
            solutionFilePath = os.path.join(self.MBD_system.MBD_folder_abs_path, self.MBD_system.solutionFilename)
            if os.path.isfile(solutionFilePath):

                self.__automaticaly_load_solution_file(filename=solutionFilePath)
            else:
                print "Solution File not found at path %s: " % self.MBD_system.MBD_folder_abs_path
                print "Check Filename: %s" % self.MBD_system.solutionFilename

    def setMBDsystem(self, MBD_system):
        """

        :param MBD_system:
        :return:
        """
        self.MBD_system = MBD_system

    def connect_signals(self):
        """

        :return:
        """
        self.ui.simulationStartButton.clicked.connect(self.simulationStart)
        self.ui.simulationStartButton.clicked.connect(self.solver.start_solver)

        self.ui.simulationStopButton.clicked.connect(self.simulationStop)
        self.ui.simulationStopButton.clicked.connect(self.solver.stop_solver)

        if self.solver.analysis is not None:
            self.solver.analysis.finished_signal.signal_finished.connect(self.simulationFinished)

            #   signal time integration error
            self.solver.analysis.error_time_integration_signal.signal_time_integration_error.connect(self._time_integration_error)

            self.solver.analysis.refresh_signal.signal_refresh.connect(self.vtkWidget.refresh)
            self.solver.analysis.refresh_signal.signal_refresh.connect(self._update_simulation_info)

    def connect_buttons(self):
        """

        :return:
        """
        self.ui.simulationResetButton.clicked.connect(self.simulationReset)

        self.ui.backwardButton.clicked.connect(self.animation_backward)
        self.ui.forwardButton.clicked.connect(self.animation_forward)
        self.ui.playButton.clicked.connect(self.animationPlay)

    def connect_UIWidgetItems2variables(self):
        """

        :return:
        """
        self.ui.Hmax.textChanged.connect(self.__update_Hmax)
        self.ui.Hmin.textChanged.connect(self.__update_Hmin)

        self.ui.currentStep_lineEdit.textChanged.connect(self.__update_currentStep)
        self.ui.updateStep_lineEdit.textChanged.connect(self.__update_updateStep)
        self.ui.endTime.textChanged.connect(self.__update_endTime)

        self.ui.updateDt_lineEdit.textChanged.connect(self.update_dt)

        #   change integration method
        self.ui.integrationMethodComboBox.currentIndexChanged.connect(self.selectedIntegrationMethod)

        self.ui.loadSolutionFileStatus.stateChanged.connect(self.__update_loadSolutionFileWhenFinished)

        self.ui.updateDisplay_comboBox.currentIndexChanged.connect(self.selectedUpdateDiyplayMethod)

        self._parent.tree_view_widget.create_animation_file.signal_createAnimationFile.connect(self._create_animation_file)

    def resizeEvent(self, event):
        """
        
        """
        # print "resized"
        # print "w =", self.frameGeometry().width()
        # print "h =", self.frameGeometry().height()
    
    def update_dt(self):
        """
        
        """
        value, sucess = self.ui.updateDt_lineEdit.text().toFloat()
        if sucess:
            self._dt = value

    def _analysis_type_changed(self):
        """
        Function assignes new value to object attribute if it is changed by user in combo box
        :return:
        """
        self.MBD_system.analysis_type = self.ui.analysisTypeComboBox.currentText().toLower()
        
        if str(self.ui.analysisTypeComboBox.currentText()).lower() == "kinematic":
            self.ui.integrationMethodComboBox.setEnabled(False)
        else:
            self.ui.integrationMethodComboBox.setEnabled(True)

    @QtCore.pyqtSlot()
    def profile_functions(self):
        """
        Function profiles main functions that are used in numerical integration
        :return:
        """
        print "profile here"
        self.profile.enable()
        self.solver.analysis.DAE_fun.preprocessing()
        self.profile.disable()
        self.profile.print_stats()

    def _time_integration_error(self):
        """

        :return:
        """
        if self.solver.analysis._error:
            QtGui.QMessageBox.critical(self._parent, "Error!", "Hmin exceeded! Procedure failed!", QtGui.QMessageBox.Ok, QtGui.QMessageBox.NoButton,QtGui.QMessageBox.NoButton)

        if self.solver.analysis.DAE_fun.error:
            QtGui.QMessageBox.critical(self._parent, "Error!", "Violation of constraints! Procedure failed!", QtGui.QMessageBox.Ok, QtGui.QMessageBox.NoButton, QtGui.QMessageBox.NoButton)

        #   update color of progress bar
        # template_css = """
        #                 QProgressBar::chunk
        #                 {
        #                 background-color: #d7801a;
        #                 text-align: center
        #                 }
        #                 QProgressBar
        #                 {
        #                     border: 2px solid grey;
        #                     border-radius: 5px;
        #                     text-align: center;
        #                 }
        #                 """
        # css = template_css % "red"

        self.ui.simulation_progressBar.setStyleSheet(self.progress_bar_error_css_stype)

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
        self._step = int(self.ui.currentStep_lineEdit.text())
        if self._step in self.MBD_system.loaded_solution._step_num_solution_container:
            # self.MBD_system.update_coordinates_and_angles_of_all_bodies(self.q[self._step, :])
            # self.MBD_system.update_simulation_properties(time=self.t[self._step], step_num=self._step)
            self.step_num_signal.signal_step.emit(self._step)
            self._refresh(step=self._step)

    def __automaticaly_load_solution_file(self, solution_data_object_ID=None, filename=None):
        """
        Function loads solution data from input:
        if input is filename solution data is read from specified file
        if input is solution data object, solution is read from object attribute solution_data
        :param filename:                a filename of the solution data file
        :param solution_data_object:    a solution data object with attribute solution_data
        :return: None
        """
        if self.MBD_system.loadSolutionFileWhenFinished or self.MBD_system.solutionFilename is not None:
            self.load_solution_file(solution_data_object_ID, filename)

    def load_solution_file(self, solution_object_id=None, filename=None):
        """
        Function loads solution data from object (first) or from file (second)
        """
        if solution_object_id is None:
            #   assign solution data from solution data object to new variable
            if isinstance(self.solver.analysis._solution_data, SolutionData):
                solution_data = self.solver.analysis._solution_data.load_solution_data()

            #   load solution data from file specified
            elif self.solver.analysis._solution_data is None:
                solution_data = SolutionData(parent=self.MBD_system._children[self.MBD_system.dict_of_object_group_indexes["Solution"]])

                if self.MBD_system.solutionMBDFormat != "new":
                    solution_data.set_format(self.MBD_system.solutionMBDFormat)

                solution_data.read_file(filename)
                self.MBD_system.loaded_solution = solution_data

            else:
                raise Warning, "Solution object not constructed!"
        elif filename is not None:
            pass
        else:
            print "Solution data not loaded!"

        # if filename is not None and solution_object_id is None:
        #     solution_data = self.solver.analysis.load_simulation_solution_from_file(filename)

        #   assign a solution data object to pointer of object attribute

        for sol in self.MBD_system.solutions:
            if id(sol) == solution_object_id:
                self.MBD_system.loaded_solution = sol
                if sol.solution_data is None and filename is not None:
                    sol.read_file(sol.filename)
                    sol.loaded = True

        # self.step = solution_data[:, 0]
        # self.energy = solution_data[:, 1]
        # self.error = solution_data[:, 2]
        # self.dt = solution_data[:, 3]
        # self.t = solution_data[:, 4]
        # self.q = solution_data[:, 5:]
        self.step = self.MBD_system.loaded_solution._step_num_solution_container
        self._step = 0
        self.t = self.MBD_system.loaded_solution._t_solution_container
        self._status = "animation"
        
        self.ui.forwardButton.setEnabled(True)
        self.ui.backwardButton.setEnabled(False)
        self.ui.playButton.setEnabled(True)
        self.ui.currentStep_lineEdit.setEnabled(True)
        
        self.ui.currentStep_lineEdit.setText(str(int(self._step)))
        self.ui.currentStep_lineEdit.setMaxLength(4)
        
        self.ui.solutionFileLoaded_display.setText(filename)
        self.ui.numberOfSteps_lineEdit.setText(str(int(len(self.MBD_system.loaded_solution._step_num_solution_container)-2)))
        
        # self.status_signal.signal_status.emit("Animation")

        if self.MBD_system.restoreInitialConditionsWhenFinished:
            self.simulationReset()

        #   tab widget of simulation control widget
        self.ui.tabWidget.setCurrentIndex(1)

    def _refresh(self, step=None):
        """
        Function updates values of vector q from solution object to display it
        """
        if step is None:
            step = self._step

        self.ui.currentStep_lineEdit.setText(str(int(step)))

        self.MBD_system.time = self.MBD_system.loaded_solution._t_solution_container[int(step)]
        self.MBD_system.step_num = step

        t = self.MBD_system.loaded_solution._t_solution_container[int(step)]
        q = self.MBD_system.loaded_solution._q_solution_container[int(step)]
        self.MBD_system.update_coordinates_and_angles_of_all_bodies(t, q, step=step)
        self.vtkWidget.refresh(step=int(step))

        #   step signal
        self.step_num_signal.signal_step.emit(int(step))

        #    energy data signal
        _energy = self.MBD_system.loaded_solution._mechanical_energy_solution_container[step]
        _energy_delta = _energy - self.MBD_system.loaded_solution._mechanical_energy_solution_container[int(step)-1]
        self.energy_signal.signal_energy.emit(_energy, _energy_delta)

    def _refresh_measure_graph(self):

        for measure in self.MBD_system.measures:
            measure._paintGL()

    def _update_cpu_memory_data(self):
        """
        
        """
        p = psutil.Process(self.solver.pid)
        #   update CPU stats
        self.ui.cpu_doubleSpinBox.setValue(p.cpu_percent())
        #   update RAM stats
        _memory = psutil.virtual_memory()
        self.ui.memory_doubleSpinBox.setValue(_memory.percent)

    def animation_forward(self):
        """
        Go to next solution time step
        """
        if self.ui.updateDisplay_comboBox.currentText() == "step":
            if (self._step + self._delta_step) <= int(self.step[-1]):
                self._step += int(self._delta_step)
                self._refresh()
                self.ui.backwardButton.setEnabled(True)
            else:
                self.ui.forwardButton.setEnabled(False)

        else:
            print "animation_forward()-TODO"

    def animation_backward(self):
        """
        Return to prevouos solution time step
        """
        if self.ui.updateDisplay_comboBox.currentText() == "step":
            if (self._step - self._delta_step) >= self.step[0]:
                self._step -= int(self._delta_step)
                self._refresh()
                self.ui.forwardButton.setEnabled(True)
            else:
                self.ui.backwardButton.setEnabled(False)

    def take_snapshot(self):
        """
        Create a snapshot
        """
        captured_figure = self.opengl_widget.takeSnapShot()
        captured_figure.save(self.solver.analysis.screenshot_filename_abs_path + '.png', 'png')

    def selectedIntegrationMethod(self, int):
        """
        Assign a selected integration method to object attribute
        """
        self.MBD_system.integrationMethod = str(self.ui.integrationMethodComboBox.currentText())

    def selectedUpdateDiyplayMethod(self):
        """

        :return:
        """
        self._update_display_type = self.ui.updateDisplay_comboBox.currentText()
        self.t = 0.

    def __update_updateStep(self):
        """

        :return:
        """
        self.solver.analysis.update_opengl_widget_every_Nth_step = int(self.ui.updateStep_lineEdit.text())
        self._delta_step = int(self.ui.updateStep_lineEdit.text())
        self.MBD_system.updateEveryIthStep = self._delta_step

    def __update_Hmax(self):
        """

        :return:
        """
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
        """

        :return:
        """
        self.endTime = float(self.ui.endTime.text())
        self.MBD_system.t_n = self.endTime

    def setWindowFlags(self, flags):
        super(SimulationControlWidget, self).setWindowFlags(flags)

    def _update_simulation_info(self):
        """

        :return:
        """
        #   simulation progress
        if not np.isnan(self.solver.analysis.progress):
            self.ui.simulation_progressBar.setValue(int(self.solver.analysis.progress * 100))

        #   elapsed time
        elapsed_time = time.time() - self.solver.analysis.start_time_simulation_info_UTC
        m, s = divmod(elapsed_time, 60)
        h, m = divmod(m, 60)
        _time = QtCore.QTime()
        _time.setHMS(h, m, s)
        self.ui.simulationElapsedTime_timeEdit.setTime(_time)

        #   time to finish
        if self.solver.analysis.progress > 0:
            finish_time = (elapsed_time / self.solver.analysis.progress) - elapsed_time
            m, s = divmod(finish_time, 60)
            h, m = divmod(m, 60)
            _time = QtCore.QTime()

            if np.isnan(h):
                h = 0
            if np.isnan(m):
                m = 0
            if np.isnan(s):
                s = 0

            _time.setHMS(int(h), int(m), int(s))
            self.ui.simulationTimeToFinish_timeEdit.setTime(_time)

    def simulationStart(self):
        """

        :return:
        """
        if hasattr(self.solver.analysis, "running"):
            if not self.solver.analysis.running:
                self._data_timer.start(1000)

                self.job_info_time_started = time.clock()
                self.solver.analysis.running = True
                self.solver.analysis.stopped = False
                self.solver.analysis.finished = False

        self.ui.simulationStartButton.setEnabled(False)
        self.ui.simulationResetButton.setEnabled(False)
        self.ui.simulationStopButton.setEnabled(True)

        #   set simulation start time to display it in a widget
        self.solver.analysis.start_time_simulation_info = QtCore.QDateTime.currentDateTime()
        self.ui.simulationStartTime_dateTimeEdit.setDateTime(self.solver.analysis.start_time_simulation_info)

        #   start simulation timer
        self.elapsed_time = QtCore.QElapsedTimer()
        self.elapsed_time.start()

    def simulationStop(self):
        """

        :return:
        """
        if hasattr(self.solver.analysis, "running"):
            if self.solver.analysis.running:
                self.solver.analysis.stopped = True

        self.ui.simulationStartButton.setEnabled(True)
        self.ui.simulationStopButton.setEnabled(False)
        self.ui.simulationResetButton.setEnabled(True)

    def simulationFinished(self):
        """

        :return:
        """
        self.ui.simulationStopButton.setEnabled(False)
        self.ui.simulationResetButton.setEnabled(True)
        self.job_info_time_finished = time.clock()

        self._update_simulation_info()
        self.ui.simulation_progressBar.setValue(100)

        QtGui.QMessageBox.information(self._parent, "Finished", "Simulation finished successfully!", QtGui.QMessageBox.Ok, QtGui.QMessageBox.NoButton, QtGui.QMessageBox.NoButton)

    def simulationReset(self):
        """

        :return:
        """
        self.solver.analysis.restore_initial_condition()
        self.ui.simulationResetButton.setEnabled(False)
        self.ui.simulationStartButton.setEnabled(True)
        self.signal_simulation_status.signal_simulation_status.emit("Ready")

        self.ui.simulation_progressBar.setValue(0)
        _time = QtCore.QTime()
        _time.setHMS(0, 0, 0)
        self.ui.simulationElapsedTime_timeEdit.setTime(_time)
        self.ui.simulationTimeToFinish_timeEdit.setTime(_time)

        self.vtkWidget.refresh(step=0, h=self.MBD_system.Hmax)

    def animationPlay(self):
        """
        Function plays animation when solution data is loaded
        :return:
        """
        print "animationPlay()"
        if self._update_display_type == "step":
            self.animationPlay_step()

        if self._update_display_type == "dt":
            self.animationPlay_dt()

    def animationPlay_step(self):
        """
        
        :return: 
        """
        print "animationPlay_step()"
        print "self.ui.updateStep_lineEdit.text() =", self.ui.updateStep_lineEdit.text()
        self._delta_step = int(self.ui.updateStep_lineEdit.text())
        for _step in xrange(0, len(self.step), self._delta_step):
            print "_step =", _step
            self._step = int(_step)
            self._refresh()

            time.sleep(self.ui.playbackSpeed_doubleSpinBox.value()*1E-2)

    def animationPlay_dt(self):
        """
        
        :return: 
        """
        print "animationPlay_dt()"
        print "self.t[-1] =", self.t[-1]
        print "self._dt =", self._dt
        t = np.arange(0., self.t[-1], self._dt)
        print "t =", t
        for i in xrange(0, len(t)):
            print i
            indx = np.argmin(abs(self.t - t[i]))
            self._step = self.step[indx]
            self._refresh()

            time.sleep(self.ui.playbackSpeed_doubleSpinBox.value()*1E-2)

    def _create_animation_file(self):
        """
        Function creates video file of animation
        :return:
        """
        #   video maker thread to create video
        self.movie_maker = MovieMaker(parent=self)
        self.movie_maker.setRenderWindow(self.vtkWidget.GetRenderWindow())
        self.movie_maker.setInputConnection()
        self.movie_maker.fps = self._parent.preferences_widget.ui.fps_spinBox.value()
        self.movie_maker.fps = 24

        self.movie_maker.movieWriter.SetRate(self.movie_maker.fps)
        self.movie_maker.movieWriter.Start()

        if self._update_display_type == "step":
            self._delta_step = int(self.ui.updateStep_lineEdit.text())
            for _step in xrange(0, len(self.MBD_system.loaded_solution._step_num_solution_container), int(self._delta_step)):
                #   assign step and repaint visualization widget
                self._refresh(step=_step)

                #   get image of current vtk widget scene
                self.movie_maker.windowToImageFilter.Modified()
                self.movie_maker.movieWriter.Write()

                #   wait
                time.sleep(1E-3)

        if self._update_display_type == "dt":
            t = np.arange(0, self.t[-1], self._dt)

            for i in xrange(0, len(t)):
                #   find index of
                step = np.argmin(abs(self.t - t[i]))

                #   assign step and repaint visualization widget
                self._refresh(step=step)

                #   get image of current vtk widget scene
                self.movie_maker.windowToImageFilter.Modified()
                self.movie_maker.movieWriter.Write()

                #   wait
                time.sleep(1E-3)

        #   starts video maker in new thread
        self.movie_maker.movieWriter.End()


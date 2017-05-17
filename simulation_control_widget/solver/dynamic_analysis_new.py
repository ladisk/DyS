"""

created by: lskrinjar
date of creation: 10/05/2016
time of creation: 19:24
"""
import sys
import copy
import datetime
import fractions
import inspect
import itertools
import logging
import os
import time


import numpy as np
from PyQt4 import QtCore, QtGui


from DAE_fun import DAEfun
from MBD_system.check_filename import check_filename
from MBD_system.q2R_i import q2R_i
from MBD_system.q2dR_i import q2dR_i
from MBD_system.q2dtheta_i import q2dtheta_i
from MBD_system.solution_data.solution_data import SolutionData
from signals import ErrorTimeIntegrationSignal
from signals import SolutionSignal
from signals import SignalSimulationStatus
from signals import StepSignal
from signals import RefreshSignal
from signals import FinishedSignal
from signals import SaveScreenshot
from signals import SolutionFilenameSignal
from signals import EnergySignal


class DynamicAnalysis(object):
    """
    classdocs
    """
    __simulation_id = itertools.count(0)


    def __init__(self, MBD_system, parent=None):
        """
        Constructor
        """
        # super(DynamicAnalysis, self).__init__(parent)

        #   parent
        self._parent = parent

        #   states of simulation
        self.running = False
        self.stopped = False
        self.finished = False
        self.failed = False
        self._error = False
        self.errorControl = True

        #   signals
        self.step_signal = StepSignal()
        self.finished_signal = FinishedSignal()
        self.refresh_signal = RefreshSignal()
        self.energy_signal = EnergySignal()
        self.error_time_integration_signal = ErrorTimeIntegrationSignal()
        self.filename_signal = SolutionFilenameSignal()
        self.save_screenshot_signal = SaveScreenshot()
        self.solution_signal = SolutionSignal()
        self.signal_simulation_status = SignalSimulationStatus()

        #   display active
        self._active_scene = True

        #   MBD system object
        self.MBD_system = MBD_system
        #    copy MBD object
        self._MBD_system = copy.copy(self.MBD_system)
        #    DAE fun object
        self.DAE_fun = DAEfun(self._MBD_system, parent=self)

        #   define simulation attributes
        self.simulation_id = None
        self.h_contact = None
        self.t_0 = 0.
        self.t_n = None
        self.stepsNumber = None
        self.absTol = None
        self.relTol = None
        self.absError = None
        self.relError = None
        self.progress = 0
        self._energy_t = 0.
        self._energy_delta = 0.
        self.h_min_level = None
        self.h_max_level = None
        self.h_contact_level = None

        #    simulation settings - properties
        # self.integrationMethod = "RKF"  #    self.MBD_system.integrationMethod

        #    update opengl widget every N-th step
        self._dt = self.MBD_system._dt
        self.step = 0
        self.t = 0.
        self.update_opengl_widget_step_count = 0
        self.update_opengl_widget_every_Nth_step = self.MBD_system.updateEveryIthStep
        self.t_constant_time_step = None

        self.save_screenshot_step_count = 0
        self.save_screenshot_every_Nth_step = 2 * self.update_opengl_widget_every_Nth_step

        #    info properties
        self.start_time_simulation_info_in_sec_UTC = []
        self.start_time_simulation_info = []
        self.end_time_simulation_info_in_sec_UTC = []
        self.end_time_simulation_info = []
        self.simulation_time_info_in_sec_UTC = []
        self.simulation_time_info = []
        self.stop_time_simulation_info_in_sec_UTC = []
        self.stop_time_simulation_info = []

        self.FLAG = 1
        self.FLAG_contact = 0
        self.FLAG_contact_finish= 0

        #   level of integration step size - time
        self.h_level = 0
        self.dh = 0
        self._steps_after_contact = 0
        self._dh_constant_steps = 0

        #   solution data settings
        #   write to file during simulation
        self._append_to_file_every_Nth_step = 20
        self._append_to_file_step = 0
        self._append_to_file_during_simulation = True
        self._append_to_file_keep_last_Nth_steps_in_memory = 10

    def _solution_containers(self):
        """
        Function initializes solution containers object attributes to save solution data during numerical integration process
        :return:
        """
        #    solution container variable
        self.q0 = self._MBD_system.evaluate_q0()
        self.q_sol_matrix = np.array([self.q0])
        self.t_vector = np.array([0])
        self.R_vector = np.array([0])
        self.step_counter_vector = np.array([0]).astype(int)
        self.step_size_vector = np.array([self._MBD_system.Hmax])
        self.h_level_vector = np.array([0]).astype(int)
        self.energy_data_vector = np.array([self._mechanical_energy(self.q0)])
        self.contact_status_vector = np.array([0]).astype(int)

    def start_solve(self):
        """
        Start solver (time integration process)
        """
        #   set step size
        self.Hmin = self._MBD_system.Hmin
        self.Hmax = self._MBD_system.Hmax
        self.Hcontact = self._MBD_system.Hcontact
        self.Houtput = self._MBD_system.Houtput

        #   set step level
        self.h_min_level = 0
        self.h_max_level = int(np.trunc(np.log2(self.Hmax / self.Hmin)))#np.trunc
        self.h_contact_level = np.trunc(np.log2(self.Hmax / self.Hcontact))
        print "self.h_contact_level =", self.h_contact_level
        self.h_output_level = np.trunc(np.log2(self.Hmax / self.Houtput))

        #   solution id
        self.simulation_id = self.__simulation_id.next()

        #   create solution data object to store data at each simulation start
        self._solution_data = SolutionData(MBD_item=self.MBD_system)
        #   clear temp solution file
        self._solution_data._clear_file()

        self.MBD_system.solutions.append(self._solution_data)

        self.signal_simulation_status.signal_simulation_status.emit("Running")
        logging.getLogger("DyS_logger").info("Simulation started, simulation id: %s"%self.simulation_id)

        self.running = True
        self.stopped = False
        self.finished = False

    def stop_solve(self):
        """
        Stop solver (time integration process) - interput by user
        """
        #   get solution data to one variable
        solution_data = self.collect_solution_data()
        #   append solution data to file
        self._solution_data._append_data_to_temp_file(solution_data[0:-1, :])

        self.stopped = True
        self.running = False

        self.stop_time_simulation_info_in_sec_UTC = time.time()
        self.stop_time_simulation_info = datetime.datetime.fromtimestamp(self.stop_time_simulation_info_in_sec_UTC).strftime("%H:%M:%S %d.%m.%Y")

        self.signal_simulation_status.signal_simulation_status.emit("Stopped")
        logging.getLogger("DyS_logger").info("Simulation stopped by user at simulation time:%s", self.t)

    def finished_solve(self):
        """
        Method sets some parameters when integration process is finished
        :return:
        """
        self._refresh(t=self.t, step=self.step)
        if self.finished:
            self.solution_signal.solution_data.emit(id(self._solution_data), self._solution_data._name)
            logging.getLogger("DyS_logger").info("Simulation finished successfully!")
            self.finished_signal.signal_finished.emit()
            self.signal_simulation_status.signal_simulation_status.emit("Finished")

        #   set booleans for states
        self.finished = True
        self.running = False

    def simulation_failed(self):
        """

        :return:
        """
        print "simulation_failed()"
        self.error_time_integration_signal.signal_time_integration_error.emit()
        self.signal_simulation_status.signal_simulation_status.emit("Failed")
        if self._error or self.failed:
            logging.getLogger("DyS_logger").info("Simulation failed!")
        if self.DAE_fun.error:
            logging.getLogger("DyS_logger").info("Violation of kinematic constraint!")

    def _info(self, t, w):
        """

        :return:
        """
        #    update coordinates and angles of all bodies
        self._update_coordinates_and_angles_of_all_bodies(t, w)

        #    update display
        self.refresh(t, w)

        # self.step += 1
        # step_counter = np.append(step_counter, self.step)

        # self.update_opengl_widget_step_count += 1
        self.save_screenshot_step_count += 1

        #    step number signal
        self.step_signal.signal_step.emit(self.step)

        #    energy data signal
        self.energy_signal.signal_energy.emit(self._energy_t, self._energy_delta)

        #   simulation progress
        if self.t_n is not None:
            self.progress = t / self.t_n
        else:
            self.progress = self.step / self.stepsNumber

    def _track_data(self, h, t, q):
        """

        :return:
        """
        #    save values at next time step to matrix of solution (q_matrix)
        if self.FLAG_contact in [0, 1]:
            #   append solution at i-th time to solution matrix
            self.q_sol_matrix = np.vstack((self.q_sol_matrix, np.array(q)))

            self.step_size_vector = np.append(self.step_size_vector, h)
            self.step += 1
            step_num = self.step

            #   append current step number to step counter history array
            self.step_counter_vector = np.append(self.step_counter_vector, self.step)
            #    add to vector
            self.t_vector = np.append(self.t_vector, t)
            self.R_vector = np.append(self.R_vector, self.absError)
            #   append current mechanical energy to array
            self.energy_data_vector = np.append(self.energy_data_vector, self._energy_t)

            #   contact status
            self.contact_status_vector = np.append(self.contact_status_vector, self.FLAG_contact)

            #   step size level
            self.h_level_vector = np.append(self.h_level_vector, self.h_level)

            for contact in self.MBD_system.contacts:
                contact._track_data(self.step, h, t)

            for force in self.MBD_system.forces:
                force._track_data(self.step, h, t)

            for spring in self.MBD_system.springs:
                spring._track_data(self.step, h, t)

            for joint in self.MBD_system.joints:
                joint._track_data(self.step, h, t)

            for measure in self.MBD_system.measures:
                measure._track_data(self.step, h, t, q)

            #   check is step number is reached to append data to file
            if self._append_to_file_step == self._append_to_file_every_Nth_step-1 and self._append_to_file_during_simulation:
                #   collect solution data to one variable
                solution_data = self.collect_solution_data()
                #   append solution data to file
                self._solution_data._append_data_to_temp_file(solution_data[0:self._append_to_file_keep_last_Nth_steps_in_memory, :])
                # print "solution data appended!"
                #   remove all but last value in array/matrix
                self._reset_solution_data_after_append_to_temp_file()

                #   counter
                self._append_to_file_step = 0
            else:
                self._append_to_file_step += 1

        #    reduce step size and continue integration from previous step
        elif self.absError > self.absTol:
            #    go to one before last solution
            q = self.q_sol_matrix[-1, :]
            #   remove last stored value in matrix
            # self.q_sol_matrix = np.delete(self.q_sol_matrix, -1, axis=0)
            #    store time that is too large
            self._t_FLAG1 = t
            t = self.t_vector[-1]  # , 0

#             #    reduce step size
#             h = 0.5 * h
#             #    level os step size
#             self.t_level += 1
            # self.step -= 1
            # step_num  = self.step
            # step_counter[-1] = self.step

            self.energy_data_vector[-1] = self._energy_t

        return h, t, q

    def _evaluate_contacts(self, t, q):
        """
        Function solves contacts. Finds overlap pairs of (sub)AABB and calculates actual distances between node
        and line (edge). If distance is below userdefined limit contact is present and contact equations are solved.
        :type self: object
        Args:
            q    vector of positions and velocities (translational and rotational)
        Returns:
            q_   vector of new calculated velocities (positions are equal before and after contact
            status - 0 - contact is not detected continue with integration to next time step
                     1 - contact is detected (has happened) return to previous time step solution and split time step
                         size
                     2 - contact is detected (will happen) split time step size
        Raises:
            None
        """
        #   predefine array of zeros
        self.contact_status_list = np.zeros(len(self.MBD_system.contacts), dtype="int")
        #    check for contact of every contact pair
        for i, contact in enumerate(self.MBD_system.contacts):
            # print "contact._contact_point_found =", contact._contact_point_found, "t =", t
            if not contact._contact_point_found:
                #    function recursively loops through all potential AABB objects in AABB tree of each body
                #    and if two AABB objects overlap, a new overlap object is created
                #    that has both overlapping AABB objects
                #   general contact
                if contact._contact_type in ["general", "roughness profile"]:
                    #   adds simulation data to contact object
                    #   reset to empty list a contact object attribute of list of overlap pairs,
                    #   as this is not possible in the next line (function) due to the recursion
                    contact.AABB_list_of_overlap_pairs = []
                    #   check for overlap and if overlap is present build a overlap pair object
                    contact.check_for_overlap(q, contact.AABB_i, contact.AABB_j)

                    #   if list of overlap object pairs is not empty, check for contact
                    if contact.AABB_list_of_overlap_pairs != []:
                        self.contact_status_list[i] = contact.check_for_contact(self.step, t, q)
                    else:
                        contact.no_overlap(self.step, t, q)

                #   clearance joints of contact
                elif contact._contact_type.lower() in contact._contact_types:
                    self.contact_status_list[i] = contact.check_for_contact(self.step, t, q)

                else:
                    self.FLAG = 0
                    QtGui.QMessageBox.critical(self._parent, "Error!", "Contact type not correct: %s"%contact._type+". \nDefine correct contact type!", QtGui.QMessageBox.Ok)
                    raise ValueError, "Contact type not correct: %s"%contact._type, "\n Define correct contact type!"
            else:
                pass
                #    adds simulation data to contact object as this is not possible in the next line because of the recursion
                self.contact_status_list[i] = contact.contact_update(self.step, t, q)

        #    no contacts
        if (self.contact_status_list == 0).all():
            FLAG_contact = 0

            callerframerecord = inspect.stack()[0]      # 0 represents this line
                                                        # 1 represents line at caller
            frame = callerframerecord[0]
            # print "NO CONTACT"
            # print "simulation sleep"

            info = inspect.getframeinfo(frame)
            # print "file:", info.filename
            # print "function:", info.function+"()"
            # time.sleep(100)

            if self.contact_status_vector[-1] == 1 and self.FLAG_contact_finish == 1:
                FLAG_contact = -1
                self.FLAG_contact_finish -= 1

            return FLAG_contact

        #    contact has already happened - reduce time step
        if (self.contact_status_list == -1).any():
            FLAG_contact = -1  # 0-no contact, 1-contact has already happened - reduce time step
            return FLAG_contact

        #    contact
        if (self.contact_status_list == 1).any():
            FLAG_contact = 1
            self.FLAG_contact_finish = 1
            # logging.getLogger("DyS_logger").info("Contact detected at simulation time:%s \nbody_i=%s \nbody_j=%s", t, self.MBD_system.bodies[contact.body_id_i]._name, self.MBD_system.bodies[contact.body_id_j]._name)
            # logging.getLogger("DyS_logger").info("Contact detected at simulation time:%s "%t)
            #   repaint opengl widget

#                 self.repaintGL_signal.signal_repaintGL.emit()

            #   solve contacts - evaluate contact forces
            # self.DAE_fun.solve_contacts(t, q)

            # logging.getLogger("DyS_logger").info("Contact calculation finished")
            return FLAG_contact

        else:
            FLAG_contact = 0

        return FLAG_contact

    def _mechanical_energy(self, q):
        """

        """
        energy = self.MBD_system.evaluate_mechanical_energy(q=q)
        return energy

    def check_time_step_TODO(self, t, h, w):
        """

        :param t:
        :param h:
        :param w:
        :return:
        """
        error = False
        #   continue with current time step - no error
        if (self.absError <= self.absTol) and not (self.DAE_fun.check_C(w, t)):
            # print "NO ERROR"
            self.dh = 0
            if self.FLAG_contact == 0:
                # print "FLAG = 0"
                if self.contact_status_vector[-1] == 1:
                    self.h_level = int(max(self.h_level_vector))

                if self.h_level > 0:
                    self.dh = -1

            if self.FLAG_contact == +1:
                # print "FLAG = 1"
                self.h_level = self.h_contact_level

            if self.FLAG_contact == -1:
                # print "FLAG = -1"
                self.dh = +1

        else:
            # print "ERROR"
            # if self.FLAG_contact == 0:
            #     self.dh = 0
            #
            # if self.FLAG_contact == +1:
            #     print "FLAG = 1"
            #     self.h_level = self.h_contact_level
            if self.FLAG_contact in [-1, +1]:
                # print "self.FLAG_contact =", self.FLAG_contact
                self.dh = +1

        # print "self.h_level =", self.h_level
        # print "self.dh =", self.dh
        self.h_level += self.dh
        # print "self.h_level =", self.h_level
        h = self._evaluate_time_step(self.h_level)

        #   time step size is limited by max value
        if h > self.Hmax:
            h = self.Hmax
            self.h_level = 0

        #   value of time step is smaller than the limit - break computation and raise error
        if h < self.Hmin:
            error = True
            print "Hmin =", self.Hmin
            print "Suggested value is h =", h
            print "absError =", self.absError
            print "absTol =", self.absTol
            print "C(q, t) =", self.DAE_fun.evaluate_C(w, t)
            raise ValueError, "Hmin exceeded!"

        return h, error

    def check_time_step(self, t, h, w):
        """
        Function checks time step size
        """
        # print "self.contact_status_vector[-1]@check_time_step()=", self.contact_status_vector[-1]
        error = False
        self.dh = 0

        #    if at least one contact is present in the MBD system
        if self.FLAG_contact == 1:
            self.dh = 0
            # print 'h_level=', self.h_level
            #   condition to increase speed for calculationg contact conditions
            if self._dh_constant_steps == 0:
                if self.h_level - 1 >= self.h_contact_level:
                    for h_lev in np.arange(self.h_level):
                        if np.around((self.t/self.Hmax*2**(h_lev)), decimals=5)%1 == 0:
                            # dh = -self.h_level + h_lev
                            dh = -1
                            self.h_level += dh
                            if self.h_level < self.h_contact_level:
                                self.h_level = self.h_contact_level
                            break
                else:
                    self.h_level = self.h_contact_level
            else:
                self._dh_constant_steps -= 1

        #    at least one contact has current initial penetration depth larger than initial minimum specified depth by user
        elif self.FLAG_contact == -1 or self.absError > self.absTol:
            if self.h_level + 1 <= self.h_max_level:
                self.dh = 1
                self._dh_constant_steps = 8

            else:
                self.dh = 0
                error = True

                print "Hmin =", self.Hmin
                print "Suggested value is h =", h
                print "absError =", self.absError
                print "absTol =", self.absTol
                print "C(q, t) =", self.DAE_fun.evaluate_C(w, t)

                print Warning, "Hmin exceeded!"

        #   contact is finished and the integration is started from previous time step with smalled step size
        elif self.contact_status_vector[-1] == 1:
            self.h_level = int(max(self.h_level_vector))
            # print 'ELIF h_level', self.h_level

        #    no contact is present at current time step, integration step size can be increased
        else:
            if self._dh_constant_steps == 0:
                if self.h_level - 1 >= self.h_min_level:
                    # print "np.arange(self.h_level) =", np.arange(self.h_level)
                    for h_lev in np.arange(self.h_level):
                        if np.around((self.t/self.Hmax*2**(h_lev)), decimals=5)%1 == 0:
                            # dh = self.h_level - h_lev
                            dh = 1

                            if self.step < self._steps_after_contact:
                                _n = self.step

                            else:
                                _n = self._steps_after_contact

                            _status = self.contact_status_vector[-_n-1:-1]

                            if (_status == np.zeros(_n)).all():
                                self.dh = -dh
                            else:
                                self.dh = 0

                            break
            else:
                self._dh_constant_steps -= 1

        self.h_level += self.dh

        h = self._evaluate_time_step(self.h_level)

        #   time step size is limited by max value
        if h > self.Hmax:
            h = self.Hmax
            self.h_level = 0

        if h < self.Hmin and self.FLAG_contact != 1:
            h = self.Hmin
            error = True
            print "-----------------------------------------"
            for contact in self.MBD_system.contacts:
                print "contact id =", contact.contact_id
                print "delta =", contact._delta

            print "Hmin =", self.Hmin
            print "Suggested value is h =", h
            print "absError =", self.absError
            print "absTol =", self.absTol
            print "C(q, t) =", self.DAE_fun.evaluate_C(w, t)

            print Warning, "Hmin exceeded!"
            #   multiple of time step
#             m = t * (2.**self.t_level) / self.Hmax

#             next_level = fractions.gcd((2**self.t_level), m)
            # print "next_level =", next_level
            # time.sleep(100)

        return h, error

    def _evaluate_time_step(self, level):
        """

        """
        h = self.Hmax / (2**level)
        return h

    def write_simulation_solution_to_file(self):
        """
        Function writes solution data to solution data object
        """
        print "write_simulation_solution_to_file()"
        #   write solution data and info to solution data object
        self._solution_data.set_filetype(self.MBD_system._solution_filetype)

        self.solution_filename = 'solution_data_%02d'%self.simulation_id#+self.MBD_system._solution_filetype

        #    add to temp file
        data = self.collect_solution_data()
        self._solution_data._append_data_to_temp_file(data)

        #    read from temp file
        if os.path.isfile(self._solution_data._temp_filename):
            self._solution_data._solution_data_from_temp_file()

        if self.MBD_system._save_options == "save to new":
            #   check filename if already exists
            self.solution_filename = check_filename(self.solution_filename)

        #   set solution data filename
        self._solution_data.setName(self.solution_filename)

        #   add data to object attribute
        self._solution_data.finished = True

        # print "self.MBD_system._solution_save_options =", self.MBD_system._solution_save_options
        if self.MBD_system._solution_save_options != "discard":
            #   write data to file
            self._solution_data.write_to_file()

            #   emit signal
            # self.solution_signal.solution_data.emit(id(self._solution_data), self._solution_data._filename)

        #   write contact values of every simulation time step to file
        for contact in self.MBD_system.contacts:
            if contact._solution_save_options != "discard":
                contact.write_to_file()

        for joint in self.MBD_system.joints:
            if joint._solution_save_options != "discard":
                joint.write_to_file()

    def collect_solution_data(self):
        """

        :return:
        """
        #   join in one ndarray (matrix)
        solution_data = np.hstack((np.array([self.step_counter_vector]).T,
                                   np.array([self.energy_data_vector]).T,
                                   np.array([self.R_vector]).T,
                                   np.array([self.step_size_vector]).T,
                                   np.array([self.t_vector]).T,
                                   np.array([self.h_level_vector]).T,
                                   np.array([self.contact_status_vector]).T,
                                   self.q_sol_matrix))

        return solution_data

    def _reset_solution_data_after_append_to_temp_file(self):
        """
        Function removes all elements in array except last item (row)
        :return:
        """
        for data_name in ["step_counter_vector", "energy_data_vector", "R_vector", "step_size_vector", "t_vector", "h_level_vector", "contact_status_vector"]:
            data = getattr(self, data_name)
  
            new_data = np.array(data[self._append_to_file_keep_last_Nth_steps_in_memory:])

            setattr(self, data_name, new_data)

        #   solution matrix
        new_data = self.q_sol_matrix[self._append_to_file_keep_last_Nth_steps_in_memory:,:]
        self.q_sol_matrix = new_data

    def save_solution_data(self):
        """

        :return:
        """
        if os.path.isfile(self._solution_data._temp_filename):
            self._solution_data._solution_data_from_temp_file()
        # self._solution_data.add_data(data)

    def load_simulation_solution_from_file(self, filename):
        """
        Function load solution_data from file.
        """
        # self.finished_signal.signal_finished.emit("Animation")
        self.signal_simulation_status.signal_simulation_status.emit("Animation")
        if os.path.exists(filename):
            self.solution_data = np.loadtxt(str(filename), skiprows=2)
            return self.solution_data

    def _update_coordinates_and_angles_of_all_bodies(self, t, q):
        """
        Function updates MBD_system data before display
        """
        # self.MBD_system.update_coordinates_and_angles_of_all_bodies(q)
        self.MBD_system.update_vtk_data(t, q)

    def refresh(self, t, q):
        """
        Update visualization widget at preset time steps
        """
        if self._active_scene:
            #   update display on every i-th step
            if self._parent._update_display_type == "step":
                self.update_opengl_widget_step_count += 1
                if self.update_opengl_widget_step_count == self.update_opengl_widget_every_Nth_step or self.finished or self.stopped:
                    self._refresh(t, self.step)

                    self.time = t

                    self.save_updated_screenshot_to_file(self.time)

                    #    reset step counter to update opengl widget
                    self.update_opengl_widget_step_count = 0

            #   update display on simulation time interval dt
            if self._parent._update_display_type == "dt":
                _dt = t - self._dt
                if _dt >= self._parent._dt:
                    self._refresh(t, self.step)

                    #   set new value
                    self._dt = t

        return None

    def _refresh(self, t, step):
        """

        :return:
        """
        #    update visualization widget
        self.MBD_system.update_simulation_properties(time=t, step_num=step)

        #   emit signal
        if self._active_scene:
            self.refresh_signal.signal_refresh.emit()

    def save_updated_screenshot_to_file(self, t):
        """

        """
        if self.MBD_system.save_screenshots:
            if (self.save_screenshot_step_count == self.save_screenshot_every_Nth_step) or (t == self.t_0) or (self.finished):

                self.screenshot_filename_abs_path = os.path.join(self.MBD_system.saved_screenshots_folder_abs_path, str("step=" + "%010d" % self.step))
                self.save_screenshot_signal.signal_saveScreenshot.emit()
                self.save_screenshot_step_count = 0

    def restore_initial_condition(self):
        """
        Function restores initial conditions and displays it with opengl widget.
        """
        self.step = 0
        self.MBD_system._restore_initial_conditions()
        self.step_signal.signal_step.emit(self.step)
        self.energy_signal.signal_energy.emit(0., 0.)

        self.MBD_system.time = 0.
        self.MBD_system.step_num = 0

        self.refresh_signal.signal_refresh.emit()
        self.refresh(t=0, q=self.MBD_system.q0)


"""
Created on 13. mar. 2014

@author: luka.skrinjar
"""

import copy
import datetime
import inspect
import itertools
import logging
import os
import time

import numpy as np
from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import *

from DAE_fun import DAEfun
from MBD_system.check_filename import check_filename
from MBD_system.q2R_i import q2R_i
from MBD_system.q2dR_i import q2dR_i
from MBD_system.q2dtheta_i import q2dtheta_i
from MBD_system.solution_data.solution_data import SolutionData
from signals import ErrorTimeIntegrationSignal
from signals import SolutionSignal


class stepSignal(QtCore.QObject):
    signal_step = QtCore.pyqtSignal(int, name='signal_step')

class RepaintGLSignal(QtCore.QObject):
    signal_repaintGL = QtCore.pyqtSignal()

class FinishedSignal(QtCore.QObject):
    signal_finished = QtCore.pyqtSignal(str)

class SaveScreenshot(QtCore.QObject):
    signal_saveScreenshot = QtCore.pyqtSignal()

class solutionFilenameSignal(QtCore.QObject):
    signal_filename = QtCore.pyqtSignal(str, name='')

class EnergySignal(QtCore.QObject):
    signal_energy = QtCore.pyqtSignal(float, float, name='energy')


class SolveDynamicAnalysis(QObject):    # Thread, QObject
    """
    classdocs
    """
    __simulation_id = itertools.count(0)

    def __init__(self, MBD_system, parent=None):
        """
        Constructor
        """
        super(SolveDynamicAnalysis, self).__init__(parent)

        #   parent
        self._parent = parent

        #   states of simulations
        self.running = False
        self.stopped = False
        self.finished = False
        self.failed = False
        self._no_error = False

        #   signals
        self.step_signal = stepSignal()
        self.finished_signal = FinishedSignal()
        self.repaintGL_signal = RepaintGLSignal()
        self.energy_signal = EnergySignal()
        self.error_time_integration_signal = ErrorTimeIntegrationSignal()
        self.filename_signal = solutionFilenameSignal()
        self.save_screenshot_signal = SaveScreenshot()
        self.solution_signal = SolutionSignal()

        #   MBD system object
        self.MBD_system = MBD_system
        #    copy MBD object
        self._MBD_system = copy.copy(self.MBD_system)
        #    DAE fun object
        self.DAE_fun = DAEfun(self._MBD_system, parent=self)

        #    simulation settings - properties
        # self.integrationMethod = "RKF"  #    self.MBD_system.integrationMethod

        #    update opengl widget every N-th step
        self.step = 0
        self.update_opengl_widget_step_count = 0
        self.update_opengl_widget_every_Nth_step = self.MBD_system.updateEveryIthStep

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

    def _solution_containers(self):
        """

        :return:
        """
        #    solution container variable
        self.q0 = self._MBD_system.evaluate_q0()
        self.q_sol_matrix = np.array([self.q0])
        self.t_vector = np.array([[0]])
        self.R_vector = np.array([0])
        self.step_counter = np.array([0])
        self.step_size = np.array([self._MBD_system.Hmax])
        self.energy_data = np.array([self._mechanical_energy(self.q0)])

    def solve(self):
        """
        Solves system of ode with dormand-prince order 5 method with runge-kutta algorithm
        """
        #    copy MBD object
        self._MBD_system = copy.copy(self.MBD_system)
        #    ode fun object
        self.DAE_fun = DAEfun(self._MBD_system, parent=self)

        self.start_solve()

        self.start_time_simulation_info_in_sec_UTC = time.time()
        self.start_time_simulation_info = datetime.datetime.fromtimestamp(self.start_time_simulation_info_in_sec_UTC).strftime("%H:%M:%S %d.%m.%Y")
        print "Simulation started: ", self.start_time_simulation_info

        #    create mass and inverse mass matrix (only once)
        self.DAE_fun.evaluate_M()
        #    create array of initial conditions for differential equations
        if not self.MBD_system.q0_created:
            self.MBD_system.evaluate_q0()

        #   get initial conditions
        q0 = self._MBD_system.evaluate_q0()
        self.FLAG = 1

        #   solution containers
        self._solution_containers()

        self.simulation_id = 0
        #    0 - no contact
        #    1 - contact detected
        #    -1 - contact already happened - reduce integration step
        self.step = 0
        self.h_contact = self._MBD_system.Hmin

        if self.MBD_system.integrationMethod.title() == "ABAM":
            pass

        elif self.MBD_system.integrationMethod.title() == "Runge-Kutta":
            self.solve_ODE_RK(t_0=self._MBD_system.t_0, t_n=self._MBD_system.t_n, q0=q0, absTol=self._MBD_system.absTol, relTol=self._MBD_system.relTol, Hmax=self._MBD_system.Hmax, Hmin=self._MBD_system.Hmin)

        elif self.MBD_system.integrationMethod.title() == "Euler":
            self.solve_ODE_Euler(t_0=self._MBD_system.t_0, t_n=self._MBD_system.t_n, q0=q0, Hmax=self._MBD_system.Hmax)

        else:
            raise ValueError, "Selected integration method not supported! Method: %s"%self.MBD_system.integrationMethod

    def solve_ODE_Euler(self, t_0, t_n, q0, Hmax):
        """
        Euler method for time integration
        :param t_0:
        :param t_n:
        :param q0:
        :param h:
        :return:
        """
        self.t_0 = self._dt = t_0
        self.t_n = t_n
        t = t_0
        w = q0

        self.Hmax = h = Hmax
        self.R = 0

        while self.FLAG == 1:

            if self.stopped:
                self.updateGL(t=t, q=w)
                self.stop_solve()
                self.FLAG = 0

            if t >= self.t_n:
                self.FLAG = 0
                self.finished = True
                self.updateGL(t=t, q=w)

            self.t = t = t + h
            # print "-------------------------------------"
            # print "t =", t, "step =", self.step
            w = w + h * self.DAE_fun.evaluate_dq(h, t, w)
            # print "q =", w
            #    solve contacts here
            # time.sleep(100)
            self.FLAG_contact = self.__evaluate_contacts(t, w)
            # print "self.FLAG_contact =", self.FLAG_contact
            #   evaluate mechanical energy of system
            self._energy_t = self._mechanical_energy(w)
            #   evaluate difference of mechanical energy of a system between this and prevouis time step
            self._energy_delta = self._energy_t - self.energy_data[-1]

            # self.__evaluate_contacts(h, t, w)

            h, t, w = self._track_data(h, t, w)
            # print "track data"

            self._info(t, w)

            #   check time step
            h = self.check_time_step(h, w)

        if self.finished or self.failed:
            self.finished_solve()

        #    save data to file
        if self.MBD_system._save_options.lower() != "discard":
            self.write_simulation_solution_to_file()

    def solve_ODE_RK(self, t_0, t_n, q0, absTol, relTol, Hmax, Hmin):
        """
        RKF - Runge-Kutta-Fehlberg method for time integration
        Based on Numerical Analysis 9th ed Burden & Faires
        Args:
        t_0 -
        t_n -
        q0 -
        absTol -
        Hmax -
        Hmin -
        """
        self.t_0 = self._dt = t_0
        self.t_n = t_n
        t = t_0
        w = q0

#         logging.getLogger("DyS_logger").info("Size of MBD system in bytes %s" % sys.getsizeof(self.MBD_system))
        logging.getLogger("DyS_logger").info("Simulation started with initial conditions q0:\n%s" % q0)

        self.Hmax = Hmax
        self.Hmin = Hmin

        h = Hmax
        h_contact = self._MBD_system.Hmin
        self._t_FLAG1 = Hmax

#        print "absTol =", absTol
#        self.update_opengl_widget_every_Nth_step = 1*((t_n - t_0)/Hmax)
#        print "self.update_opengl_widget_every_Nth_step =", self.update_opengl_widget_every_Nth_step

        self.save_updated_GL_screenshot_to_file(self.t_0)
        np.set_printoptions(precision=20, threshold=1000, edgeitems=True, linewidth=1000, suppress=False, formatter={'float': '{: 10.9e}'.format})

        while self.FLAG == 1:
            # print "-------------------------------------"
            # print "step =", self.step, "contact =", self.FLAG_contact
            # print "t =", t,
            if self.stopped:
                self.updateGL(t=t, q=w)
                self.stop_solve()
                self.FLAG = 0

            if self.FLAG_contact == 1:
                h = h_contact
            elif self.FLAG_contact == 0:
                h = Hmax
            # print "h =", h,
            # print "self.FLAG_contact =", self.FLAG_contact
            # print "--------------------------------"
            # print "h =", h
            # print "t =", t
            # print "w =", w
            # print "self.DAE_fun.evaluate_dq(h, t, w) =", self.DAE_fun.evaluate_dq(h, t, w)
            K1 = h * self.DAE_fun.evaluate_dq(h, t, w)
            K2 = h * self.DAE_fun.evaluate_dq(h, t + (1 / 4.) * h, w + (1 / 4.) * K1)
            K3 = h * self.DAE_fun.evaluate_dq(h, t + (3 / 8.) * h, w + (3 / 32.) * K1 + (9 / 32.) * K2)
            K4 = h * self.DAE_fun.evaluate_dq(h, t + (12 / 13.) * h, w + (1932 / 2197.) * K1 - (7200 / 2197.) * K2 + (7296 / 2197.) * K3)
            K5 = h * self.DAE_fun.evaluate_dq(h, t + h, w + (439 / 216.) * K1 - 8 * K2 + (3680 / 513.) * K3 - (845 / 4104.) * K4)
            K6 = h * self.DAE_fun.evaluate_dq(h, t + (1 / 2.) * h, w - (8 / 27.) * K1 + 2 * K2 - (3544 / 2565.) * K3 + (1859 / 4104.) * K4 - (11 / 40.) * K5)

            self.R = (1. / h) * np.linalg.norm((1 / 360.) * K1 - (128 / 4275.) * K3 - (2197 / 75240.) * K4 + (1 / 50.) * K5 + (2 / 55.) * K6)
            self.R = absTol
            #    if calculated difference is less the absTolerance limit accept the calculated solution
            if self.R <= absTol or self.FLAG_contact == 1:
                self.t = t = t + h
                # print "i =", self.step, "t =", t
                # self._parent.ui.simulation_progressBar.setValue(int(1/self.t))
                # print "t =", t, "h =", h,
                w = w + (25. / 216.) * K1 + (1408. / 2565.) * K3 + (2197. / 4104.) * K4 - (1. / 5.) * K5

#                 if self.step > 220:
#                     print self.step, w
#                     print "-----------------------------------"
                #    check time step size
#                t = self.check_time_step(t, w)
#                print "t_out =", t
#                 print "t = %1.3e" % t  # , "absTol =%10.3E" % R_, "h =%10.3E" % h, "step =", self.step
#                 __print_options = np.get_printoptions()
#                 np.set_printoptions(precision=10, threshold=1000, edgeitems=True, linewidth=1000, suppress=False, formatter={'float': '{: 10.9e}'.format})
#                 np.set_printoptions(**__print_options)

                #    solve contacts here
                # print "t =", t
                # time.sleep(100)
                self.FLAG_contact = self.__evaluate_contacts(t, w)
                # print "self.FLAG_contact =", self.FLAG_contact
                #   evaluate mechanical energy of system
                self._energy_t = self._mechanical_energy(w)
                self._energy_delta = self._energy_t - self.energy_data[-1]

                h, t, w = self._track_data(h, t, w)

                #   evaluate difference of mechanical energy of a system between this and prevouis time step

                self._info(t, w)

                #   check time step
                h = self.check_time_step(h, w)

            #    evaluate error and change step size if needed
            # delta = 0.84 * (absTol / R) ** (1 / 4.)
            #
            #
            # if self.FLAG_contact == 0 or self.FLAG_contact == 1:
            #     if delta <= 0.1:
            #         h = 0.1 * h
            #     elif delta >= 4:
            #         h = 4 * h
            #     else:
            #         h = delta * h
            #
            #     if h > Hmax:
            #         h = Hmax
            #
            # else:
            #     #   new step size - h has to be smaller than the last calculated over impact time
            #     # print "t =", t
            #     # print "h =", h
            #     # print "t+h =", t+h
            #     # print "self._t_FLAG1 =", self._t_FLAG1
            #     if t + h >= self._t_FLAG1:
            #         h = (self._t_FLAG1 - t)/2
            #     # print "h =", h
            #
            #   if end time is reached stop/break integration process
            if t >= self.t_n:
                self.FLAG = 0
                self.finished = True
                self.updateGL(t=t, q=w)
            #
            # #    reduce step size to get to final time t_n
            # elif t + h > t_n:
            #     h = t_n - t
            #
            # #    break integration as minimum step size is exceeded
            # elif h < Hmin and R > absTol:
            #     self.failed = True
            #     self.FLAG = 0
            #     print "t =", t
            #     print "absTol =", absTol
            #     print "R =", R
            #     print "h =", h
            #     h = Hmin
            #     self.error_time_integration_signal.signal_time_integration_error.emit()
            #
            #     print "Hmin exceeded! Procedure failed!"
            #    time step after contact is constant and very small value
            if self.FLAG_contact == 1:
                h = self.Hmin
                # if h_after_contact < 10 * Hmin:
                #     h = 10 * Hmin
                # else:
                #     h = Hmax

            #    integration finished - end time reached successfully
            if self.finished or self.failed or self.stopped:
                self.finished_solve()

        #    save data to file
        self.write_simulation_solution_to_file()

    def _info(self, t, w):
        """

        :return:
        """
        #    update coordinates and angles of all bodies
        self.__update_coordinates_and_angles_of_all_bodies(w)

        #    update opengl display
        self.updateGL(t, w)

        # self.step += 1
        # step_counter = np.append(step_counter, self.step)

        self.update_opengl_widget_step_count += 1
        self.save_screenshot_step_count += 1

        #    step number signal
        self.step_signal.signal_step.emit(self.step)

        #    energy data signal
        self.energy_signal.signal_energy.emit(self._energy_t, self._energy_delta)

    def _track_data(self, h, t, q):
        """

        :return:
        """
        #    save values at next time step to matrix of solution (q_matrix)
        if self.FLAG_contact == 0 or self.FLAG_contact == 1:
            #   append solution at i-th time to solution matrix
            self.q_sol = self.q_sol_matrix = np.vstack((self.q_sol_matrix, np.array(q)))

            self.step_size = np.vstack((self.step_size, h))
            self.step += 1
            step_num = self.step

            #   append current step number to step counter history array
            self.step_counter = np.append(self.step_counter, self.step)
            #    add to vector
            self.t_vector = np.vstack((self.t_vector, t))
            self.R_vector = np.vstack((self.R_vector, self.R))
            #   append current mechanical energy to array
            self.energy_data = np.append(self.energy_data, self._energy_t)

            for contact in self.MBD_system.contacts:
                contact._track_data(self.step, h, t, q)

        #    reduce step size and continue integration from previous step
        elif self.FLAG_contact == -1:
            #    go to one before last solution
            q = self.q_sol_matrix[-1, :]
            #    store time that is too large
            self._t_FLAG1 = t
            t = self.t_vector[-1, 0]  # , 0

            #    reduce step size
            h = 0.5 * h
            # self.step -= 1
            # step_num  = self.step
            # step_counter[-1] = self.step

            self.energy_data[-1] = self._energy_t

        return h, t, q

    def __evaluate_contacts(self, t, q):
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
        # print "---------------------------------"
        # print "t =", t, "step =", self.step,
        self.__update_coordinates_and_angles_of_all_bodies(q)

        self.contact_status_list = []
        #    check for contact of every contact pair
        for contact in self.MBD_system.contacts:
            contact.data_tracker(t, self.step)
            # print "contact._contact_point_found =", contact._contact_point_found, "t =", t
            if not contact._contact_point_found:
                #    function recursively loops through all potential AABB objects in AABB tree of each body
                #    and if two AABB objects overlap, a new overlap object is created
                #    that has both overlapping AABB objects
                #   general contact
                if contact._type == "General":
                    #   adds simulation data to contact object
                    #   reset to empty list a contact object attribute of list of overlap pairs,
                    #   as this is not possible in the next line (function) due to the recursion
                    contact.AABB_list_of_overlap_pairs = []
                    #   check for overlap and if overlap is present build a overlap pair object
                    contact.check_for_overlap(q, contact.AABB_i, contact.AABB_j)

                    #   if list of overlap object pairs is not empty, check for contact
                    if contact.AABB_list_of_overlap_pairs is not []:
                        status = contact.check_for_contact(q)
                        self.contact_status_list.append(status)
                    else:
                        contact.no_overlap()

                #   revolute clearance joint contact
                elif contact._type.lower() == "revolute clearance joint" or contact._type.lower() == "contact sphere-sphere" or contact._type.lower() == "contact plane-sphere" or contact._type.lower() == "pin-slot clearance joint linear":
                    status = contact.check_for_contact(q)
                    # print "t =", t, "status(check_for_contact) =", status
                    self.contact_status_list.append(status)

                else:
                    self.FLAG = 0
                    QtGui.QMessageBox.critical(self._parent, "Error!", "Contact type not correct: %s"%contact._type+". \nDefine correct contact type!", QtGui.QMessageBox.Ok)
                    raise ValueError, "Contact type not correct: %s"%contact._type, "\n Define correct contact type!"
            else:
                #    adds simulation data to contact object as this is not possible in the next line because of the recursion
                # print "contact.update_status()"
                # print "contact_update() @ solve_ODE"
                status = contact.contact_update(self.step, t, q)
                # print "status (from contact_update) =", status
                self.contact_status_list.append(status)

            self.contact_status_list = np.array(self.contact_status_list)
            # print "self.contact_status_list =", self.contact_status_list

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
                return FLAG_contact

            #    contact
            if (self.contact_status_list == 1).any():
                FLAG_contact = 1
                # logging.getLogger("DyS_logger").info("Contact detected at simulation time:%s \nbody_i=%s \nbody_j=%s", t, self.MBD_system.bodies[contact.body_id_i]._name, self.MBD_system.bodies[contact.body_id_j]._name)
                # logging.getLogger("DyS_logger").info("Contact detected at simulation time:%s "%t)
                #   repaint opengl widget

#                 self.repaintGL_signal.signal_repaintGL.emit()

                #   solve contacts - construct contact forces
                # self.DAE_fun.solve_contacts(t, q)

                # logging.getLogger("DyS_logger").info("Contact calculation finished")
                return FLAG_contact

            #    contact has already happened - reduce time step
            if (self.contact_status_list == -1).any():
                FLAG_contact = -1  # 0-no contact, 1-contact has already happened - reduce time step
                return FLAG_contact

        else:
            FLAG_contact = 0
            return FLAG_contact

    def _mechanical_energy(self, q):
        """

        """
        #    predefine zero array
        _energy = np.zeros(self.MBD_system.number_of_bodies)

        #   energy of all bodies
        for i, body in enumerate(self.MBD_system.bodies):
            _q = q2R_i(q, body.body_id)
            _dq = np.append(q2dR_i(q, body.body_id), q2dtheta_i(q, body.body_id))

            #   body energy
            # body_energy = 0.5 * (body.mass * (_dq**2) + body.J_zz * (_omega**2)) + (body.mass * self.MBD_system.gravity * _q[1])
            body_energy = body.mechanical_energy(q=_q, dq=_dq, gravity=self.MBD_system.gravity)
            _energy[i] = body_energy

        #   energy of normal contact forces
        _energy_of_contacts = np.zeros(len(self.MBD_system.contacts))
        # for i, contact in enumerate(self.MBD_system.contacts):
        #     _energy_of_contacts[i] = contact.mechanical_energy()

        #   total mechanical energy
        energy = np.sum(_energy) + np.sum(_energy_of_contacts)

        return energy

    def check_time_step(self, h, q):
        """

        """
        # _t_max = []
        # for contact in self.MBD_system.contacts:
        #     if contact.AABB_list_of_overlap_pairs != []:
        #         for body_id in contact.contact_body_id_list:
        #
        #             _t = self.MBD_system.bodies[body_id].skin_thickness / np.linalg.norm(q2dR_i(q, body_id) + q2dtheta_i(q, body_id) * self.MBD_system.bodies[body_id].uP_i_max, ord=2)
        #             _t_max.append(_t)
        #
        #         t_max = min(_t_max)
        #         if t > t_max:
        #             return t_max
        #         else:
        #             return t
        #     else:
        #         return t
        # if h < self.Hmin:
        #     h = self.Hmin
        if self.FLAG_contact == 1:
            h = self.h_contact
        elif self.FLAG_contact == -1:
            h = h
        else:
            h = self.Hmax

        return h

    def write_simulation_solution_to_file(self):
        """
        Function writes solution data to solution data object
        """
        print "write_simulation_solution_to_file()"
        solution_data = self.collect_solution_data()

        #   write solution data and info to solution data object
        self._solution_data.set_filetype(self.MBD_system._solution_filetype)

        self.solution_filename = 'solution_data_%02d'%self.simulation_id#+self.MBD_system._solution_filetype

        if self.MBD_system._save_options == "save to new":
            #   check filename if already exists
            self.solution_filename = check_filename(self.solution_filename)

        #   set solution data filename
        self._solution_data.setName(self.solution_filename)

        #   add data to object attribute
        self._solution_data.add_data(solution_data)
        self._solution_data.finished = True

        # print "exe1"
        # print "self.MBD_system._solution_save_options =", self.MBD_system._solution_save_options
        if self.MBD_system._solution_save_options != "discard":
            print "exe2"
            #   write data to file
            self._solution_data.write_to_file()

            #   emit signal
            # self.solution_signal.solution_data.emit(id(self._solution_data), self._solution_data._filename)

            #   write contact values of every simulation time step to file
            for contact in self.MBD_system.contacts:
                if contact._solution_save_options != "discard":
                    contact.write_to_file()

    def collect_solution_data(self):
        """

        :return:
        """
        #   reshape data to join in one matrix
        #   step
        _step = np.array([self.step_counter]).T
        #   energy
        _energy_data = np.array([self.energy_data]).T
        #   add all solution data to one variable
        solution_data = np.hstack((_step, _energy_data, self.R_vector, self.step_size, self.t_vector, self.q_sol_matrix))

        return solution_data

    def save_solution_data(self):
        """

        :return:
        """
        data = self.collect_solution_data()
        self._solution_data.add_data(data)

    def start_solve(self):
        """
        Start solver (time integration process)
        """
        self.simulation_id = self.__simulation_id.next()

        #   create solution data object to store data at each simulation start
        self._solution_data = SolutionData(MBD_system=self.MBD_system)#parent=self.MBD_system.Solution
        print "self._solution_data._name =", self._solution_data._name
        self.MBD_system.solutions.append(self._solution_data)

        logging.getLogger("DyS_logger").info("Simulation started, simulation id: %s"%self.simulation_id)

        self.running = True
        self.stopped = False
        self.finished = False

    def stop_solve(self):
        """
        Stop solver (time integration process)
        """
        self.stopped = True
        self.running = False

        self.stop_time_simulation_info_in_sec_UTC = time.time()
        self.stop_time_simulation_info = datetime.datetime.fromtimestamp(self.stop_time_simulation_info_in_sec_UTC).strftime("%H:%M:%S %d.%m.%Y")

        logging.getLogger("DyS_logger").info("Simulation stopped by user at simulation time:%s", self.t)

    def finished_solve(self):
        """
        Method sets some parameters when integration process is finished
        :return:
        """
        if self.finished:
            self.save_solution_data()
            self.solution_signal.solution_data.emit(id(self._solution_data), self._solution_data._name)
            logging.getLogger("DyS_logger").info("Simulation finished successfully!")
            self.finished_signal.signal_finished.emit("Finished")
        else:
            logging.getLogger("DyS_logger").info("Simulation failed!")

        #   set booleans
        self.finished = True
        self.running = False

    def load_simulation_solution_from_file(self, filename):
        """
        Function load solution_data from file.
        """
        self.finished_signal.signal_finished.emit("Animation")
        if os.path.exists(filename):
            self.solution_data = np.loadtxt(str(filename), skiprows=2)
            return self.solution_data

    def __update_coordinates_and_angles_of_all_bodies(self, q):
        """
        Function updates MBD_system data before display
        """
        self.MBD_system.update_coordinates_and_angles_of_all_bodies(q)

    def updateGL(self, t, q):
        """
        Update opengl widget at preset time steps
        """
        #   update display on every i-th step
        if self._parent._update_display_type == "step":
            if self.update_opengl_widget_step_count == self.update_opengl_widget_every_Nth_step or self.finished or self.stopped:
                self._updateGL()

                self.time = t

                self.save_updated_GL_screenshot_to_file(self.time)

                #    reset step counter to update opengl widget
                self.update_opengl_widget_step_count = 0

        #   update display on simulation time interval dt
        if self._parent._update_display_type == "dt":
            _t = t - self._dt
            if _t >= self._parent._dt:
                self._updateGL(t, self.step)

                #   set new value
                self._dt = t


        return None

    def _updateGL(self, t, step):
        """

        :return:
        """
        #    update opengl widget
        self.MBD_system.update_simulation_properties(time=t, step_num=step)
        #   emit signal
        self.repaintGL_signal.signal_repaintGL.emit()

    def save_updated_GL_screenshot_to_file(self, t):
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
        print "restore_initial_condition()"
        self.step = 0
        self.MBD_system._restore_initial_conditions()

        self.repaintGL_signal.signal_repaintGL.emit()
        self.updateGL(t=0, q=self.MBD_system.q0)
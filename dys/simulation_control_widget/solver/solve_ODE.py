'''
Created on 13. mar. 2014

@author: luka.skrinjar
'''
import operator
import copy
import logging
import datetime
import os
from pprint import pprint
import sys
import time
import numpy as np
import itertools
import matplotlib.pyplot as plt
from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import *
from xlwt import Workbook
from dxfwrite import DXFEngine as dxf


from ODE_fun import ODEfun
from MBD_system.q2dR_i import q2dR_i
from MBD_system.q2dtheta_i import q2dtheta_i
from MBD_system.q2R_i import q2R_i
from MBD_system import convert_bytes_to_
from signals import ErrorTimeIntegrationSignal

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
    

class SolveODE(QObject):  # Thread, QObject
    """
    Integrator
    """
    __simulation_id = itertools.count(0)
    
    def __init__(self, MBD_system=[], parent=None):
        super(SolveODE, self).__init__(parent)

        self._parent = parent

        self.running = False
        self.stopped = False
        self.finished = False
        self.failed = False
        self.step_signal = stepSignal()
        self.finished_signal = FinishedSignal()
        self.repaintGL_signal = RepaintGLSignal()
        self.energy_signal = EnergySignal()
        self.error_time_integration_signal = ErrorTimeIntegrationSignal()
        
        self.filename_signal = solutionFilenameSignal()
        
        self.save_screenshot_signal = SaveScreenshot()

        
        self.MBD_system = MBD_system
        #    copy MBD objec
        self._MBD_system = copy.copy(self.MBD_system)
        #    ode fun object
        self.ode_fun = ODEfun(MBD_system=self._MBD_system, parent=self)


        self.FLAG_contact = 0
        self.simulation_id = 0

        #    simulation settings - properties
        self.integrationMethod = "RKF"  #    self.MBD_system.integrationMethod
        
        
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
    
        
    def solve_ODE(self):
        """
        Solves system of ode with dormand-prince order 5 method with runge-kutta algorithm
        in:
            ode_fun - ode function
            q0 - initial values (conditions) array q0_i = [Rx_i, Ry_i, theta_i, dRx_i, dRy_i, dtheta_i]
        """
        
        #    copy MBD objec
        self._MBD_system = copy.copy(self.MBD_system)
        #    ode fun object
        self.ode_fun = ODEfun(MBD_system=self._MBD_system, parent=self)
        
        self.start_solve()
        
        self.start_time_simulation_info_in_sec_UTC = time.time()
        self.start_time_simulation_info = datetime.datetime.fromtimestamp(self.start_time_simulation_info_in_sec_UTC).strftime("%H:%M:%S %d.%m.%Y")
        print "Simulation started: ", self.start_time_simulation_info
        
        #    create mass and inverse mass matrix (only once)
        self.ode_fun.create_M()
        #    create array of initial conditions for differential equations
        if not self.MBD_system.q0_created:
            self.MBD_system.create_q0()
            
        q0 = self.initial_conditions()

        if self.integrationMethod == "RKF":

            self.solve_ODE_RKF(t_0=self._MBD_system.t_0, t_n=self._MBD_system.t_n, q0=q0, absTol=self._MBD_system.absTol, relTol=self._MBD_system.relTol, Hmax=self._MBD_system.Hmax, Hmin=self._MBD_system.Hmin)

        elif self.integrationMethod == "ABAM":
            print "SELECTED METHOD: ABAM" 
        elif self.integrationMethod == "Newton-Rhapson":
            print "SELECTED METHOD: Newton-Rhapson" 
        else:
            None


    def solve_ODE_RKF(self, t_0, t_n, q0, absTol, relTol, Hmax, Hmin):
        """
        RKF - Runge-Kutta-Fehlberg 
        Based on Numerical Analysis 9th ed Burden & Faires
        Args:
        t_0 - 
        t_n - 
        q0 - 
        absTol - 
        Hmax - 
        Hmin - 
        """
        self.t_0 = t_0
        self.t_n = t_n
        t = t_0
        w = q0

#         logging.getLogger("DyS_logger").info("Size of MBD system in bytes %s" % sys.getsizeof(self.MBD_system))
        logging.getLogger("DyS_logger").info("Simulation started with initial conditions q0:\n%s" % q0)


        h = Hmax
        h_contact = 1E-5
        self._t_FLAG1 = Hmax

        self.FLAG = 1
        self.FLAG_contact = 0  #    0 - no contact
                                #    1 - contact detected
                                #    -1 - contact already happened - reduce integration step
        
        #    solution container variable
        q_sol_matrix = np.array([w])
        t_vector = np.array([[0]])
        R_vector = np.array([0])
        step_counter = np.array([0])
        step_size = np.array([h])
        energy_data = np.array([self.__mechanical_energy(q0)])
        
        
#        print "absTol =", absTol
#        self.update_opengl_widget_every_Nth_step = 1*((t_n - t_0)/Hmax)
#        print "self.update_opengl_widget_every_Nth_step =", self.update_opengl_widget_every_Nth_step

        self.step = step_num = 0
        self.save_updated_GL_screenshot_to_file(self.t_0)
        np.set_printoptions(precision=20, threshold=1000, edgeitems=True, linewidth=1000, suppress=False, formatter={'float': '{: 10.9e}'.format})


        while self.FLAG == 1:
            # print "h =", h
            if self.stopped:
                self.update_GL_(t=t, q=w)
                self.stop_solve()
                break

            if self.FLAG_contact == 1:
                h = h_contact

            # print "self.FLAG_contact =", self.FLAG_contact
            print "--------------------------------"
            print "i =", self.step, "h =", h, "t =", t+h,
            # print "t =", t

            K1 = h * self.ode_fun.create_dq(h, t, w)
            K2 = h * self.ode_fun.create_dq(h, t + (1 / 4.) * h, w + (1 / 4.) * K1)
            K3 = h * self.ode_fun.create_dq(h, t + (3 / 8.) * h, w + (3 / 32.) * K1 + (9 / 32.) * K2)
            K4 = h * self.ode_fun.create_dq(h, t + (12 / 13.) * h, w + (1932 / 2197.) * K1 - (7200 / 2197.) * K2 + (7296 / 2197.) * K3)
            K5 = h * self.ode_fun.create_dq(h, t + h, w + (439 / 216.) * K1 - 8 * K2 + (3680 / 513.) * K3 - (845 / 4104.) * K4)
            K6 = h * self.ode_fun.create_dq(h, t + (1 / 2.) * h, w - (8 / 27.) * K1 + 2 * K2 - (3544 / 2565.) * K3 - (1859 / 4104.) * K4 - (11 / 40.) * K5)
            
            R = (1 / h) * np.linalg.norm((1 / 360.) * K1 - (128 / 4275.) * K3 - (2197 / 75240.) * K4 + (1 / 50.) * K5 + (2 / 55.) * K6)
            R = absTol
            # R_ = np.linalg.norm(R, ord=2)
            # print "R =", R_
            #    if calculated difference is less the absTolerance limit accept the calculated solution 
            if R <= absTol or self.FLAG_contact == 1:
                self.t = t = t + h
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
#                 print "t(solve_ODE.py) =", t
                w = self.solve_contacts(t, w)

                #   evaluate mechanical energy of system
                _energy = self.__mechanical_energy(w)

                #   evaluate difference of mechanical energy of a system between this and prevouis time step
                _energy_delta = _energy - energy_data[-1]
                
                #    save values at next time step to matrix of solution (q_matrix)
                if self.FLAG_contact == 0 or self.FLAG_contact == 1:
#                     print self.step, t, R_,

                    self.q_sol = q_sol_matrix = np.vstack((q_sol_matrix, np.array(w)))
                    
                    step_size = np.vstack((step_size, h))
                    self.step += 1
                    step_num = self.step
                    # print "i =", step_num,


                    step_counter = np.append(step_counter, self.step)
                    #    add to vector
                    self.t_sol = t_vector = np.vstack((t_vector, t))
                    R_vector = np.vstack((R_vector, R))
                    
                    
                    energy_data = np.append(energy_data, _energy)
                    
                #    reduce step size and continue integration from previous step
                elif self.FLAG_contact == -1:
                    #    go to one before last solution
                    w = q_sol_matrix[-1, :]
                    #    store time that is too large
                    self._t_FLAG1 = t
                    t = t_vector[-1, 0]  # , 0
                    
                    #    reduce step size
                    h = 0.5 * h
                    # self.step -= 1
                    # step_num  = self.step
                    # step_counter[-1] = self.step
                    
                    energy_data[-1] = _energy
                    
                #    update coordinates and angles of all bodies
                self.__update_coordinates_and_angles_of_all_bodies(w)

                #    update opengl display
                self.update_GL_(t, w)

                # self.step += 1
                # step_counter = np.append(step_counter, self.step)
                
                self.update_opengl_widget_step_count += 1
                self.save_screenshot_step_count += 1

                #    step number signal
                self.step_signal.signal_step.emit(self.step)
                
                #    energy data signal
                self.energy_signal.signal_energy.emit(_energy, _energy_delta)

            
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
            #    if end time is reached stop/break integration process
            if t >= self.t_n:
                self.FLAG = 0
                self.finished = True
                self.update_GL_(t=t, q=w)

            #    reduce step size to get to final time t_n
            elif t + h > t_n:
                h = t_n - t
            
            #    break integration as minimum step size is exceeded
            elif h < Hmin and R > absTol:
                self.failed = True
                self.FLAG = 0
                print "t =", t
                print "absTol =", absTol
                print "R =", R
                print "h =", h
                h = Hmin
                self.error_time_integration_signal.signal_time_integration_error.emit()

                print "Hmin exceeded! Procedure failed!"
            #    time step after contact is constant and very small value
            if self.FLAG_contact == 1:
                h = 1E-6
                # if h_after_contact < 10 * Hmin:
                #     h = 10 * Hmin
                # else:
                #     h = Hmax
                
        
        #    integration finished - end time reached successfully
        if self.finished or self.failed:
            self.finished_solve()

        self.R_vector = R_vector
        self.t_sol = t_vector
        self.q_sol = q_sol_matrix
        self.step_size = step_size
        self.step = step_counter
        self.energy_data = energy_data


        #    save data to file
        self.write_simulation_solution_to_file()


    def solve_contacts(self, t, q):
        """
        Function solves contacts. Finds overlap pairs of (sub)AABB and calculates actual distances between node 
        and line (edge). If distance is below userdefined limit contact is present and contact equations are solved.
        :type self: object
        Args:
            q    vector of positions and velocities (translational and rotational)
        Returns:
            q_   vector of new calculated velocities (positions are equal before and after contact
            status - 0 - contact is not detected continue with integration to next time step
                     1 - contact is detected (has happened) return to previous time step solution and split time step size
                     2 - contact is detected (will happen) split time step size   
        Raises:
            None 
        """
        self.__update_coordinates_and_angles_of_all_bodies(q)

        self.contact_status_list = []
        # print t, self.step
        #    check for contacts
        if self.MBD_system.contacts != []:
            
            #    check for contact of every contact pair
            for contact in self.MBD_system.contacts:
                contact.data_tracker(t, self.step)
                if not contact._contact_point_found:
                    #    function recursively loops through all potential AABB objects in AABB tree of each body
                    #    and if two AABB objects overlap, a new overlap object is created
                    #    that has both overlapping AABB objects

                    #    adds simulation data to contact object as this is not possible in the next line because of the recursion
#                     __step_num = self.step+0

                    # contact.data_tracker(t, self.step)

                    #   general contact
                    if contact.AABB_i is not None and contact.AABB_j is not None:
                        # #    adds simulation data to contact object as this is not possible in the next line because of the recursion
                        # __step_num = self.step+0
                        # contact.data_tracker(t, __step_num)

                        contact.AABB_list_of_overlap_pairs = []
                        contact.check_for_AABB_AABB_overlap(q, contact.AABB_i, contact.AABB_j)

                        if contact.AABB_list_of_overlap_pairs != []:
                            status = contact.check_for_contact(q)
                            self.contact_status_list.append(status)
                        else:
                            contact.no_overlap()

                    #   revolute clearance joint contact
                    elif contact.u_iP is not None and contact.u_jP is not None:
                        # print "contact.check_for_contact()"
                        status = contact.check_for_contact(q)
                        self.contact_status_list.append(status)
                        # print "status =", status

                else:
                    #    adds simulation data to contact object as this is not possible in the next line because of the recursion
                    # print "contact.update_status()"
                    status = contact.contact_update(self.step, t, q)
                    # print "status (from contact_update) =", status
                    self.contact_status_list.append(status)

            self.contact_status_list = np.array(self.contact_status_list)
            print "t =", t, "self.contact_status_list =", self.contact_status_list

            #    no contacts
            if (self.contact_status_list == 0).all():
                self.FLAG_contact = 0
                return q
            
            #    contact
            if (self.contact_status_list == 1).any():
                self.FLAG_contact = 1
                # logging.getLogger("DyS_logger").info("Contact detected at simulation time:%s \nbody_i=%s \nbody_j=%s", t, self.MBD_system.bodies[contact.body_id_i]._name, self.MBD_system.bodies[contact.body_id_j]._name)
                # logging.getLogger("DyS_logger").info("Contact detected at simulation time:%s "%t)
                #   repaint opengl widget

#                 self.repaintGL_signal.signal_repaintGL.emit()

                #   solve contacts - construct contact forces
                self.ode_fun.solve_contacts(t, q)

                # logging.getLogger("DyS_logger").info("Contact calculation finished")
                return q
            
            #    contact has already happened - reduce time step
            if (self.contact_status_list == -1).any():
                self.FLAG_contact = -1  # 0-no contact, 1-contact has already happened - reduce time step
                return q
        
        else:
            self.FLAG_contact = 0
            return q


    def __mechanical_energy(self, q):
        """
        
        """
        #    predefine zero array
        _energy = np.zeros(self.MBD_system.number_of_bodies)
        

        for i, body in enumerate(self.MBD_system.bodies):

            _q = q2R_i(q, body.body_id)
            _dq = np.linalg.norm(q2dR_i(q, body.body_id), ord=2)
            _omega = q2dtheta_i(q, body.body_id)
            
            body_energy = 0.5 * (body.mass * (_dq**2) + body.J_zz * (_omega**2)) + (body.mass * self.MBD_system.gravity * _q[1]) 

            _energy[i] = body_energy
        
        energy = np.sum(_energy)
        
        return energy
 

    def check_time_step(self, t, q):
        """
        
        """
        _t_max = []
        for contact in self.MBD_system.contacts:
            if contact.AABB_list_of_overlap_pairs != []:
                for body_id in contact.contact_body_id_list:
                     
                    _t = self.MBD_system.bodies[body_id].skin_thickness / np.linalg.norm(q2dR_i(q, body_id) + q2dtheta_i(q, body_id) * self.MBD_system.bodies[body_id].uP_i_max, ord=2)
                    _t_max.append(_t)
        
                t_max = min(_t_max)
                if t > t_max:
                    return t_max
                else:
                    return t
            else:
                return t


    def write_simulation_solution_to_file(self):
        """
            Function writes solution_data to txt file and saves it.
        """
        _step = np.array([self.step]).T

        _energy_data = np.array([self.energy_data]).T


        self.solution_data = np.hstack((_step, _energy_data, self.R_vector, self.step_size, self.t_sol, self.q_sol))


        print "self.MBD_system._solution_filetype =", self.MBD_system._solution_filetype
        print "self.simulation_id =", self.simulation_id
        self.solution_filename = 'solution_data_%02d'%self.simulation_id+self.MBD_system._solution_filetype
        #    order of columns: step number, energy, time, q
        __frmt = ['%5i']+['%20.16f']+['%20.16f']+['%20.16f']+['%20.16f']+['%.10E']*len(self.initial_conditions())

        __header = "i-th step mechanical  energy  \t  R \t\t\t\t\t  dt \t\t\t\t\t  time \t"

        #    add header for q
        for body in sorted(self.MBD_system.bodies, key=lambda body: body.body_id):
            _id = body.body_id

            _header = "\t\t\t\tRx_"+str(_id) +"\t\t\t\tRy_"+str(_id) +"\t\t\t\ttheta_"+str(_id)
            __header = __header + _header

        #    add header for dq
        for body in sorted(self.MBD_system.bodies, key=lambda body: body.body_id):
            _id = body.body_id

            _header = "\t\t\t\tdRx_"+str(_id) +"\t\t\t\tdRy_"+str(_id) +"\t\t\t\tomega_"+str(_id)
            __header = __header + _header

        __comments ='#Insert comments here.\n'

        np.savetxt(self.solution_filename, self.solution_data, fmt=__frmt, delimiter='\t', header = __header, comments = __comments)


        logging.getLogger("DyS_logger").info("Solution data saved to file: %s. Size is %s", self.solution_filename, convert_bytes_to_.convert_size(os.path.getsize(self.solution_filename)))

        self.filename_signal.signal_filename.emit(self.solution_filename)


        # for contact in self.MBD_system.contacts:
        #     contact.save_solution_data()


    def start_solve(self):
        """

        """
        self.simulation_id = self.__simulation_id.next()
        logging.getLogger("DyS_logger").info("Simulation started, simulation id: %s"%self.simulation_id)
        
        self.running = True
        self.stopped = False
        self.finished = False


    def stop_solve(self):
        """

        """
        self.stopped = True
        self.running = False

        self.stop_time_simulation_info_in_sec_UTC = time.time()
        self.stop_time_simulation_info = datetime.datetime.fromtimestamp(self.stop_time_simulation_info_in_sec_UTC).strftime("%H:%M:%S %d.%m.%Y")

        logging.getLogger("DyS_logger").info("Simulation stopped by user at simulation time:%s", self.t)


    def finished_solve(self):
        self.finished = True
        self.running = False

        self.finished_signal.signal_finished.emit("Finished")
        
        logging.getLogger("DyS_logger").info("Simulation finished successfully!")

    
    def load_simulation_solution_from_file(self, filename):
        '''
        Function load solution_data from file.
        '''
        self.finished_signal.signal_finished.emit("Animation")
        if os.path.exists(filename):
            self.solution_data = np.loadtxt(str(filename), skiprows=2)
            return self.solution_data
    
    
    def __update_coordinates_and_angles_of_all_bodies(self, q):
        '''
        Function updates MBD_system data before display
        '''
        self.MBD_system.update_coordinates_and_angles_of_all_bodies(q)
    
    
    def update_GL_(self, t, q):
        """
        Update opengl widget at preset time steps
        """
        
        if self.update_opengl_widget_step_count == self.update_opengl_widget_every_Nth_step or self.finished or self.stopped:
            #    update opengl widget
            self.MBD_system.update_simulation_properties(time=t, step_num=self.step)

            self.repaintGL_signal.signal_repaintGL.emit()

            self.time = t

            self.save_updated_GL_screenshot_to_file(self.time)
            
            #    reset step counter to update opengl widget
            
            self.update_opengl_widget_step_count = 0


        return None
    
    
    def save_updated_GL_screenshot_to_file(self, t):
        """
        
        """
        if self.MBD_system.save_screenshots:
            if (self.save_screenshot_step_count == self.save_screenshot_every_Nth_step) or (t == self.t_0) or (self.finished):

                self.screenshot_filename_abs_path = os.path.join(self.MBD_system.saved_screenshots_folder_abs_path, str("step=" + "%010d" % self.step))
                self.save_screenshot_signal.signal_saveScreenshot.emit()
                self.save_screenshot_step_count = 0


    def initial_conditions(self):
        #    create array of initial conditions for differential equations
        self.q0 = self._MBD_system.create_q0()
        return self.q0
    
    
    def restore_initial_condition(self):
        """
        Function restores initial conditions and displays it with opengl widget.
        """
        self.step = 0
        self.MBD_system._restore_initial_conditions()
        
        self.repaintGL_signal.signal_repaintGL.emit()
        self.update_GL_(t=0, q=self.MBD_system.q0)

"""
Created on 9. jul. 2014

@author: lskrinjar
"""
import sys
import os
import itertools
import numpy as np
from matplotlib import pyplot as plt
from pprint import pprint
import subprocess
from moviepy.editor import *

from MBD_system.MBD_system import *
from MBD_system.MBD_system_items import SolutionDataItem


class SolutionData(SolutionDataItem):
    """
    classdocs
    """
    __id = itertools.count(0)

    def __init__(self, name=None, _file=None, parent=None):
        super(SolutionData, self).__init__(name, parent)
        """
        Constructor of solution data class
        """
        self._parent = parent

        #   solution id
        self._id = self.__id.next()

        if name is None:
            self._name = "solution_data_%.2i"%self._id
        else:
            self._name = name

        #   abs path to solution file
        self._file_abs_path = _file

        #    set file type: .dat, .xlsx, .csv, .sol
        self._filetype = ".sol"

        #   load data from file if already exists
        self.solution_data = None
        if os.path.isfile(str(self._file_abs_path)):
            _path, _file = os.path.split(self._file_abs_path)
            self.filename, self._filetype = os.path.splitext(_file)
            self.solution_data = self.read_file(self._file_abs_path)

        #    simulation status
        #    options: finished, running, waiting
        self.simulation_status = "waiting"

        #   parameters
        self._parameters = []

    def _containers(self):
        """
        Function predefines containers for solution data
        :return:
        """
        #   step number
        self._step_num_solution_container = []
        #   mechanical energy of the system
        self._energy_solution_container = []
        #   numerical integration error
        self._R_solution_container = []
        #   step size
        self._step_size_solution_container = []
        #   time
        self._t_solution_container = []
        #   positions and velocities
        self._q_sol_container = []

    def add_parameters(self, parameters):
        """
        Function adds additional parameters
        """
        self._parameters.append(parameters)

    def set_filetype(self, filetype):
        """
        Function sets filetype (extension) of solution data file
        :param filetype:
        :return:
        """
        self._filetype = filetype

    def read_file(self, _file_abs_path):
        """
        Read (load) solution data from file
        """
        #   load data from file if already exists
        _path, _file = os.path.split(_file_abs_path)
        self._name, self._filetype = os.path.splitext(_file)

        if self._filetype == ".dat" or self._filetype == ".sol":
            data = self._read_ascii_file(_file_abs_path)
        elif self._filetype == ".xlsx":
            data = self._read_excel_file(_file_abs_path)
        elif self._filetype == ".csv":
            data = self._read_csv_file(_file_abs_path)
        else:
            raise ValueError, "Filetype not supported."

        self.solution_data = data

        self._step_num_solution_container = data[:, 0]
        #   mechanical energy of the system
        self._energy_solution_container = data[:, 1]
        #   numerical integration error
        self._R_solution_container = data[:, 2]
        #   step size
        self._step_size_solution_container = data[:, 3]
        #   time
        self._t_solution_container = data[:, 4]
        #   positions and velocities
        self._q_sol_container = data[:, 5::]

    def _read_ascii_file(self, _file):
        """
        Load data from file (filetype .dat, .sol), saved with numpy savetxt
        """
        data = np.loadtxt(str(_file), skiprows=2)
        return data

    def _read_excel_file(self):
        print "Under Construction: ", os.path.realpath(__file__)

    def _read_csv_file(self):
        print "Under Construction: ", os.path.realpath(__file__)
    
    def add_data(self, data):
        """
        
        """
        self.solution_data = data
    
    def write_to_file(self):
        """
        Function saves data to file
        """
        
    def _write_to_ascii_file(self, data):
        """
        Function saves data to .dat filetype
        """
        self.solution_data = data

        print "self.MBD_system._solution_filetype =", self.MBD_system._solution_filetype
        print "self.simulation_id =", self.simulation_id
        self.solution_filename = 'solution_data_%02d'%self.simulation_id+self.MBD_system._solution_filetype

        self.solution_filename = check_filename(self.solution_filename)

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

        np.savetxt(self.solution_filename, data, fmt=__frmt, delimiter='\t', header = __header, comments = __comments)
        
    def _write_to_excel_file(self):
        """
        Function saves data to .xlsx filetype
        """
        
    def _write_to_csv_file(self):
        """
        Function saves data to .csv filetype
        """

    def _plot(self, x, y, color=None, label=None):
        """

        :param x:
        :param y:
        :param color:
        :param label:
        :return:
        """
        if x == "t":
            x_data = self._t_solution_container
        elif x == "Rx":
            x_data = self._q_sol_container[:, 0]
        else:
            raise ValueError, "x not corrent!"

        if y == "dRx":        # print self._q_sol_containerx":
            y_data1 = self._q_sol_container[:, 6]
            y_data2 = self._q_sol_container[:, 9]
            plt.plot(x_data, y_data1, marker=None, color=color, label=label)
            plt.plot(x_data, y_data2, marker=None, color=color)
        elif y == "Ry":
            y_data = self._q_sol_container[:, 1]
            plt.plot(x_data, y_data, marker=None, color=color, label=label)
        elif y == "energy":
            y_data = self._energy_solution_container

            plt.plot(x_data, y_data, marker=None, color=color, label=label)

    def _save_animation_file(self, folder=None):
        """
        Function saves
        :param folder:  absolute path to folder that contains figures to be merged in animation movie file
        :return:        None
        """




if __name__ == "__main__":
    # _file = "solution_data_example.sol"
    # _file = "solution_data_01_hertz.sol"
    sol = SolutionData()

    # sol.read_file(_file)
    #
    # fig = plt.figure(num=1, figsize=(8, 6), dpi=100, facecolor='w', edgecolor='k')
    # ax = plt.subplot(111)
    # ax.ticklabel_format(style='sci',scilimits=(-4,4), axis='both')

    # plt.plot(sol._t_solution_container, sol._energy_solution_container, label="E")
    # plt.plot(sol._t_solution_container, sol._q_sol_container[:, 0], label="Rx_0")
    # plt.plot(sol._t_solution_container, sol._q_sol_container[:, 6], label="dRx_0")
    # plt.plot(sol._t_solution_container, sol._q_sol_container[:, 9], label="dRx_1")
    # plt.plot(sol._t_solution_container, sol._q_sol_container[:, 0], label="Rx_0")

    # ax.legend(loc='best', fontsize=12)
    # plt.grid()
    # plt.show()

    folder = "c:\\Users\\lskrinjar\\Dropbox\\DyS\\dynamic_systems\\0_2_2\\screenshots"
    images = []
    for _file in os.listdir(folder):
        _file_path = os.path.join(folder, _file)
        images.append(_file_path)

    print "images ="
    print images
    video = ImageSequenceClip(images, fps=24)
    video.write_videofile("test.avi",codec='mpeg4')
    print "finished"
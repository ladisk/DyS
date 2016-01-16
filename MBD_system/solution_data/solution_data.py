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
import xlsxwriter
import xlrd

from MBD_system.MBD_system import *
from MBD_system.MBD_system_items import SolutionDataItem
from MBD_system import convert_bytes_to_
from MBD_system.check_filename import check_filename


class SolutionData(SolutionDataItem):
    """
    classdocs
    """
    __id = itertools.count(0)

    def __init__(self, name=None, _file=None, MBD_system=None, parent=None):
        """
        Constructor of solution data class
        """
        super(SolutionData, self).__init__(name, parent)
        self._parent = parent

        #   MBD system pointer
        self.MBD_system = MBD_system

        #   solution id
        self._id = self.__id.next()

        if name is None:
            self._name = "solution_data_%.2i"%self._id
        else:
            self._name = name

        #   solution type
        if self._parent is None:
            self._type = "solution"
        else:
            self._type = self._parent._type

        #   abs path to solution file
        self._file_abs_path = _file

        #    set file type: .dat, .xlsx, .csv, .sol(default)
        self._filetype = ".sol"

        #   load data from file if already exists
        self.solution_data = None
        if os.path.isfile(str(self._file_abs_path)):
            _path, _file = os.path.split(self._file_abs_path)
            self.filename, self._filetype = os.path.splitext(_file)
            self.solution_data = self.read_file(self._file_abs_path)

        #    simulation status
        #    options: finished, running, waiting
        self.waiting = True
        self.running = False
        self.finished = False
        self.simulation_status = "waiting"

        #   animation status
        self.loaded = False
        self.animation = False

        #   parameters
        self._parameters = []

        #   number of bodies
        self.n_b = None

        #   number of absolute coordinates of the system
        self.n = None

        #   fixed header size
        self.__header0_size = 5

        #   header ids
        self._header_ids = []

        self.__excel_worksheet = "solution"

    def load_solution_data(self):
        self.loaded = True
        return self.solution_data

    def getFileSizeInBytes(self):
        """

        :return:
        """
        _file = self._name+self._filetype
        if os.path.isfile(_file):
            return os.path.getsize(_file)
        else:
            return 0

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
        self._file_abs_path = _file_abs_path

        #   load data from file if already exists
        _path, _file = os.path.split(_file_abs_path)
        self._name, self._filetype = os.path.splitext(_file)

        if self._filetype == ".dat" or self._filetype == ".sol":
            data = self._read_ascii_file()
        elif self._filetype == ".xlsx":
            data = self._read_excel_file()
        elif self._filetype == ".csv":
            data = self._read_csv_file()
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

        #   status
        self.loaded = True

    def _read_ascii_file(self):
        """
        Load data from file (filetype .dat, .sol), saved with numpy savetxt
        """
        data = np.loadtxt(str(self._file_abs_path), skiprows=2)
        return data

    def _read_excel_file(self):
        """

        :return:
        """
        self.__excel_settings()

        #   open excel workbook
        workbook = xlrd.open_workbook(self._file_abs_path)

        #   get solution worksheet
        worksheet = workbook.sheet_by_name(self.__excel_worksheet)

        #   step number
        _step_num_solution_container = worksheet.col_values(colx=self.__col_step_num_solution_container, start_rowx=self.__column_start_write)
        #   mechanical energy
        _energy_solution_container = worksheet.col_values(colx=self.__col_mechanical_energy_solution_container, start_rowx=self.__column_start_write)
        #   error
        _R_solution_container = worksheet.col_values(colx=self.__col_error_solution_container, start_rowx=self.__column_start_write)
        #   step size
        _step_size_solution_container = worksheet.col_values(colx=self.__col_step_size_solution_container, start_rowx=self.__column_start_write)
        #   time
        _t_solution_container = worksheet.col_values(colx=self.__col_t_solution_container, start_rowx=self.__column_start_write)

        # print np.array([_step_num_solution_container, _energy_solution_container])
        #   read headers
        self.n = len(worksheet.row_values(self.__column_start_write - 1)) - self.__header0_size
        # print self.n

        _q_sol_container = []
        for i in range(self.__header0_size, len(worksheet.row_values(self.__column_start_write - 1))):
            col = np.array(worksheet.col_values(colx=i, start_rowx=self.__column_start_write))
            _q_sol_container.append(col)

        _data = np.array([_step_num_solution_container,
                        _energy_solution_container,
                        _R_solution_container,
                        _step_size_solution_container,
                        _t_solution_container])
        print np.shape(_data), type(_data), type(np.array(_q_sol_container)), np.shape(_q_sol_container)
        # print _data.T
        # print _q_sol_container
        data = np.hstack((_data.T, np.array(_q_sol_container).T))
        return data

    def _read_csv_file(self):
        print "Under Construction: ", os.path.realpath(__file__)
    
    def add_data(self, data):
        """

        :param data:

        """
        self.solution_data = data
    
    def write_to_file(self, _file_abs_path=None):
        """
        Function saves data to file

        :param _file_abs_path:
        """
        if _file_abs_path != None:
            self._file_abs_path = _file_abs_path

        #   number of coordinates and velocities of MBD system
        [steps, n] = np.shape(self.solution_data)
        self.n = n - self.__header0_size

        #   nuimber of bodies
        self.n_b = self.n / 6.

        #   create header ids
        self._header_ids = np.arange(0, self.n_b)

        #   check filename
        #   check filename
        if self._file_abs_path is None:
            self._filename = self._name+self._filetype
        else:
            self._filename = self._file_abs_path

        if self._filetype not in self._filename:
            self._filename = self._filename + self._filetype
        self._filename = check_filename(self._filename)

        #   save to file of selected type
        if self._filetype == ".dat" or self._filetype == ".sol":
            data = self._write_to_txt_file()
        elif self._filetype == ".xlsx":
            data = self._write_to_excel_file()
        elif self._filetype == ".csv":
            data = self._read_csv_file()
        else:
            raise ValueError, "Filetype not supported."

        self.simulation_status = "saved to file"

        #   write info to log file
        logging.getLogger("DyS_logger").info("Solution data saved to file: %s. Size is %s", self._filename, convert_bytes_to_.convert_size(os.path.getsize(self._filename)))
        
    def _write_to_txt_file(self):
        """
        Function saves data to .dat filetype
        """
        #   set headers
        self._headers = self.__txt_headers()

        #    format for each column
        __frmt = ['%5i']+['%20.16f']+['%20.16f']+['%20.16f']+['%20.16f']+['%.10E']*self.n

        #   add new line at the end of comment
        self._comments = self._comments + "\n"

        #   write data to txt file
        np.savetxt(self._filename, self.solution_data, fmt=__frmt, delimiter='\t', header=self._headers, comments=self._comments)

    def __txt_headers(self):
        """
        Function creates headers for txt file
        :return:
        """
        header = "i-th step mechanical  energy  \t  R \t\t\t\t\t  dt \t\t\t\t\t  time \t"

        #   header for q of a MBD system
        for _id in self._header_ids:
            _id_str = str(int(_id))
            q_i_header = "\t\t\t\tRx_" + _id_str + "\t\t\t\tRy_" + _id_str + "\t\t\t\ttheta_" + _id_str

            header += q_i_header

        #    add header for dq
        for body in self._header_ids:
            _id_str = str(int(_id))
            dq_i_header = "\t\t\t\tdRx_"+_id_str + "\t\t\t\tdRy_" + _id_str + "\t\t\t\tomega_"+_id_str

            header += dq_i_header

        return header

    def _write_to_excel_file(self):
        """
        Function saves data to .xlsx filetype
        """
        #   set up excel file
        self.__excel_settings()

        #   create an new Excel file and add a worksheet
        workbook = xlsxwriter.Workbook(self._filename, {'nan_inf_to_errors': True})
        format_1 = workbook.add_format({'num_format': '#0.000000000000000'})
        format_2 = workbook.add_format({'num_format': '#0.00000000'})
        #   add worksheet
        worksheet = workbook.add_worksheet(self.__excel_worksheet)

        #   write comments
        self._comments = "integration method: %s"%self.MBD_system.integrationMethod
        if self._comments == "":
            try:
                comments = self._parent._parent._name
            except:
                comments = ""
        else:
            comments = self._comments
        worksheet.write(0, 0, comments)

        #   write header
        worksheet.write_row(1, 0, self._headers)

        #   write solution data to columns
        #   step number
        worksheet.write_column(self.__column_start_write, self.__col_step_num_solution_container, self.solution_data[:,0])
        #   mechanical energy
        worksheet.write_column(self.__column_start_write, self.__col_mechanical_energy_solution_container, self.solution_data[:,1])
        worksheet.set_column('B:B', 20)
        #   error
        worksheet.write_column(self.__column_start_write, self.__col_error_solution_container, self.solution_data[:,2])
        #   step size
        worksheet.write_column(self.__column_start_write, self.__col_step_size_solution_container, self.solution_data[:,3])
        #   time
        worksheet.write_column(self.__column_start_write, self.__col_t_solution_container, self.solution_data[:,4])

        #   q, dq MBD system solution data
        for i in range(self.__col_t_solution_container, self.__col_t_solution_container+self.n+1):
            worksheet.write_column(self.__column_start_write, i, self.solution_data[:,i])

        #   freeze first two rows
        worksheet.freeze_panes(2, 0)
        #   close file
        workbook.close()

    def __excel_settings(self):
        """
        Function sets some basic excel document settings
        :return:
        """
        #   set up headers
        self._headers = self.__excel_headers()

        #   row number to start writing data
        self.__column_start_write = 2

        #   columns
        self.__col_step_num_solution_container = 0
        self.__col_mechanical_energy_solution_container = 1
        self.__col_error_solution_container = 2
        self.__col_step_size_solution_container = 3
        self.__col_t_solution_container = 4

    def __excel_headers(self):
        """
        Function creates headers for excel file
        :return:
        """
        header = ["i-th step", "mechanical  energy",  "R", "dt", "time"]

        #   header for q of a MBD system
        for _id in self._header_ids:
            _id_str = str(int(_id))
            q_i_header = ["Rx_"+_id_str, "Ry_"+_id_str, "theta_"+_id_str]

            header += q_i_header

        #   header for dq of a MBD system
        for body in self._header_ids:
            _id_str = str(int(_id))
            dq_i_header = ["dRx_"+_id_str, "dRy_"+_id_str, "dtheta_"+_id_str]

            header += dq_i_header

        return header

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

    def _create_animation_file(self, folder=None):
        """
        Function saves
        :param folder:  absolute path to folder that contains figures to be merged in animation movie file
        :return:        None
        """
        print "SolutionData._create_animation_file()"

if __name__ == "__main__":
    # _file = "solution_data_example.sol"
    # _file = "solution_data_01_hertz.sol"
    sol = SolutionData()

    # sol.read_file(_file)
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

    #   create movie from saved screenshoots
    folder = "c:\\Users\\lskrinjar\\Dropbox\\DyS\\dynamic_systems\\screenshots"
    images = []
    for _file in os.listdir(folder):
        _file_path = os.path.join(folder, _file)
        images.append(_file_path)
    #
    # print "images ="
    # print images
    # video = ImageSequenceClip(images, fps=24)
    video = ImageSequenceClip(images, fps=1)
    video.write_videofile(folder+"\\"+"test.avi",codec='mpeg4')
    print "finished"

    #   save solution data to xlsx file
    # sol.setName("save2excel")
    # sol.solution_data = np.random.rand(11, 11)
    # sol.set_filetype(".xlsx")
    # pprint(vars(sol))

    # sol.write_to_file()
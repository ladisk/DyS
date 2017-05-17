"""
Created on 9. jul. 2014

@author: lskrinjar
"""
import itertools
import os
import datetime
from ast import literal_eval


import numpy as np
import xlsxwriter
import xlrd
from PyQt4 import QtCore
from matplotlib import pyplot as plt
import logging
import threading
from collections import OrderedDict

# try:
#     from moviepy.editor import *
# except:
moviepy = None

from MBD_system.MBD_system_items import SolutionDataItem
from MBD_system.check_filename import check_filename
from MBD_system import convert_bytes_to_


class SolutionData(SolutionDataItem):
    """
    classdocs
    """
    __id = itertools.count(0)

    def __init__(self, name=None, _file=None, MBD_item=None, parent=None):
        """
        Constructor of solution data class
        """
        super(SolutionData, self).__init__(name, parent)
        self._parent = parent

        #   MBD system pointer
        self.MBD_item = MBD_item

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

        self.workbook = None

        #   abs path to solution file
        self._file_abs_path = _file

        #    set file type: .dat, .xlsx, .csv, .sol(default)
        self._qfiletypes = self.qfiletypes()
        self._filetypes = [".dat", ".xlsx", ".csv", ".sol"]
        self._filetype = ".xlsx"
        if hasattr(self._parent, "_solution_filetype"):
            self._filetype = self._parent._solution_filetype

        if hasattr(self.MBD_item, "_solution_filetype"):
            self._filetype = self.MBD_item._solution_filetype

        #   saving options
        #   options:
        #   discard
        #   overwrite
        #   save to new
        self._solution_save_options = "discard"

        #   load data from file if already exists
        self.solution_data = None
        if os.path.isfile(str(self._file_abs_path)):
            _path, _file = os.path.split(self._file_abs_path)
            self._filename, self._filetype = os.path.splitext(_file)

            # #   read file
            # self.read_file(self._file_abs_path)
        else:
            self._filename = self._name

        self.filename = self._filename + self._filetype

        #   temp solution file properties
        #   filename
        self._temp_filename = self._filename + ".sol"
        #   delimiter
        self._txt_delimitier = '\t'

        #   color
        self.color = None

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
        #
        self.cols = 7

        #   number of absolute coordinates of the system
        self.n = None

        #   skip rows
        self._skip_rows = 2
        self._skip_rows_temp = 0

        self._excel_worksheet = "data"
        
        #   header ids
        self._header_ids = []
        #   fixed header
        self._header0 = ["step",
                         "Em",
                         "Ek",
                         "Ep",
                         "Ees",
                         "Error",
                         "h level",
                         "h",
                         "t",
                         "iterNum",
                         "contact"]

        #    initialize containers attributes
        self._containers()
        self.solution_container_names = ["_step_num_solution_container",
                                         "_mechanical_energy_solution_container",
                                         "_kinetic_energy_solution_container",
                                         "_potential_energy_solution_container",
                                         "_elastic_strain_energy_solution_container",
                                         "_absError_solution_container",
                                         "_h_level_solution_container",
                                         "_h_solution_container",
                                         "_t_solution_container",
                                         "_iterNumber_solution_container",
                                         "_contact_status_solution_container"]

    def _clear_file(self):
        """
        Function checks if file already has data before writing to it during simulation and clears it before simulations starts
        :return:
        """
        if os.path.isfile(self._temp_filename):
            if os.path.getsize(self._temp_filename) > 0:
                with open(self._temp_filename, "w"):
                    pass

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

    def print_solution_containers(self):
        """

        :return:
        """
        for sol_cont_name in self.solution_container_names:
            print sol_cont_name
            print getattr(self, sol_cont_name), type(getattr(self, sol_cont_name))

    def _containers(self):
        """
        Function predefines containers for solution data
        :return:
        """
        self._step_num_solution_container = []
        self._mechanical_energy_solution_container = []
        self._kinetic_energy_solution_container = []
        self._potential_energy_solution_container = []
        self._elastic_strain_energy_solution_container = []
        self._absError_solution_container = []
        self._iterNumber_solution_container = []
        self._h_solution_container = []
        self._t_solution_container = []
        self._h_level_solution_container = []
        self._contact_status_solution_container = []
        self._iterNumber_solution_container = []

        #   positions and velocities
        self._q_solution_container = []

        #   size of vector of generalized (nodal and natural) coordinates
        self.n_q = None
        self.n_dq = None

    def _set_initial_data(self, Em0, Ek0, Ep0, Ees0, q0):
        """

        :return:
        """
        self._step_num_solution_container.append(0)
        self._mechanical_energy_solution_container.append(Em0)
        self._kinetic_energy_solution_container.append(Ek0)
        self._potential_energy_solution_container.append(Ep0)
        self._elastic_strain_energy_solution_container.append(Ees0)
        self._absError_solution_container.append(0)
        self._iterNumber_solution_container.append(0)
        self._h_solution_container.append(0)
        self._t_solution_container.append(0)
        self._h_level_solution_container.append(0)
        self._contact_status_solution_container.append(0)
        self._q_solution_container.append(q0)

    def _track_data(self, stepNum, Em, Ek, Ep, Ees, absError, h, t, h_level, contact_status, iterNumber, q):
        """

        :return:
        """
        self._step_num_solution_container.append(stepNum)
        self._mechanical_energy_solution_container.append(Em)
        self._kinetic_energy_solution_container.append(Ek)
        self._potential_energy_solution_container.append(Ep)
        self._elastic_strain_energy_solution_container.append(Ees)
        self._absError_solution_container.append(absError)
        self._iterNumber_solution_container.append(iterNumber)
        self._h_solution_container.append(h)
        self._t_solution_container.append(t)
        self._h_level_solution_container.append(h_level)
        self._contact_status_solution_container.append(contact_status)
        self._q_solution_container.append(q)

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

    def read_file(self, file_abs_path=None):
        """
        Read (load) solution data from file
        """
        if file_abs_path is not None:
            self._file_abs_path = file_abs_path

        #   load data from file if already exists
        _path, self._file = os.path.split(self._file_abs_path)
        self._name, self._filetype = os.path.splitext(self._file)

        if self._filetype == ".dat" or self._filetype == ".sol":
            data = self._read_ascii_file()

        elif self._filetype == ".xlsx":
            data = self._read_excel_file()

        elif self._filetype == ".csv":
            data = self._read_csv_file()

        else:
            raise ValueError, "Filetype not supported."

        if data is None:
            print "Variable data is None!"

        else:
            self.add_data(data)
            print "Solution read from file: %s"%self._file_abs_path

    def add_data(self, data):
        """
        Function transforms solution data ndarray (matrix) to solution containers
        :return:
        """
        self.solution_data = data

        if data.ndim == 2:
            n = len(self.solution_container_names)
            for i in range(0, n):
                name = self.solution_container_names[i]

                #   integer
                if name in ["_step_num_solution_container",
                            "_step_size_solution_container",
                            "_h_level_solution_container"
                            ]:
                    col = np.array(data[:, i], dtype="float64")
                    setattr(self, name, col)

                #   float
                elif name in ["_mechanical_energy_solution_container",
                              "_kinetic_energy_solution_container",
                              "_potential_energy_solution_container",
                              "_elastic_strain_energy_solution_container",
                              "_absError_solution_container",
                              "_iterNumber_solution_container",
                              "_h_solution_container",
                              "_t_solution_container",
                              "_contact_status_solution_container"
                              ]:

                    #   check if data is empty - only first element in list is checked
                    if data[0, i] == "":
                        col = []
                        print name, "data is empty!"
                    else:
                        col = data[:, i]

                    setattr(self, name, np.array(col, dtype="float64"))

                else:
                    raise ValueError, "Wrong solution container name:%s" %name

            # self._step_num_solution_container = data[:, 0]
            # #   mechanical energy of the system
            # self._energy_solution_container = np.array(data[:, 1], dtype="float64")
            # #   numerical integration error
            # self._R_solution_container = np.array(data[:, 2], dtype="float64")
            # #   step size
            # self._step_size_solution_container = np.array(data[:, 3])
            # #   time
            # self._t_solution_container = np.array(data[:, 4]).astype("float64")
            # #   h level
            # self._dt_level_solution_container = np.array(data[:, 5], dtype="float64")
            # #   contact status in MBD system
            # self._contact_status_solution_container = np.array(data[:, 6])
            #
            # #   iter number
            # self._iterNumber_solution_container = np.array(data[:, 7])

            #   positions and velocities
            self._q_solution_container = np.array(data[:, len(self._header0)+1::], dtype="float64")

            #   number of absolute nodal coordinates
            self.n = len(self._q_solution_container[0])
            self.n_q = self.n_dq = self.n / 2

            #   status
            self.loaded = True

    def _read_ascii_file(self, filename=None):
        """
        Load data from file (filetype .dat, .sol), saved with numpy savetxt
        """
        if filename is None:
            filename = self._file_abs_path
        print filename
        data = np.loadtxt(str(filename), skiprows=0)
        return data

    def _read_excel_file(self):
        """

        :return:
        """
        self._excel_settings()

        #   open excel workbook
        self.workbook = xlrd.open_workbook(self._file_abs_path)

        #   get solution worksheet
        worksheet = self.workbook.sheet_by_name(self._excel_worksheet)

        #   step number
        _step_num_solution_container = worksheet.col_values(colx=self._col_step_num_solution_container, start_rowx=self._column_start_write)

        #   mechanical energy
        _mechanical_energy_solution_container = worksheet.col_values(colx=self._col_mechanical_energy_solution_container, start_rowx=self._column_start_write)
        # print "_mechanical_energy_solution_container =", type(_mechanical_energy_solution_container[0])
        #   kinetic energy
        _kinetic_energy_solution_container = worksheet.col_values(colx=self._col_kinetic_energy_solution_container, start_rowx=self._column_start_write)
        # print "_kinetic_energy_solution_container =", _kinetic_energy_solution_container
        # print "_kinetic_energy_solution_container =", len(_kinetic_energy_solution_container)
        # print "_kinetic_energy_solution_container =", type(_kinetic_energy_solution_container[0])
        #   potential energy
        _potential_energy_solution_container = worksheet.col_values(colx=self._col_potential_energy_solution_container, start_rowx=self._column_start_write)

        #   elastic strain energy
        _elastic_strain_energy_solution_container = worksheet.col_values(colx=self._col_elastic_strain_energy_solution_container, start_rowx=self._column_start_write)

        #   error
        _absError_solution_container = worksheet.col_values(colx=self._col_absError_solution_container, start_rowx=self._column_start_write)

        #   h level
        _h_level_solution_container = worksheet.col_values(colx=self._col_h_level_solution_container, start_rowx=self._column_start_write)

        #   step size
        _h_solution_container = worksheet.col_values(colx=self._col_h_solution_container, start_rowx=self._column_start_write)

        #   time
        _t_solution_container = worksheet.col_values(colx=self._col_t_solution_container, start_rowx=self._column_start_write)

        #   iter number
        _iterNumber_solution_container = worksheet.col_values(colx=self._col_iterNumber_solution_container, start_rowx=self._column_start_write)

        #   contact status
        _contact_status_solution_container = worksheet.col_values(colx=self._col_contact_status_solution_container, start_rowx=self._column_start_write)

        #   read headers
        self.n = len(worksheet.row_values(self._column_start_write - 1)) - self._header0_size
        # print self.n

        _q_sol_container = []
        for i in range(self._header0_size - 1, len(worksheet.row_values(self._column_start_write - 1))):
            col = np.array(worksheet.col_values(colx=i, start_rowx=self._column_start_write))
            _q_sol_container.append(col)

        _data = np.array([_step_num_solution_container,
                          _mechanical_energy_solution_container,
                          _kinetic_energy_solution_container,
                          _potential_energy_solution_container,
                          _elastic_strain_energy_solution_container,
                          _absError_solution_container,
                          _h_level_solution_container,
                          _h_solution_container,
                          _t_solution_container,
                          _iterNumber_solution_container,
                          _contact_status_solution_container])

        data = np.hstack((_data.T, np.array(_q_sol_container).T))
        return data

    def _read_csv_file(self):
        print "Under Construction: ", os.path.realpath(__file__)
        return None

    def _write_to_file(self, _file_abs_path=None):
        """

        :return:
        """
        if _file_abs_path is not None:
            self._file_abs_path = _file_abs_path

        #   number of coordinates and velocities of MBD system
        if self.MBD_item.typeInfo() == "MBDsystem":
            if self.solution_data.ndim == 2:
                [steps, self.cols] = np.shape(self.solution_data)
            else:
                steps = 0
                self.cols = len(self.solution_data)
                self.solution_data = np.array([self.solution_data])

        #   check filename
        if self._file_abs_path is None:
            self._filename = self._name + self._filetype
        else:
            self._filename = self._file_abs_path

        if self._filetype not in self._filename:
            self._filename = self._filename + self._filetype

        if type(self._filename) == QtCore.QString:
            self._filename = str(self._filename)
        else:
            if self._solution_save_options == "save to new":
                self._filename = check_filename(self._filename)

        #   save to file of selected type
        if self._filetype == ".dat" or self._filetype == ".sol":
            self._write_to_txt_file()

        elif self._filetype == ".xlsx":
            self._write_to_excel_file()

        elif self._filetype == ".csv":
            self._read_csv_file()

        else:
            raise ValueError, "Filetype not supported."

        self.simulation_status = "saved to file"

        #   write info to log file
        logging.getLogger("DyS_logger").info("Solution data saved to file: %s. Size is %s", self._filename, convert_bytes_to_.convert_size(os.path.getsize(self._filename)))

    def write_to_file(self, _file_abs_path=None):
        """
        Function saves data to file

        :param _file_abs_path:
        """
        write_to_file_thread = SaveSolutionThread(sol_obj=self)
        write_to_file_thread.start()

    @staticmethod
    def qfiletypes():
        """
        Function returns filetypes of solution file that are supported in form of to be used in QFileDialog
        :return:
        """
        filetypes = "Text File (*.sol);;Text File (*.dat);;MS Excel (*.xlsx);;CSV (*.csv)"
        return filetypes
        
    def _write_to_txt_file(self):
        """
        Function saves data to .dat filetype
        """
        #   set headers
        self._headers = self._txt_headers()

        #    format for each column
        __frmt = ['%5i']+['%20.16f']+['%20.16f']+['%20.16f']+['%20.16f']+['%.10E']*self.n

        #   add new line at the end of comment
        self._comments = self._comments + "\n"

        #   write data to txt file
        np.savetxt(self._filename, self.solution_data, fmt=__frmt, delimiter='\t', header=self._headers, comments=self._comments)

    def _txt_headers(self):
        """
        Function creates headers for txt file
        :return:
        """
        self._txt_delimitier = '\t'
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
        self._excel_settings()

        #   create an new Excel file and add a worksheet
        self.workbook = xlsxwriter.Workbook(self._filename, {'nan_inf_to_errors': True})
        format_1 = self.workbook.add_format({'num_format': '#0.000000000000000'})
        format_2 = self.workbook.add_format({'num_format': '#0.00000000'})
        #   add worksheet
        solution_worksheet = self.workbook.add_worksheet(self._excel_worksheet)

        #   write comments
#         self._comments = "integration method: %s"%self.MBD_system.integrationMethod
        if self._comments == "":
            try:
                comments = self._parent._parent._name
            except:
                comments = ""
        else:
            comments = self._comments

        solution_worksheet.write(0, 0, comments)

        #   write index of vector q of MBD system
        indx = np.arange(0, self.n, dtype=int)
        solution_worksheet.write_row(0, self._header0_size, indx)

        #   write header
        solution_worksheet.write_row(1, 0, self._headers)

        #   write solution data to columns information columns
        for i, cont in enumerate(self.solution_container_names):
            solution_worksheet.write_column(self._column_start_write, i, self.solution_data[:,i])

        #   q, dq MBD system solution data
        for i in range(len(self.solution_container_names), self.cols):
            solution_worksheet.write_column(self._column_start_write, i, self.solution_data[:,i])

        #   freeze first two rows
            solution_worksheet.freeze_panes(2, 0)

        #   write MBD system data
        self._write_to_excel_file_MBD_system_properties()

        #   write group data
        self._write_to_excel_file_groups()

        #   write data of every body
        self._write_to_excel_file_bodies()

        #   set width of first column
        for worksheet in self.workbook.worksheets():
            if worksheet is not solution_worksheet:
                worksheet.set_column(0, 0, 30.)

        #   close file
        self.workbook.close()

    def _write_to_excel_file_groups(self):
        """

        :return:
        """
        for abs_path in self.MBD_item.abs_paths:
            if abs_path is not None:
                if os.path.isfile(abs_path) and abs_path != self.MBD_item.abs_path_to_bodies:
                    group_worksheet = self.workbook.add_worksheet(os.path.splitext(os.path.basename(abs_path))[0])

                    with open(abs_path, 'r') as file_:
                        for i, line in enumerate(file_):
                            if "FILE-START" in line:
                                pass
                            elif "FILE-END" in line:
                                break
                            else:
                                group_worksheet.write_string(i - 1, 0, line)

                        file_.close()

    def _write_to_excel_file_bodies(self):
        """

        :return:
        """
        for body in self.MBD_item.bodies:
            body_worksheet = self.workbook.add_worksheet(body._name)
            for i, (key, val) in enumerate(body._dict.iteritems()):
                body_worksheet.write_string(i, 0, key)
                body_worksheet.write_string(i, 1, str(val))

    def _write_to_excel_file_MBD_system_properties(self):
        """

        :return:
        """
        MBD_worksheet = self.workbook.add_worksheet("MBD_system")
        for i, (key, val) in enumerate(self.MBD_item.MBD_file_properties_dict.iteritems()):
            MBD_worksheet.write_string(i, 0, key)
            if isinstance(val, np.float):
                MBD_worksheet.write_number(i, 1, val)

            elif isinstance(val, bool):
                MBD_worksheet.write_boolean(i, 1, val)

            else:
                if isinstance(val, np.ndarray):
                    val = np.array2string(val)
                MBD_worksheet.write_string(i, 1, val)

        #   add computation time of simulation
        MBD_worksheet.write_string(i+1, 0, "Computation time")
        dt = self.MBD_item.end_time_simulation_info_in_sec_UTC - self.MBD_item.start_time_simulation_info_in_sec_UTC
        dt = datetime.timedelta(seconds=dt)
        MBD_worksheet.write_string(i+1, 1, str(dt))

    def _excel_settings(self):
        """
        Function sets some basic excel document settings
        :return:
        """
        #   set up headers
        self._headers = self._excel_headers()

        #   row number to start writing data
        self._column_start_write = 2

        #   number of column
        self._col_step_num_solution_container = 0
        self._col_mechanical_energy_solution_container = 1
        self._col_kinetic_energy_solution_container = 2
        self._col_potential_energy_solution_container = 3
        self._col_elastic_strain_energy_solution_container = 4
        self._col_absError_solution_container = 5
        self._col_h_level_solution_container = 6
        self._col_h_solution_container = 7
        self._col_t_solution_container = 8
        self._col_iterNumber_solution_container = 9
        self._col_contact_status_solution_container = 10

    def _excel_headers(self):
        """
        Function creates headers for excel file
        :return:
        """
        header = []
        self._header0_size = len(self._header0)

        for _header in self._header0:
            header.append(_header)

        #   number of bodies
        self.n_b = int((self.cols - self._header0_size) / 6)

        #   header for q of a MBD system
        if self.MBD_item is not None:
            for body in self.MBD_item.bodies:
                header += body._excel_header_q()

            #   header for dq of a MBD system
            for body in self.MBD_item.bodies:
                # header += dq_i_header
                header += body._excel_header_dq()

        return header

    def _write_to_csv_file(self):
        """
        Function saves data to .csv filetype
        """

    def _append_data_to_temp_file(self, data=None):
        """
        Function appends data to created file
        :param data:
        :return:
        """
        if data is None:
            data = self.collect_solution_data()

        #   open file to append data
        _file = open(self._temp_filename, "a")

        #   save to file
        np.savetxt(_file, data, delimiter=self._txt_delimitier, newline='\n')

        #   close file
        _file.close()

        #   remove all but last value in array/matrix
        self._reset_solution_data_after_append_to_temp_file()

    def collect_solution_data(self):
        """

        :return:
        """
        #   join in one ndarray (matrix)
        solution_data = np.hstack((np.array([self._step_num_solution_container]).T,
                                   np.array([self._mechanical_energy_solution_container]).T,
                                   np.array([self._kinetic_energy_solution_container]).T,
                                   np.array([self._potential_energy_solution_container]).T,
                                   np.array([self._elastic_strain_energy_solution_container]).T,
                                   np.array([self._absError_solution_container]).T,
                                   np.array([self._iterNumber_solution_container]).T,
                                   np.array([self._h_solution_container]).T,
                                   np.array([self._t_solution_container]).T,
                                   np.array([self._h_level_solution_container]).T,
                                   np.array([self._contact_status_solution_container]).T,
                                   self._q_solution_container))

        return solution_data

    def _reset_solution_data_after_append_to_temp_file(self):
        """
        Remove all but last value in solution containers
        :return:
        """
        for i in range(0, len(self.solution_container_names)):
            val = [getattr(self, self.solution_container_names[i])[-1]]
            setattr(self, self.solution_container_names[i], val)

        self._q_solution_container = [self._q_solution_container[-1]]

    def _solution_data_from_temp_file(self):
        """

        :return:
        """
        #   read data from temp file to store it in solution file (with filetype) specified by user
        data = np.loadtxt(str(self._temp_filename), skiprows=0)

        self.add_data(data)

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
            x_data = self._q_solution_container[:, 0]
        else:
            raise ValueError, "x not corrent!"

        if y == "dRx":        # print self._q_sol_containerx":
            y_data1 = self._q_solution_container[:, 6]
            y_data2 = self._q_solution_container[:, 9]
            plt.plot(x_data, y_data1, marker=None, color=color, label=label)
            plt.plot(x_data, y_data2, marker=None, color=color)
        elif y == "Ry":
            y_data = self._q_solution_container[:, 1]
            plt.plot(x_data, y_data, marker=None, color=color, label=label)
        elif y == "energy":
            y_data = self._mechanical_energy_solution_container

            plt.plot(x_data, y_data, marker=None, color=color, label=label)


class SaveSolutionThread(threading.Thread):

    def __init__(self, sol_obj):
        threading.Thread.__init__ (self)

        self.sol_obj = sol_obj

    def run(self):
        """

        :param sol_obj:
        :return:
        """
        self.sol_obj._write_to_file()


if __name__ == "__main__":
    # filename = "solution_data.xlsx"
    filename = "solution_data_00_rigid.xlsx"
    # _file = "solution_data_01_hertz.sol"
    sol = SolutionData()

    sol.read_file(filename)

    # for val in sol._h_solution_container:
    #     print val

    # for i in range(0, len(sol._t_solution_container)):
    #     print sol._t_solution_container[i], sol._q_solution_container[i, 1]
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
    # folder = "c:\\Users\\lskrinjar\\Dropbox\\DyS\\dynamic_systems\\screenshots"
    # images = []
    # for _file in os.listdir(folder):
    #     _file_path = os.path.join(folder, _file)
    #     images.append(_file_path)
    #
    # video.write_videofile(folder+"\\"+"test.avi",codec='mpeg4')
    # print "finished"

    #   save solution data to xlsx file
    # sol.setName("save2excel")
    # sol.solution_data = np.random.rand(11, 11)
    # sol.set_filetype(".xlsx")
    # pprint(vars(sol))

    # sol.write_to_file()
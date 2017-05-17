"""
Created on 29. mar. 2016

@author: luka.skrinjar
"""
import xlrd
import xlsxwriter
import numpy as np

from solution_data import SolutionData


class SolutionDataJoint(SolutionData):
    """
    classdocs
    """

    def __init__(self, name=None, MBD_item=None, parent=None):
        """
        Constructor
        """
        super(SolutionDataJoint, self).__init__(name=name, MBD_item=MBD_item, parent=parent)
        
    def _containers(self):
        """
        Function predefines containers for solution data
        :return:
        """
        #   fixed header size
        #    0 - step number
        #    1 - simulation time
        self._header0_size = 1
        
        self._step_num_solution_container = []
        self._t_solution_container = []
        
        self._Lx_solution_container = []
        self._Ly_solution_container = []

        self._solution_container_names_list = ["_step_num_solution_container",
                                               "_t_solution_container",
                                               "_Lx_solution_container",
                                               "_Ly_solution_container",
                                               "_Qc_i_x_solution_container",
                                               "_Qc_i_y_solution_container",
                                               "_Qc_i_zz_solution_container",
                                               "_Qc_j_x_solution_container",
                                               "_Qc_j_y_solution_container",
                                               "_Qc_j_zz_solution_container"]
    
    def _excel_headers(self):
        """
        Function creates headers for excel file
        :return:
        """
        header = ["i-th step", "time", "Lx", "Ly", "Q_c_i_x", "Q_c_i_y", "Q_c_i_zz", "Q_c_j_x", "Q_c_j_y", "Q_c_j_zz"]
        
        return header

    def add_data(self, data):
        """

        :param data:

        """
        self.solution_data = data
        
        self._step_num_solution_container = data[:, 0]
        self._t_solution_container = data[:, 1]
        
        self._Lx_solution_container = data[:, 2]
        self._Ly_solution_container = data[:, 3]

        self._Qc_i_x_solution_container = data[:, 4]
        self._Qc_i_y_solution_container = data[:, 5]
        self._Qc_i_zz_solution_container = data[:, 6]
        self._Qc_j_x_solution_container = data[:, 7]
        self._Qc_j_y_solution_container = data[:, 8]
        self._Qc_j_zz_solution_container = data[:, 9]
        
    def _write_to_excel_file(self):
        """
        Function saves data to .xlsx filetype
        """
        #   set up excel file
        self._excel_settings()

        #   create an new Excel file and add a worksheet
        workbook = xlsxwriter.Workbook(self._filename, {'nan_inf_to_errors': True})
        format_1 = workbook.add_format({'num_format': '#0.000000000000000'})
        format_2 = workbook.add_format({'num_format': '#0.00000000'})
        #   add worksheet
        worksheet = workbook.add_worksheet(self._excel_worksheet)

        #   write comments
#         self._comments = "integration method: %s"%self.MBD_system.integrationMethod
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
        # print "self.solution_data ="
        # print self.solution_data
        # print np.shape(self.solution_data)
        # print "self.solution_data[:,0] ="
        # print self.solution_data[:,0]
        worksheet.write_column(self._column_start_write, self._col_step_num_solution_container, self.solution_data[:,0].astype(int))
        worksheet.set_column('B:B', 20)
        #   time
        worksheet.write_column(self._column_start_write, 1, self._t_solution_container)

        #   lagrange multipliers
        worksheet.write_column(self._column_start_write, 2, self._Lx_solution_container)
        worksheet.write_column(self._column_start_write, 3, self._Ly_solution_container)

        #   vector Q_c on each body in joint
        for i, Q_c in enumerate([self._Qc_i_x_solution_container,
                                self._Qc_i_y_solution_container,
                                self._Qc_i_zz_solution_container,
                                self._Qc_j_x_solution_container,
                                self._Qc_j_y_solution_container,
                                self._Qc_j_zz_solution_container]):
            worksheet.write_column(self._column_start_write, 4+i, Q_c)

        #   freeze first two rows
        worksheet.freeze_panes(2, 0)

        #   close file
        workbook.close()

    def _excel_settings(self):
        """

        :return:
        """
        self.__workbook_name = "data"
        #   set up headers
        self._headers = self._excel_headers()

        #   row number to start writing data
        self._column_start_write = self.__column_start = 2

        #   column's indexes
        self._col_step_num_solution_container = 0
        self._col_t_solution_container = 1
        self._col_Lx_solution_container = 2
        self._col_Ly_solution_container = 3
        self._col_Qc_i_x_solution_container = 4
        self._col_Qc_i_y_solution_container = 5
        self._col_Qc_i_zz_solution_container = 6
        self._col_Qc_j_x_solution_container = 7
        self._col_Qc_j_y_solution_container = 8
        self._col_Qc_j_zz_solution_container = 9

        self._cols_containers_list = [self._col_step_num_solution_container,
                                       self._col_t_solution_container,
                                       self._col_Lx_solution_container,
                                       self._col_Ly_solution_container,
                                       self._col_Qc_i_x_solution_container,
                                       self._col_Qc_i_y_solution_container,
                                       self._col_Qc_i_zz_solution_container,
                                       self._col_Qc_j_x_solution_container,
                                       self._col_Qc_j_y_solution_container,
                                       self._col_Qc_j_zz_solution_container]

    def _read_excel_file(self):
        """
        Function reads excel file of joint solution data
        :return:
        """
        #   set default properties for excel format of data
        self._excel_settings()

        #   create workbook object
        workbook = xlrd.open_workbook(self._file_abs_path)
        #   get data sheet from workbook object
        worksheet = workbook.sheet_by_name(self.__workbook_name)

        #   size of first column
        self.__N_col0 = len(worksheet.col(0))-2

        for sol_container_name, col_sol_container in zip(self._solution_container_names_list, self._cols_containers_list):
            #   read column and returns values as a list
            _list = worksheet.col_values(col_sol_container,
                                         start_rowx=self.__column_start,
                                         end_rowx=self.__N_col0)

            # print _list
            #   convert list to numpy array
            __solution_container = np.array(_list)

            #   set values to attribute of object
            setattr(self, sol_container_name, __solution_container)

        return None
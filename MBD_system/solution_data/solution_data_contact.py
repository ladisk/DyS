"""
Created on 29. mar. 2016

@author: luka.skrinjar
"""
import os
from collections import OrderedDict
import xlrd
import xlsxwriter
import numpy as np
from matplotlib import pyplot as plt


from MBD_system.fix_string import fix_string
from MBD_system.solution_data.solution_data import SolutionData


class SolutionDataContact(SolutionData):
    """
    classdocs
    """

    def __init__(self, name=None, MBD_item=None, parent=None):
        """
        Constructor
        """
        super(SolutionDataContact, self).__init__(name=name, MBD_item=MBD_item, parent=parent)

        #   filename
        if self.MBD_item is None:
            self._name = name
        else:
            self._name = self._filename = self.MBD_item._name.lower()

        #   filename with type
        self.filename = self._filename + self._filetype
        self.contact_model = None

        #   read columns
        self.read_uP = True
        self.read_rP = True

        self._containers()

    def _get_label(self):
        """
        
        :return: 
        """
        if "-" not in self.contact_model:
            label = self.contact_model.capitalize()

        else:
            label = self.contact_model.title().replace("C", "c")

        return label

    def reset(self):
        """

        :return:
        """
        self._containers()

    def _containers(self):
        """
        Function predefines containers for solution data
        :return:
        """
        #   fixed header size
        #    0 - step number
        #    1 - simulation time
        self._header0_size = 1
        
        #    solution containers over integration times
        #    dictionary is shaped as
        #    key - solution container name
        #    val - name (string) of MBD_item attribute
        self._solution_container_dict = OrderedDict({"_step_num_solution_container": "_step",
                                        "_status_container": "status",
                                        "_step_size_solution_container": "_h",
                                        "_t_solution_container": "_t",
                                        "_distance_solution_container": "_delta",
                                        "_dqn_solution_container": "_dq_n",
                                        "_dqt_solution_container": "_dq_t",
                                        "_Fn_solution_container": "Fn",
                                        "_Ft_solution_container": "Ft",
                                        "_F_solution_container": "F",
                                        "_u_P_solution_container": "u_P_LCS_list",
                                        "_r_P_solution_container": "r_P_GCS_list"})

        self.solution_container_names = ["_step_num_solution_container",
                                        "_status_container",
                                        "_step_size_solution_container",
                                        "_t_solution_container",
                                        "_distance_solution_container",
                                        "_dqn_solution_container",
                                        "_dqt_solution_container",
                                        "_Fn_solution_container",
                                        "_Ft_solution_container",
                                        "_F_solution_container",
                                        "_u_P_solution_container",
                                        "_r_P_solution_container"]

        #   create object attributes for names in list
        for sol_cont in self.solution_container_names:
            setattr(self, sol_cont, [])

    def _track_data_main(self, t, step):
        """
        Function tracks main data during integration
        :return:    None
        """

    def _track_data(self):
        """
        Function tracks object specific data at each time step during integration
        :return:    None
        """
        for i, attr in enumerate(self.solution_container_names):
            #   get value of parent attribute
            val = getattr(self.MBD_item, self._solution_container_dict[attr])

            #   append value to list
            lst = getattr(self, attr)

            if attr == "_distance_solution_container":
                if type(val) is list:
                    lst.extend(val)
                else:
                    lst.append(val)

            else:
                lst.append(val)

        # print "------------------------"
        # # print getattr(self.MBD_item, self._solution_container_dict[0])
        # print "attr =", attr

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

        #   add values to object attributes in for loop as arrays are 1D
        for container, solution_data in zip(self._solution_container_names_list[0:-2], self.solution_data):
            setattr(self, container, solution_data)

        #   separately as arrays are 2D (matrix) shaped
        self._u_P_solution_container = self.solution_data[9:13, :]

        self._r_P_solution_container = self.solution_data[13::, :]

    def _write_to_excel_file(self):
        """
        Function saves data to .xlsx filetype
        :return:
        """
        self._excel_settings()

        #   column headers
        _header = ["i-th step", "status", "dt", "time", "delta", "dq_n", "dq_t", "Fn", "Ft", "Fx", "Fy", "uPi_x", "uPi_y", "uPj_x", "uPj_y", "rPi_x", "rPi_y", "rPj_x", "rPj_y"]

        #   create an new Excel file and add a worksheet
        workbook = xlsxwriter.Workbook(self.filename, {'nan_inf_to_errors': True})
        format_1 = workbook.add_format({'num_format': '#0.000000000000000'})
        format_2 = workbook.add_format({'num_format': '#0.00000000'})
        worksheet = workbook.add_worksheet(self.__workbook_name)

        #   write comments
        worksheet.write(0, 0, self.MBD_item.contact_model._type)
        #   write header
        worksheet.write_row(1, 0, _header)

        #   write solution data to columns
        #   step number
        worksheet.write_column(self.__column_start_write, self.__col_step_num_solution_container, self._step_num_solution_container)

        #   contact status
        worksheet.write_column(self.__column_start_write, self.__col_status_container, self._status_container)

        #   step size
        worksheet.write_column(self.__column_start_write, 2, self._step_size_solution_container, format_1)
        worksheet.set_column('C:C', 20)

        #   time
        worksheet.write_column(self.__column_start_write, self.__col_t_solution_container, self._t_solution_container, format_1)
        worksheet.set_column('D:D', 20)

        #   distance if no contact or delta when contact is present
        worksheet.write_column(self.__column_start_write, self.__col_distance_solution_container, np.array(self._distance_solution_container).flatten(), format_1)
        worksheet.set_column('E:E', 20)

        #   dqn
        worksheet.write_column(self.__column_start_write, self.__col_dqn_solution_container, self._dqn_solution_container, format_1)
        worksheet.set_column('F:F', 22)

        #   dqt
        worksheet.write_column(self.__column_start_write, self.__col_dqt_solution_container, self._dqt_solution_container, format_1)
        worksheet.set_column('G:G', 22)

        #   Fn
        worksheet.write_column(self.__column_start_write, self.__col_Fn_solution_container, self._Fn_solution_container, format_2)
        worksheet.set_column('H:H', 16)

        #   Ft
        worksheet.write_column(self.__column_start_write, self.__col_Ft_solution_container, self._Ft_solution_container, format_2)
        worksheet.set_column('I:I', 16)

        #   Fx
        worksheet.write_column(self.__column_start_write, self.__col_F_solution_container, np.array(self._F_solution_container)[:, 0], format_2)
        worksheet.set_column('J:J', 16)
        #   Fy
        worksheet.write_column(self.__column_start_write, self.__col_F_solution_container + 1, np.array(self._F_solution_container)[:, 1], format_2)
        worksheet.set_column('K:K', 16)

        #   uPi and uPj of contact force in LCS
        self._u_P_solution_container = np.array(self._u_P_solution_container).reshape(len(self._u_P_solution_container), 4)
        for i in range(0, 4):
            worksheet.write_column(self.__column_start_write, self.__col_uP_solution_container+i, self._u_P_solution_container[:, i], format_1)
        worksheet.set_column('L:O', 20)

        #   rPi and rPj of contact force in GCS
        self._r_P_solution_container = np.array(self._r_P_solution_container).reshape(len(self._r_P_solution_container), 4)
        for i in range(0, 4):
            worksheet.write_column(self.__column_start_write, self.__col_rP_solution_container+i, self._r_P_solution_container[:, i], format_1)
        worksheet.set_column("P:S", 22)

        #   freeze first two rows
        worksheet.freeze_panes(2, 0)

        #   close file
        workbook.close()

        print "Contact solution data of %s saved to file %s"%(self.MBD_item._name, self.filename)

    def _excel_settings(self):
        """

        :return:
        """
        self.__workbook_name = "data"

        #   row number to start writing data
        self.__column_start_write = self.__column_start = 2

        self.__col_step_num_solution_container = 0
        self.__col_status_container = 1
        self.__col_step_size_solution_container = 2
        self.__col_t_solution_container = 3
        self.__col_distance_solution_container = 4
        self.__col_dqn_solution_container = 5
        self.__col_dqt_solution_container = 6
        self.__col_Fn_solution_container = 7
        self.__col_Ft_solution_container = 8
        self.__col_F_solution_container = 9
        self.__col_uP_solution_container = 11
        self.__col_rP_solution_container = 15

        self.__cols_containers_list = [self.__col_step_num_solution_container,
                                        self.__col_status_container,
                                        self.__col_step_size_solution_container,
                                        self.__col_t_solution_container,
                                        self.__col_distance_solution_container,
                                        self.__col_dqn_solution_container,
                                        self.__col_dqt_solution_container,
                                        self.__col_Fn_solution_container,
                                        self.__col_Ft_solution_container,
                                        self.__col_F_solution_container,
                                        self.__col_uP_solution_container,
                                        self.__col_rP_solution_container]

    def read_file(self, file_abs_path=None):
        """


        :param filename:
        :return:
        """
        #   load data from file if already exists
        _path, _file = os.path.split(file_abs_path)
        self.filename, self._filetype = os.path.splitext(_file)

        if self._filetype == ".xlsx" or self._solution_filetype == ".xls":
            self._read_excel_file(file_abs_path)

    def _read_excel_file(self, file_path):
        """

        :return:
        """
        #   set default properties for excel format of data
        self._excel_settings()
        #   create workbook object
        workbook = xlrd.open_workbook(file_path)
        #   get data sheet from workbook object
        worksheet = workbook.sheet_by_name(self.__workbook_name)

        #   get contact model
        cell_value = worksheet.cell(0, 0).value
        if ":" in cell_value:
            __contact_model_type_info = cell_value.strip()[cell_value.index(":")+2::]
            self.contact_model = __contact_model_type_info.encode('ascii','ignore')

        else:
            self.contact_model = cell_value

            #   size of first column
        self.__N_col0 = len(worksheet.col(0))-2

        # for i, (sol_container_name, col_sol_container) in enumerate(zip(self.solution_container_names, self.__cols_containers_list)):
        #     #   read column and returns values as a list
        #     if sol_container_name == "_F_solution_container":
        #         _list_Fx = worksheet.col_values(col_sol_container+1,
        #                                      start_rowx=self.__column_start,
        #                                      end_rowx=None)
        #
        #         _list_Fy = worksheet.col_values(col_sol_container+2,
        #                                      start_rowx=self.__column_start,
        #                                      end_rowx=None)
        #         _list = [_list_Fx, _list_Fy]
        #     else:
        #         try:
        #             _list = worksheet.col_values(col_sol_container,
        #                                          start_rowx=2)
        #         except:
        #             pass
        #     if sol_container_name == "_status_container":
        #         _dtype = np.int32
        #     else:
        #         _dtype = np.float64
        #
        #     #   convert list to numpy array
        #     __solution_container = np.array(_list, dtype=_dtype)
        #
        #     #   set values to attribute of object
        #     setattr(self, sol_container_name, __solution_container)

        #   step number
        array = np.array(worksheet.col_values(self.__col_step_num_solution_container, start_rowx=2), dtype=np.int32)
        setattr(self, self.solution_container_names[0], array)

        #   contact status
        array = np.array(worksheet.col_values(self.__col_status_container, start_rowx=2))
        setattr(self, self.solution_container_names[1], array)

        #   step size
        array = np.array(worksheet.col_values(self.__col_step_size_solution_container, start_rowx=2), dtype=np.float64)
        setattr(self, self.solution_container_names[2], array)

        #   time
        array = np.array(worksheet.col_values(self.__col_t_solution_container, start_rowx=2), dtype=np.float64)
        setattr(self, self.solution_container_names[3], array)

        #   distance
        array = np.array(worksheet.col_values(self.__col_distance_solution_container, start_rowx=2), dtype=np.float64)
        setattr(self, self.solution_container_names[4], array)

        #   dqn
        array = np.array(worksheet.col_values(self.__col_dqn_solution_container, start_rowx=2), dtype=np.float64)
        setattr(self, self.solution_container_names[5], array)

        #   dqt
        array = np.array(worksheet.col_values(self.__col_dqt_solution_container, start_rowx=2), dtype=np.float64)
        setattr(self, self.solution_container_names[6], array)

        #   Fn
        array = np.array(worksheet.col_values(self.__col_Fn_solution_container, start_rowx=2), dtype=np.float64)
        setattr(self, self.solution_container_names[7], array)

        #   Ft
        array = np.array(worksheet.col_values(self.__col_Ft_solution_container, start_rowx=2), dtype=np.float64)
        setattr(self, self.solution_container_names[8], array)

        #   Fx
        Fx = worksheet.col_values(self.__col_F_solution_container, start_rowx=2)
        array_x = np.array(Fx)#.replace('', np.nan)
        Fy = np.array(worksheet.col_values(self.__col_F_solution_container + 1, start_rowx=2))
        array_y = np.array(Fy)
        setattr(self, self.solution_container_names[9], np.array([array_x, array_y]))

        #   uPi and uPj of contact force in LCS
        if self.read_uP:
            uP = []
            for i in range(0, 4):
                col = self.__col_uP_solution_container + i
                array = np.array(worksheet.col_values(col, start_rowx=2))
                uP.append(array)
            setattr(self, self.solution_container_names[10], np.array(uP))

        #   rPi and rPj of contact force in GCS
        if self.read_rP:
            rP = []
            for i in range(0, 4):
                array = np.array(worksheet.col_values(self.__col_rP_solution_container + i, start_rowx=2))
                rP.append(array)
            setattr(self, self.solution_container_names[11], np.array(rP))

    def plot(self, x, y, color=None, label=None, str_id=None):
        """

        :param x: options = (delta)
        :param y: options = (Fn, dqn)
        :param color:
        :return:
        """
        _contact_indices = self._contact()

        label, color = self._plot_info(color)
        label = fix_string(label)
        if x == "delta":
            x_data = abs(self._distance_solution_container[_contact_indices])[1::]
        elif x == "t":
            x_data = self._t_solution_container[_contact_indices] - self._t_solution_container[_contact_indices][0]
        else:
            raise ValueError, "x not corrent!"

        if y == "Fn":
            y_data = abs(self._Fn_solution_container[_contact_indices]) / 0.015
        elif y == "dqn":
            y_data = self._dqn_solution_container[_contact_indices][1::]
        else:
            raise ValueError, "y not corrent!"

        # plt.plot(x_data, y_data, marker=None, color=color, label=r"$\textrm{%s}$"%label)
        plt.plot(x_data, y_data, marker=None, color=color, label=label)
        # plt.loglog(x_data, y_data, basex=10, color=color, label=label)
        #   add label to graph
        # plt.text(abs(min(self._distance_solution_container[_contact_indices])), max(y_data), str_id+" "+label)

    def _contact(self):
        """

        :return:
        """
        #   indices where contact is present, value is 1
        _contact_indices = np.nonzero(self._status_container)#_Fn_solution_container, _status_container
        return _contact_indices

    def _plot_info(self, color=None):
        """

        :return:
        """
        #   label
        label = self.contact_model

        #   plot
        if color is None:
            color = self._color_GL
        else:
            color = color
        return label, color

    def _write_to_txt_file(self):
        """
        Function saves data to .txt file (.dat, .sol)
        :return:
        """
        #   set headers
        self._headers = self._txt_headers()

        #   format headers for each column
        __frmt = "i-th step\tstatus\t\tdt\t\t\t\t\t\ttime \t\t\t\t\t delta \t\t\t\t\t dq_n \t\t\t\t\t Fn \t\t\t\t\t Ft"

        #   add new line at the end of comment
        self._comments = self._comments + "\n"

        #   save to file
        np.savetxt(self.filename, self.solution_data, fmt=['%5i', '%5i', '%20.10f', '%20.10f', '%20.10f', '%20.10f', '%20.10f', '%20.10f'], delimiter='\t', header=_header, comments = self._comments)


if __name__ == "__main__":
    # filename = "C:\\Users\\lskrinjar\\Dropbox\\DyS\\dynamic_systems\\0_7_3_0_contact_models\\all_contact_models\\cr=0.7\\contact_0_hertz.xlsx"
    filename = "pin_slot_clearance_joint_0.xlsx"
    _sol_obj = SolutionDataContact()
    _sol_obj.read_file(filename)
    # print _sol_obj._status_container
    # print _sol_obj._distance_solution_container
    # print _sol_obj._dqn_solution_container
    # print _sol_obj.solution_container_names

    # print getattr(_sol_obj, _sol_obj.solution_container_names[0])[0]
    for i in xrange(0, 10):
        print _sol_obj._step_num_solution_container[i],
        print _sol_obj._status_container[i],
        print _sol_obj._step_size_solution_container[i],
        print _sol_obj._t_solution_container[i],
        print _sol_obj._distance_solution_container[i],
        print _sol_obj._dqn_solution_container[i],
        print _sol_obj._dqt_solution_container[i],
        print _sol_obj._Fn_solution_container[i],
        print _sol_obj._Ft_solution_container[i],
        print _sol_obj._F_solution_container[i, :],
        print _sol_obj._u_P_solution_container[i, :],
        print _sol_obj._r_P_solution_container[i, :],
        print ""


    plt.plot(abs(_sol_obj._distance_solution_container[_sol_obj._contact()])[1::], _sol_obj._dqn_solution_container[_sol_obj._contact()][1::])
    plt.show()

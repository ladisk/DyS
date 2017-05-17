__author__ = 'lskrinjar'

import os

import matplotlib as mpl
import seaborn
from matplotlib import pyplot as plt

from MBD_system.fix_string import fix_string
from MBD_system.solution_data.solution_data import SolutionData
from MBD_system.solution_data.solution_data_contact import SolutionDataContact


class MBD_solution_analysis(object):
    """

    :return:
    """
    def __init__(self, abs_folder_path, file_extension):
        #   predefined empty list
        self.mbd_solutions = []
        for item in os.listdir(abs_folder_path):
        # for root, dirs, files in os.walk(abs_folder_path, topdown=False):
        #     for i, _file in enumerate(files):
            #   filename and extension
            filename, _file_extension = os.path.splitext(item)

            #   check extension
            if _file_extension == file_extension:
                _file = os.path.join(abs_folder_path, item)
                #   create solution data object to read and load solution data to object
                _sol_data = SolutionData(_file=_file)
                _sol_data.read_file()

                self.mbd_solutions.append(_sol_data)

    def plot(self, x, y, colors=None, labels=None):
        """

        """
        if colors is None:
            colors = [None] * len(self.mbd_solutions)

        if labels is None:
            labels = [None] * len(self.mbd_solutions)

        #    plot the data
        for i, (sol, color, label) in enumerate(zip(self.mbd_solutions, colors, labels)):
            if label is None:
                try:
                    _indx = sol._name.index('0')+3
                    label = sol._name[_indx::]
                    label = fix_string(label)
                except:
                    label = sol._name

            print label
            sol._plot(x, y, color=color, label=label)
                # plt.plot(sol._t_solution_container, sol._q_sol_container[:, 6], color=color, linestyle="-", label=_label)
                # plt.plot(sol._t_solution_container, sol._q_sol_container[:, 9], color=color, linestyle="--", label=_label)
        # Shink current axis by 20%
        # box = ax.get_position()
        # ax.set_position([box.x0, box.y0, box.width * 0.65, box.height])
        # plt.grid()
        # plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        # plt.show()


class Contact_solution_analysis(object):
    """

    :return:
    """
    def __init__(self, abs_folder_path, file_extension):

        #   predefined empty list
        self._solutions = []
        for item in os.listdir(abs_folder_path):
        # for root, dirs, files in os.walk(abs_folder_path, topdown=False):
        #     for i, _file in enumerate(files):
                #   filename and extension
            filename, _file_extension = os.path.splitext(item)
            #   check extension
            if _file_extension == file_extension:
                _file = os.path.join(abs_folder_path, item)
                #   create contact object
                _sol_obj = SolutionDataContact()
                #   read solution data from external file and save to object
                _sol_obj.read_file(_file)

                #   append contact object to list
                self._solutions.append(_sol_obj)

    def plot(self, x, y, colors=None):
        """

        :return:
        """
        if colors is None:
            colors = [None] * len(self._solutions)

        for i, (sol, color) in enumerate(zip(self._solutions, colors)):
            if sol.contact_model == "Kelvin-Voigt":
                pass
                # sol.plot_Fn_delta(color)
            else:
                # sol.plot_Fn_delta(color)
                # sol.plot_Fn_time(color)
                sol.plot(x, y, color=color, str_id=str(i))

if __name__ == "__main__":
    #   relative path to folder with data
    # _folder_rel_path = "dynamic_systems/0_7_3_0_contact_models/all_contact_models/lankarani-nikravesh_cr=0-1/mbd"
    _folder_rel_path = "dynamic_systems/0_7_3_0_contact_models/all_contact_models/cr=0.7"
    # _folder_rel_path = "dynamic_systems/0_7_3_0_contact_models/all_contact_models/cr=1"
    # _folder_rel_path = "dynamic_systems/0_7_3_0_contact_models/all_contact_models"
    # _folder_rel_path = "dynamic_systems/0_7_3_0_contact_models/all_contact_models/int_euler"
    # _folder_rel_path = "dynamic_systems/0_7_3_0_contact_models/all_contact_models/int_methods"
    # _folder_rel_path = "dynamic_systems/0_7_3_0_contact_models_cylinder/all_contact_models_cylinder/mbd"
    # _folder_rel_path = "dynamic_systems/0_7_3_0_contact_models_cylinder/all_contact_models_cylinder/contact"
    #   absolute path to project folder
    project_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd())))
    #   absolute path to folder with files
    _folder_abs_path = os.path.join(project_dir, _folder_rel_path).replace('\\', '/')

    # _file = "c:/Users/lskrinjar/Dropbox/DyS/DyS_project/MBD_system/contact/solution_data_contact.xlsx"
    solutions = Contact_solution_analysis(_folder_abs_path, ".xlsx")
    # solutions = MBD_solution_analysis(_folder_abs_path, ".sol")
    _colors = [[0, 0, 0],
               [.7, .3, .1],
               [1, 0, 0],
               [0, 1, 0],
               [0, 0, 1],
               [.5, .5, .5],
               [.2, 0.2, 0.2],
               [0, 0, .5],
               [.7, .3, .1],
               [0.1, 0.1, .9],
               [0.4, 0.1, 0.3]]

    _filename = 'contact_models_Fn_delta_cylinder.png'
    mpl.use("pgf")

    # mpl.rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern Sans serif']})
    # mpl.rc('font',**{'family':'serif','serif':['Times']})
    # mpl.rcParams['text.usetex'] = True #True, False

    mpl.rcParams['text.usetex'] = False
    seaborn.set_style(style='white')

    fig = plt.figure(num=1, figsize=(7, 5), dpi=100, facecolor='w', edgecolor='k')
    ax = plt.subplot(111, aspect="auto")
    ax.ticklabel_format(style='sci',scilimits=(-4,4), axis='both')

    #   tex settings

    # matplotlib.rcParams['text.latex.unicode'] = True

    # ax.xaxis.major.formatter._useMathText = True
    # cont.plot_Fn_delta()
    # x = "t"#delta, t
    # y = "energy"#dqn, Fn, dRx
    # solutions.plot(x, y, _colors)

    x = "delta"#delta, t
    y = "dqn"#dqn, Fn, dRx
    solutions.plot(x, y, _colors)

    # box = ax.get_position()
    # ax.set_position([box.x0, box.y0, box.width * 0.5, box.height])

    # Put a legend to the right of the current axis
    # ax.legend(loc='upper left', bbox_to_anchor=(1, 0.5))
    leg = ax.legend(loc='best', fontsize=12, frameon=True)#, bbox_to_anchor=(1, 0.5))
    leg.get_frame().set_alpha(1)
    plt.grid()
    # plt.legend(loc='upper right', bbox_to_anchor=(1, 1))
    # plt.show()
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.tick_params(axis='both', which='minor', labelsize=12)

    # plt.xlabel(r'\ $\delta \left[m\right]$', fontsize=14)
    # plt.ylabel(r'\ $\dot{\delta}_n \left[\frac{m}{s}\right]$', fontsize=14)
    # plt.ylabel(r'\ $\dot{R}$ [m/s]', fontsize=12)
    # plt.ylabel(r'\ $E$ [J]', fontsize=12)
    # plt.xlabel(r'\ t [s]', fontsize=12)
    # plt.ylabel(r'\ $F$ [N]', fontsize=12)
    # y_label.set_rotation(0)

    # plt.xlim([0.003333, 0.00336])
    # plt.ylim([0, 1.5])
    # plt.savefig(_filename, format='png')

    # pdf.close()
    plt.show()
__author__ = 'lskrinjar'

import os

from pprint import pprint
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib import rcParams
from matplotlib.backends.backend_pdf import PdfPages


from MBD_system.solution_data.solution_data import SolutionData
from MBD_system.contact.contact import Contact
from MBD_system.contact_model.analysis_of_contact_models import MBD_solution_analysis
from MBD_system.contact_model.analysis_of_contact_models import Contact_solution_analysis


if __name__ == "__main__":
    #   relative path to folder with data
    _folder_rel_path = "dynamic_systems/0_7_3_0_plane_sphere"
    #   absolute path to project folder
    project_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd())))
    #   absolute path to folder with files
    _folder_abs_path = os.path.join(project_dir, _folder_rel_path).replace('\\', '/')
    _folder_abs_path = "c://Users//lskrinjar//Dropbox//DyS//dynamic_systems//0_7_3_0_plane_sphere"
    # _file = "c:/Users/lskrinjar/Dropbox/DyS/DyS_project/MBD_system/contact/solution_data_contact.xlsx"
    # solutions = Contact_solution_analysis(_folder_abs_path, ".xlsx")
    solutions = MBD_solution_analysis(_folder_abs_path, ".sol")
    _colors = [[.6, .6, .6],
               [.1, .8, .1]]

    # _filename = 'contact_models_Fn_delta.pdf'
    # pdf = PdfPages(_filename)
    fig = plt.figure(num=1, figsize=(6, 4), dpi=100, facecolor='w', edgecolor='k')
    ax = plt.subplot(111, aspect="auto")
    ax.ticklabel_format(style='sci',scilimits=(-4,4), axis='both')

    #   tex settings
    rcParams['text.usetex'] = True

    rc('font', family='serif')
    # ax.xaxis.major.formatter._useMathText = True
    # cont.plot_Fn_delta()
    x = "Rx"#delta, t
    y = "Ry"#dqn, Fn, dRx
    solutions.plot(x, y, _colors, [r"\ brez kontrole $\delta_{0}$", r"\ $0 \leq \delta_{0} \leq \delta_{0;TOL}$"])

    # box = ax.get_position()
    # ax.set_position([box.x0, box.y0, box.width * 0.5, box.height])

    # Put a legend to the right of the current axis
    # ax.legend(loc='upper left', bbox_to_anchor=(1, 0.5))
    ax.legend(loc='best', fontsize=12)#, bbox_to_anchor=(1, 0.5))
    plt.grid()
    # plt.legend(loc='upper right', bbox_to_anchor=(1, 1))
    # plt.show()
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.tick_params(axis='both', which='minor', labelsize=12)

    # plt.xlabel(r'\ $\delta$ [m]', fontsize=12)
    # plt.ylabel(r'\ $dq_n$ [m/s]', fontsize=12)
    plt.xlabel(r'\ x[m]', fontsize=12)
    plt.ylabel(r'\ y[m]', fontsize=12)

    # plt.ylabel(r'\ $F$ [N]', fontsize=12)
    # y_label.set_rotation(0)

    plt.xlim([0, 0.03])
    # plt.ylim([0, 1.5])
    # plt.savefig(pdf, format='pdf')

    # pdf.close()
    plt.show()
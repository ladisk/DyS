__author__ = 'lskrinjar'

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.backends.backend_pdf import PdfPages

from contact_model import ContactModel
from MBD_system.fix_string import fix_string

models = ["lankarani-nikravesh", "hunt-crossley", "herbert-mcwhannell", "lee-wang", "flores et al", "zhiying-qishao", "gonthier et al", "hu-guo"]

dcr = 1E-3
c_r = np.arange(0+dcr, 1+dcr, dcr)
K = 1
_dq0_n = 1
contact_models = []
for i, model in enumerate(models):
    print i, model
    contact_model = ContactModel(_type=model)
    contact_model.ksi_solution_data = []
    contact_model.K = K
    contact_model._dq0_n = _dq0_n

    for j in range(0, len(c_r)):
        contact_model.c_r = c_r[j]
        _ksi = contact_model.hysteresis_damping_factor()
        contact_model.ksi_solution_data.append(_ksi)

    contact_models.append(contact_model)


_colors = [[0, 0, 0],
               [.7, .3, .1],
               [1, 0, 0],
               [0, 1, 0],
               [0, 0, 1],
               [.5, .5, .5],
               [.2, 0.2, 0.2],
               [0, 0, .5],
               [.7, .3, .1],
               [0.1, 0.1, .9]]

sb.set_style(style='white')
fig = plt.figure(num=1, figsize=(8, 6), dpi=100, facecolor='w', edgecolor='k')
ax = plt.subplot(111, aspect="auto")
ax.ticklabel_format(style='sci',scilimits=(-4,4), axis='both')

#   tex settings
mpl.rcParams['text.usetex'] = True
mpl.rc('font', family='serif')



for contact_model, _color in zip(contact_models, _colors):
    _label = fix_string(contact_model._type)
    plt.plot(c_r, contact_model.ksi_solution_data, label=_label, color=_color)

leg = ax.legend(loc='best', fontsize=12, frameon=True)#, bbox_to_anchor=(1, 0.5))
leg.get_frame().set_alpha(1)
plt.grid()
# plt.legend(loc='upper right', bbox_to_anchor=(1, 1))
# plt.show()
plt.tick_params(axis='both', which='major', labelsize=12)
plt.tick_params(axis='both', which='minor', labelsize=12)

# plt.xlabel(r'\ $\delta$ [m]', fontsize=12)
# plt.ylabel(r'\ $dq_n$ [m/s]', fontsize=12)
plt.ylabel(r'\ $\chi$', fontsize=18)
plt.xlabel(r'\ $c_r$', fontsize=18)
# plt.ylabel(r'\ $F$ [N]', fontsize=12)
# y_label.set_rotation(0)

# plt.xlim([0, 20])
plt.ylim([0, 20])
# plt.savefig(pdf, format='pdf')

# pdf.close()
plt.show()

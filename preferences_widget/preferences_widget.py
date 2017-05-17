"""

created by: lskrinjar
date of creation: 28/01/2016
time of creation: 14:00
"""
import numpy as np
from PyQt4 import QtGui, QtCore

from MBD_system.array2string import array2string
from preferences_widget_ui import Ui_Form

class PreferencesWidget(QtGui.QWidget):
    """
    classdocs
    """
    def __init__(self, simulation_control_widget, parent=None):
        """
        Constructor
        """
        super(PreferencesWidget, self).__init__(parent=parent)

        self._parent = parent

        self.ui = Ui_Form()
        self.ui.setupUi(self)

        #   pointer to simnlation control widget as attribute
        self.simulation_control_widget = simulation_control_widget

        #   pointer to open gl widget as attribute
        self.vtkWidget = self.simulation_control_widget.vtkWidget

        #   update display type
        self.update_display_type = self.ui.updateDisplay_comboBox.currentText()
        
        #    force scale
        self._force_scale = 1

        #   signals
        #   list through different widgets in stacked widget
        self.ui.listWidget.currentRowChanged.connect(self._options)
        #   close when canceled is clicked
        self.ui.buttonBox.button(QtGui.QDialogButtonBox.Cancel).clicked.connect(self._cancel)
        #   save options
        self.ui.buttonBox.button(QtGui.QDialogButtonBox.Save).clicked.connect(self._save)
        #   focus on save button
        self.ui.buttonBox.button(QtGui.QDialogButtonBox.Save).setFocus()
        #   edit
        self.ui.editBackgroundColor_pushButton.clicked.connect(self._edit_background_color)

    def _edit_background_color(self):
        """

        :return:
        """
        color = QtGui.QColorDialog.getColor()
        if color.isValid():
            rgb = color.getRgbF()

            if np.linalg.norm(np.array(rgb)) > 1.2:
                text_color = np.array([0., 0., 0.])

            else:
                text_color = np.array([1., 1., 1.])

            self.vtkWidget.renderer.SetBackground(rgb[0], rgb[1], rgb[2])
            self.vtkWidget.setTextColor(text_color)
            self.vtkWidget.update()

    def _cancel(self):
        """

        """
        self.close()

    def _options(self, i):
        """

        :param i:
        :return:
        """
        self.ui.stackedWidget.setCurrentIndex(i)

    def _save(self):
        """

        :return:
        """
        # self.opengl_widget._show_filename = self.ui.filename_checkBox.checkState()
        # self.opengl_widget._show_simulationTime = self.ui.simulationTime_checkBox.checkState()
        # self.opengl_widget._show_simulationStepNumber = self.ui.simulationStepNumber_checkBox.checkState()
        # self.opengl_widget._show_timeAndDate = self.ui.timeAndDate_checkBox.checkState()

        #   when everything is saved close widget
        self.close()

    def _show(self):
        """

        :return:
        """
        #   set stacked widget to 0
        self.ui.stackedWidget.setCurrentIndex(0)

        #   main
        self._show_main()

        #   simulation
        self._show_simulation()

        #   visualization
        self._show_visualization()

        #   when all attributes are read display widget
        self.show()

    def _show_main(self):
        """

        :return:
        """
        #   show GCS
        # self.ui.showGCS_checkBox.setChecked(self.opengl_widget.GCS._visible)

    def _show_simulation(self):
        """

        :return:
        """

    def _show_visualization(self):
        """

        :return:
        """
        #   background color
        # self.ui.backgroundColor_lineEdit.setText(array2string(self.opengl_widget.background_color[0:3]))

        #   display info checkboxes - removed first item in qgroupbox (it is qlayout)
        # checkboxes = iter(self.ui.infoDisplay_groupBox.children())
        # next(checkboxes)
        # for checkBox, state in zip(checkboxes, self.vtkWidget._show_info):
        #     # checkBox = self.ui.infoDisplay_groupBox.children().itemAt(i)
        #     checkBox.setChecked(state)
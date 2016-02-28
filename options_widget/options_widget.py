"""

created by: lskrinjar
date of creation: 28/01/2016
time of creation: 14:00
"""
from PyQt4 import QtGui

from MBD_system.array2string import array2string
from options_widget_ui import Ui_Form

class OptionsWidget(QtGui.QWidget):
    """
    classdocs
    """
    def __init__(self, simulation_control_widget, parent=None):
        """
        Constructor
        """
        super(OptionsWidget, self).__init__(parent=parent)

        self._parent = parent

        self.ui = Ui_Form()
        self.ui.setupUi(self)

        #   pointer to simnlation control widget as attribute
        self.simulation_control_woidget = simulation_control_widget

        #   pointer to open gl widget as attribute
        self.opengl_widget = self.simulation_control_woidget.opengl_widget

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
            self.opengl_widget.qglClearColor(color)
            self.opengl_widget.updateGL()

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
        self.opengl_widget._show_filename = self.ui.filename_checkBox.checkState()
        self.opengl_widget._show_simulationTime = self.ui.simulationTime_checkBox.checkState()
        self.opengl_widget._show_simulationStepNumber = self.ui.simulationStepNumber_checkBox.checkState()
        self.opengl_widget._show_timeAndDate = self.ui.timeAndDate_checkBox.checkState()

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
        self.ui.showGCS_checkBox.setChecked(self.opengl_widget.GCS._visible)

    def _show_simulation(self):
        """

        :return:
        """

    def _show_visualization(self):
        """

        :return:
        """
        #   opengl background color
        self.ui.backgroundColor_lineEdit.setText(array2string(self.opengl_widget.background_color[0:3]))

        #   display info checkboxes - removed first item in qgroupbox (it is qlayout)
        checkboxes = iter(self.ui.infoDisplay_groupBox.children())
        next(checkboxes)
        for checkBox, state in zip(checkboxes, self.opengl_widget._show_info):
            # checkBox = self.ui.infoDisplay_groupBox.children().itemAt(i)
            checkBox.setChecked(state)
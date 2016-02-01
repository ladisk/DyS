"""

created by: lskrinjar
date of creation: 28/01/2016
time of creation: 14:00
"""
import os
from pprint import pprint
from PyQt4 import QtCore, QtGui


from options_widget_ui import Ui_Form


class OptionsWidget(QtGui.QWidget):
    """
    classdocs
    """
    def __init__(self, parent=None):
        """
        Constructor
        """
        super(OptionsWidget, self).__init__(parent=parent)

        self._parent = parent

        self.ui = Ui_Form()
        self.ui.setupUi(self)

        #todo - to create functional preferences widget see qstackwidget in sandbox/pyqt/qstackedwidget/example_2.py

        #    signals
        self.ui.buttonBox.button(QtGui.QDialogButtonBox.Cancel).clicked.connect(self._cancel)

    def _cancel(self):
        """

        """
        self.close()
'''
Created on 12. okt. 2015

@author: luka.skrinjar
'''
import os
from pprint import pprint
from PyQt4 import QtCore, QtGui
from MBD_system_ui import Ui_Form
from array2string import array2string
from string2array import string2array

class MBDSystemWidget(QtGui.QWidget):
    """
    classdocs
    """
    def __init__(self, MBD_system, parent=None):
        """
        Constructor
        """
        super(MBDSystemWidget, self).__init__(parent=parent)
        self._parent = parent

        self.ui = Ui_Form()
        self.ui.setupUi(self)

        self.setWindowFlags(self.windowFlags() | QtCore.Qt.Window)
        self.setParent(self._parent)
        self.setWindowModality(QtCore.Qt.WindowModal)

        self.MBD_system = MBD_system
        
        #    show widget
        self.raise_()

        #   signals
        self.ui.cancel_pushButton.clicked.connect(self.close)
        self.ui.save_pushButton.clicked.connect(self._save)
        self.ui.save_pushButton.setFocus()
        
    def _edit(self, item):
        """
        
        """
        print "item =", item
        self.item = item

        #    filename
        filename, filetype = os.path.splitext(self.MBD_system.filename_)
        self.ui.filename_lineEdit.setText(str(filename))
        
        #    project folder
        self.ui.project_folder_lineEdit.setText(str(self.MBD_system.MBD_folder_abs_path_))

        #   solution filetype
        _index = self.ui.solutionFiletype_comboBox.findText(self.item._solution_filetype)
        self.ui.solutionFiletype_comboBox.setCurrentIndex(_index)

        #    gravity - unit vector
        self.ui.gravity_vector_lineEdit.setText(array2string(self.MBD_system.gravity_vector))
        
        #    gravity - magnitude
        self.ui.gravity_value_lineEdit.setText(str(self.MBD_system.gravity))

        #   display GCS
        if self.item.GCS_visible: #self.ui.showGCS_checkBox.isChecked():
            self.ui.showGCS_checkBox.setChecked(QtCore.Qt.Checked)
        else:
            self.ui.showGCS_checkBox.setChecked(QtCore.Qt.UnChecked)

        self.show()

    def _save(self):
        """

        :return:
        """
        print "self.item =", self.item

        self.item.GCS_visible = self.ui.showGCS_checkBox.isChecked()

        #   save solution filetype
        self.item._solution_filetype = self.ui.solutionFiletype_comboBox.currentText()

        self.close()
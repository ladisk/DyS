"""
Created on 3. mar. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
"""
import os
from pprint import pprint

import numpy as np
from PyQt4 import QtCore, QtGui

from MBD_system.Ai_ui_P import Ai_ui_P_vector
from MBD_system.array2string import array2string
from MBD_system.body.body import Body
from MBD_system.string2array import string2array
from body_ui import Ui_Form

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s


class BodyWidget(QtGui.QWidget):#QtGui.QDialog
    """
    control panel interface
    """
    
    def __init__(self, group_item, parent=None):
        """
        Constructor
        """
        super(BodyWidget, self).__init__(parent=parent)
        self._parent = parent

        self.ui = Ui_Form()
        self.ui.setupUi(self)

        # self.setFocusPolicy(QtCore.Qt.NoFocus)#self.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint) #self.setFocusPolicy(QtCore.Qt.NoFocus)
        # self.setFocusPolicy(QtCore.Qt.StrongFocus)
        # self.setWindowFlags(self.windowFlags() | QtCore.Qt.WindowStaysOnTopHint | QtCore.Qt.FramelessWindowHint)
        self.setWindowFlags(self.windowFlags() | QtCore.Qt.Window)
        self.setParent(self._parent)
        self.setWindowModality(QtCore.Qt.WindowModal)
        # widget.show()

        self.group_item = group_item
        self.item = None
        
        self.ui.save_pushButton.setAutoDefault(True)

        #    set validators to limit user input
        __validator_dbl = QtGui.QDoubleValidator()
        __validator_int = QtGui.QIntValidator()
        
        self.ui.R_lineEdit.setValidator(QtGui.QDoubleValidator())
        
        #    body id
        self.ui.bodyID_lineEdit.setValidator(__validator_int)

        #    signals
        self.ui.cancel_pushButton.clicked.connect(self._cancel)
        self.ui.save_pushButton.clicked.connect(self._save)
        self.ui.save_pushButton.setFocus()

        self.ui.load_geometry_pushButton.clicked.connect(self._load_stl_file)
        self.ui.load_data_pushButton.clicked.connect(self._load_dat_file)


        self.ui.bodyID_lineEdit.setText(str(len(self.group_item._children)))
        
        
        self._move()

        #    show widget
        self.raise_()
        self.show()
        
    def closeEvent(self, event):
        """

        :param event:
        :return:
        """
        self._parent._parent.simulation_control_widget.opengl_widget._repaintGL()

    def _move(self):
        """
        Move widget (dialog window) to screen center
        """
        frameGm = self.frameGeometry()
        screen = QtGui.QApplication.desktop().screenNumber(QtGui.QApplication.desktop().cursor().pos())
        centerPoint = QtGui.QApplication.desktop().screenGeometry(screen).center()
        frameGm.moveCenter(centerPoint)
        self.move(frameGm.topLeft())

    def _cancel(self):
        """
        
        """
        self.close()

    def _save(self):
        """
        
        """
        _name = str(self.ui.name_lineEdit.text())
        
#         try:
        _mass = float(self.ui.m_lineEdit.text())
        _J_zz = float(self.ui.J_lineEdit.text())
        
        _R = np.array(string2array(self.ui.R_lineEdit.text()), dtype="float32")
        _theta = np.array(string2array(self.ui.theta_lineEdit.text()), dtype="float32")
        
        
        
        if self.ui.transformCS_comboBox.currentText() == "CAD":
            __dR = self.item.CM_CAD_LCS - Ai_ui_P_vector(self.item.CM_CAD_LCS, np.deg2rad(_theta[2]))
        else:
            __dR = np.zeros(3)

        _dR = string2array(self.ui.dR_lineEdit.text())
        _dtheta = string2array(self.ui.dtheta_lineEdit.text())
        
        _color = string2array(self.ui.color_lineEdit.text())
        _transparent = float(self.ui.transparent_lineEdit.text())
        _display_style = str(self.ui.display_style_comboBox.currentText()).lower()

        #    update data to selected object
        if self.item is not None:
            self.item._name = _name
            self.item.mass = _mass
            self.item.J_zz = _J_zz
            self.item.R = _R - __dR
            self.item.theta = np.deg2rad(_theta)
            self.item.dR = _dR
            self.item.dtheta = np.deg2rad(_dtheta)
            self.item.display_style = _display_style

            #   update vbo if it is changed
            if any(_color != self.item.color_GL) or _transparent != self.item.transparent_GL:
                self.item.color_GL = _color
                self.item.transparent_GL = _transparent
                self.item._update_VBO()

        #    create new object
        else:
            _item = Body(body_name=_name, parent=self.group_item)
            
            pos = len(self.parent_node._children)
            self.parent_node._parent.forces.append(_item)
            self._parent.ui.treeView.model().insertRow(pos, _item, self.parent_node)        

        self.close()
#         except:
#             QtGui.QMessageBox.warning(self, "Warning!",
#                 "Input not correct!",
#                 QtGui.QMessageBox.Cancel, QtGui.QMessageBox.NoButton,
#                 QtGui.QMessageBox.NoButton)

    def _load_stl_file(self):
        """

        :return:
        """
        self._load_file_ = QtGui.QFileDialog()
        self._load_file_.setDirectory(os.path.curdir)
        self.stl_filename, self.file_type = self._load_file_.getOpenFileNameAndFilter(self, 'Open file', os.path.curdir, ("STL (*.stl)"))

    def _load_dat_file(self):
        """

        :return:
        """
        self._load_file_ = QtGui.QFileDialog()
        self._load_file_.setDirectory(os.path.curdir)
        self.dat_filename, self.file_type = self._load_file_.getOpenFileNameAndFilter(self, 'Open file', os.path.curdir, ("dat (*.dat)"))
        
    def _edit(self, item):
        """
        Set object data to display in qwidget elements
        :return:
        """
        self.item = item
        pprint(vars(item))

        #   check if attribute is type QVariant and change it to string
        if type(item._name) == QtCore.QVariant:
            self.ui.name_lineEdit.setText(item._name.toString())
        else:
            self.ui.name_lineEdit.setText(item._name)

        self.ui.bodyID_lineEdit.setText(str(item.body_id))

        self.ui.m_lineEdit.setText(str(item.mass))

        self.ui.J_lineEdit.setText(str(item.J_zz))

        self.ui.R_lineEdit.setText(array2string(self.item.R))

        self.ui.dR_lineEdit.setText(array2string(item.dR))

        self.ui.theta_lineEdit.setText(array2string(np.rad2deg(item.theta)))

        self.ui.dtheta_lineEdit.setText(array2string(np.rad2deg(item.dtheta)))

        self.ui.centerOfMass_lineEdit.setText(array2string(item.CM_CAD_LCS))

        self.ui.geometrySTLfile_lineEdit.setText(os.path.abspath(item.geometry_filename))

        self.ui.color_lineEdit.setText(array2string(item.color_GL))

        self.ui.transparent_lineEdit.setText(str(item.transparent_GL))
        
        _index = self.ui.display_style_comboBox.findText(QtCore.QString(item.display_style.capitalize()))
        self.ui.display_style_comboBox.setCurrentIndex(_index)

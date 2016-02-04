"""
Created on 3. mar. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
"""
import time
from pprint import pprint

import numpy as np
from PyQt4 import QtCore, QtGui


from joint_ui import Ui_Form
from MBD_system.joint.joint import Joint
from MBD_system.array2string import array2string
from MBD_system.string2array import string2array


try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s


class JointWidget(QtGui.QDialog):
    """
    control panel interface
    """
    
    def __init__(self, parent_node=None, parent=None):
        """
        Constructor
        """
        super(JointWidget, self).__init__(parent)
        self._parent = parent


        self.ui = Ui_Form()
        self.ui.setupUi(self)
        
        self.parent_node = parent_node
        
        
        self.ui.save_pushButton.setAutoDefault(True)

        #    set validators to limit user input
        __validator_dbl = QtGui.QDoubleValidator()
        __validator_int = QtGui.QIntValidator()
        
        
        #    validators
        self.ui.jointID_lineEdit.setValidator(__validator_int)
        # self.ui.bodyID_lineEdit.setValidator(__validator_int)
        # self.ui.uPi_lineEdit.setValidator(__validator_dbl)
        # self.ui.Fx_lineEdit.setValidator(__validator_dbl)
        # self.ui.Fy_lineEdit.setValidator(__validator_dbl)
        # self.ui.Mz_lineEdit.setValidator(__validator_dbl)

        
        
        #    signals
        self.ui.cancel_pushButton.clicked.connect(self._cancel)
        self.ui.save_pushButton.clicked.connect(self._save)
        self.ui.save_pushButton.setFocus()
        
        self.ui.bodyIDi_lineEdit.textChanged.connect(self._check_if_ground)
        self.ui.bodyIDj_lineEdit.textChanged.connect(self._check_if_ground)
        
        
        self.ui.jointID_lineEdit.setText(str(len(self.parent_node._children)))
        
        
        #   move to center of screen
        self._move()

        #    show widget if new item is being created
        if parent_node._typeInfo == "group":
            self._show()


    def _check_if_ground(self):
        """
        
        """
        if self.ui.bodyIDi_lineEdit.text() == "ground" or self.ui.bodyIDi_lineEdit.text() == "-1":
            self.ui.uPi_lineEdit.setEnabled(False)
        else:
            self.ui.uPi_lineEdit.setEnabled(True)


        if self.ui.bodyIDj_lineEdit.text() == "ground" or self.ui.bodyIDj_lineEdit.text() == "-1":
            self.ui.uPj_lineEdit.setEnabled(False)
        else:
            self.ui.uPj_lineEdit.setEnabled(True)


    def _show(self):
        """

        :return:
        """
        self.show()


    def _move(self):
        frameGm = self.frameGeometry()
        screen = QtGui.QApplication.desktop().screenNumber(QtGui.QApplication.desktop().cursor().pos())
        centerPoint = QtGui.QApplication.desktop().screenGeometry(screen).center()
        frameGm.moveCenter(centerPoint)
        self.move(frameGm.topLeft())


    def _cancel(self):
        """
        
        """
        self.close()

    
    def _save(self, item=None):
        """
        
        """
        _name = self.ui.name_lineEdit.text()
        _type = self.ui.jointTypecomboBox.currentText()

        # try:
        if self.ui.bodyIDi_lineEdit.text() == "ground":
            pass
        else:
            _body_id_i = int(self.ui.bodyIDi_lineEdit.text())

        if self.ui.bodyIDj_lineEdit.text() == "ground":
            _body_id_j = self.ui.bodyIDj_lineEdit.text()
        else:
            _body_id_j = int(self.ui.bodyIDj_lineEdit.text())


        _uPi = string2array(self.ui.uPi_lineEdit.text())
        _uPj = string2array(self.ui.uPj_lineEdit.text())


        #    update data to selected object
        if self.item is not None:
            self.item._name = _name
            self.item.body_id_i = _body_id_i
            self.item.body_id_j = _body_id_j

            #   update position of marker if vector uP is changed
            for _uP, _u_P, _id, body_id in zip([_uPi, _uPj], [self.item.u_iP, self.item.u_jP], [0, 1], self.item.body_id_list):

                if (_uP != _u_P).any() and (body_id != "ground"):
                    self.item.markers[_id]._update_node(np.array(np.append(_uP, self.item.z_dim), dtype='float32'))



        else:
            _type = str(_type).lower()
            _item = Joint(_type, _body_id_i, _body_id_j)

            pos = len(self.parent_node._children)
            self.parent_node._parent.forces.append(_item)
            self._parent.ui.treeView.model().insertRow(pos, _item, self.parent_node)

        self.close()
        
        # except:
        #     QtGui.QMessageBox.warning(self, "Warning!",
        #         "Input not correct!",
        #         QtGui.QMessageBox.Cancel, QtGui.QMessageBox.NoButton,
        #         QtGui.QMessageBox.NoButton)


    def _edit(self, item=None):
        """

        :return:
        """
        self.ui.name_lineEdit.setText(item._name)

        self.ui.jointID_lineEdit.setText(str(item.joint_id))

        #   set type line edit
        _index = self.ui.jointTypecomboBox.findText(QtCore.QString(item.joint_type.capitalize()))
        self.ui.jointTypecomboBox.setCurrentIndex(_index)

        #   body id i
        self.ui.bodyIDi_lineEdit.setText(str(item.body_id_i))
        #   body id j
        self.ui.bodyIDj_lineEdit.setText(str(item.body_id_j))

        #   uPi
        self.ui.uPi_lineEdit.setText(array2string(item.u_iP_LCS))
        #   uPj
        self.ui.uPj_lineEdit.setText(array2string(item.u_jP_LCS))

        if item.joint_type == "prismatic":
            self.ui.uQi_lineEdit.setText(array2string(item.u_iQ_LCS))


        self.item = item

        self._show()
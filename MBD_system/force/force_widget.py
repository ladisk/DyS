"""
Created on 28. jul. 2015

@author: lskrinjar (email: skrinjar.luka@gmail.com)
"""
import datetime
from pprint import pprint

import numpy as np
from PyQt4 import QtCore, QtGui


from force_ui import Ui_Form
from MBD_system.force.force import Force
from MBD_system.array2string import array2string
from MBD_system.string2array import string2array

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s


class ForceWidget(QtGui.QDialog):
    """
    control panel interface
    """
    
    def __init__(self, parent_node=None, parent=None):
        """
        Constructor
        """
        super(ForceWidget, self).__init__(parent)
        self._parent = parent


        self.ui = Ui_Form()
        self.ui.setupUi(self)
        
        self.parent_node = parent_node
        
        
        self.ui.save_pushButton.setAutoDefault(True)

        #    set validators to limit user input
        __validator_dbl = QtGui.QDoubleValidator()
        __validator_int = QtGui.QIntValidator()
        
        
        #    validators
        self.ui.forceID_lineEdit.setValidator(__validator_int)
        self.ui.bodyID_lineEdit.setValidator(__validator_int)
        self.ui.uPi_lineEdit.setValidator(__validator_dbl)
#         self.ui.Fx_lineEdit.setValidator(__validator_dbl)
#         self.ui.Fy_lineEdit.setValidator(__validator_dbl)
#         self.ui.Mz_lineEdit.setValidator(__validator_dbl)

        
        
        #    signals
        self.ui.cancel_pushButton.clicked.connect(self._cancel)
        self.ui.save_pushButton.clicked.connect(self._save)
        self.ui.save_pushButton.setFocus()
        
        
        # pprint(vars(self.group_item))
        self.ui.forceID_lineEdit.setText(str(len(self.parent_node._children)))
        
        #   move to center of screen
        self._move()

        #    show widget if new item is being created
        if parent_node._typeInfo == "group":
            self._show()
        #    signals
#         self.setAttribute(QtCore.Qt.WA_DeleteOnClose)


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
        _force_name = self.ui.name_lineEdit.text()
        # try:
        _body_id = int(self.ui.bodyID_lineEdit.text())


        try:
            _Fx = float(self.ui.Fx_lineEdit.text())
        except:
            _Fx = self.ui.Fx_lineEdit.text()
        
        try:
            _Fy = float(self.ui.Fy_lineEdit.text())
        except:
            _Fy = self.ui.Fy_lineEdit.text()
        
        try:
            _Mz = float(self.ui.Mz_lineEdit.text())
        except:
            _Mz = self.ui.Mz_lineEdit.text()
            
        _u_iP_f = string2array(self.ui.Mz_lineEdit.text())

        if self.item is not None:
            self.item.body_id = _body_id
            self.item.force_name = _force_name
            self.item.Fx = _Fx
            self.item.Fy = _Fy
            self.item.Mz = _Mz
            self.item.u_iP_f = _u_iP_f
        else:
            _force = Force(body_id=_body_id, force_name=_force_name, Fx=_Fx, Fy=_Fy, Mz=_Mz, u_iP_f=_u_iP_f)

            pos = len(self.parent_node._children)

            self.parent_node._parent.forces.append(_force)
            self._parent.ui.treeView.model().insertRow(pos, _force, self.parent_node)

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
        self.ui.forceID_lineEdit.setText(str(item.force_id))
        self.ui.bodyID_lineEdit.setText(str(item.body_id))
        self.ui.uPi_lineEdit.setText(array2string(item.u_iP_f))
        self.ui.Fx_lineEdit.setText(str(item.Fx))
        self.ui.Fy_lineEdit.setText(str(item.Fy))
        self.ui.Mz_lineEdit.setText(str(item.Mz))


        self.item = item

        self._show()
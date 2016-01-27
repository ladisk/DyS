'''
Created on 28. jun. 2015

@author: lskrinjar (email: skrinjar.luka@gmail.com)
'''
import datetime
from pprint import pprint

import numpy as np
from PyQt4 import QtCore, QtGui, uic, Qt


from contact_ui import Ui_Form
from MBD_system.contact.contact import Contact
from MBD_system.array2string import array2string
from MBD_system.string2array import string2array

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s


class ContactWidget(QtGui.QDialog):
    """
    control panel interface
    """
    
    def __init__(self, parent_node=None, parent=None):
        '''
        Constructor
        '''
        super(ContactWidget, self).__init__(parent)
        self._parent = parent


        self.ui = Ui_Form()
        self.ui.setupUi(self)
        
        self.parent_node = parent_node
        
        
        self.ui.save_pushButton.setAutoDefault(True)

        #    set validators to limit user input
        __validator_dbl = QtGui.QDoubleValidator()
        __validator_int = QtGui.QIntValidator()
        
        
        #    validators
        self.ui.contactID_lineEdit.setValidator(__validator_int)

        #   body id list at custom context menu
        # self.ui.bodyIDi_lineEdit.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        # self.ui.bodyIDi_lineEdit.customContextMenuRequested.connect(self._contextMenuEvent)
        # self.ui.uPi_lineEdit.setValidator(__validator_dbl)
        # self.ui.Fx_lineEdit.setValidator(__validator_dbl)
        # self.ui.Fy_lineEdit.setValidator(__validator_dbl)
        # self.ui.Mz_lineEdit.setValidator(__validator_dbl)

        #   validators for line edit user input
        self.ui.z_dim_lineEdit.setValidator(__validator_dbl)
        
        
        #    signals
        self.ui.cancel_pushButton.clicked.connect(self._cancel)
        self.ui.save_pushButton.clicked.connect(self._save)
        self.ui.save_pushButton.setFocus()
        self.ui.contactTypecomboBox.currentIndexChanged.connect(self._type_changed)
        
        
        # pprint(vars(self.group_item))
        self.ui.contactID_lineEdit.setText(str(len(self.parent_node._children)))
        
        #   move to center of screen
        self._move()

        #    show widget if new item is being created
        if parent_node is not None:
            self._show()
        #    signals
#         self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
    def _contextMenuEvent(self, event):
        """

        :return:
        """
        menu = QtGui.QMenu(self)
        pprint(vars(self._parent))
        for body in self._parent.bodies:
            menu.addAction(str(body.body_id)+"("+body._name+")")

        menu.exec_(event.globalPos())


    def _type_changed(self, _index):
        """

        :param type:
        :return:
        """
        _type = self.ui.contactTypecomboBox.itemText(_index)

        if _type == "General":
            self.ui.uPi_lineEdit.setEnabled(False)
            self.ui.uPj_lineEdit.setEnabled(False)
            self.ui.Ri_lineEdit.setEnabled(False)
            self.ui.Rj_lineEdit.setEnabled(False)

        if _type == "Revolute Clearance Joint":
            self.ui.uPi_lineEdit.setEnabled(True)
            self.ui.uPj_lineEdit.setEnabled(True)
            self.ui.Ri_lineEdit.setEnabled(True)
            self.ui.Rj_lineEdit.setEnabled(True)


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
        _type = self.ui.contactTypecomboBox.currentText()

#         try:
        _body_id_i = int(self.ui.bodyIDi_lineEdit.text())
        _body_id_j = int(self.ui.bodyIDj_lineEdit.text())
        
        _coef_of_friction_static = float(self.ui.frictionCoefStatic_lineEdit.text())
        _coef_of_friction_dynamic = float(self.ui.frictionCoefDynamic_lineEdit.text())
        _coef_of_restitution = float(self.ui.restitutionCoef_lineEdit.text())
        
        if self.item._type == "Revolute Clearance Joint":
            _uPi = string2array(self.ui.uPi_lineEdit.text())
            _uPj = string2array(self.ui.uPj_lineEdit.text())
            
            _Ri = float(self.ui.Ri_lineEdit.text())
            _Rj = float(self.ui.Rj_lineEdit.text())

        if self.item._type == "General":
            pass
            

        if self.item is not None:
            self.item.body_id_i = _body_id_i
            self.item.body_id_j = _body_id_j

            self.item.friction_model.coef_of_friction_static = _coef_of_friction_static
            self.item.friction_model.coef_of_friction_dynamic = _coef_of_friction_dynamic
            self.item.contact_model.coef_of_restitution = _coef_of_restitution
            
            
            if self.item._type == "Revolute Clearance Joint":
                self.item.u_iP = _uPi
                self.item.u_jP = _uPj
                self.item.R_i = _Ri
                self.item.R_j = _Rj
            
            
        else:
            _item = Contact(body_id_i=_body_id_i, body_id_j=_body_id_j)

            pos = len(self.parent_node._children)
            self.parent_node._parent.forces.append(_item)
            self._parent.ui.treeView.model().insertRow(pos, _item, self.parent_node)

        self.close()

#         except:
#             QtGui.QMessageBox.warning(self, "Warning!",
#                 "Input not correct!",
#                 QtGui.QMessageBox.Cancel, QtGui.QMessageBox.NoButton,
#                 QtGui.QMessageBox.NoButton)


    def _edit(self, item):
        """

        :return:
        """
        pprint(vars(item))
        self.ui.name_lineEdit.setText(item._name)

        self.ui.contactID_lineEdit.setText(str(item.contact_id))

        #   set type line edit
        _index = self.ui.contactTypecomboBox.findText(QtCore.QString(item._type.title()))
        self.ui.contactTypecomboBox.setCurrentIndex(_index)

        #   body id i
        self.ui.bodyIDi_lineEdit.setText(str(item.body_id_i))
        #   body id j
        self.ui.bodyIDj_lineEdit.setText(str(item.body_id_j))

        #    contact model
        _index = self.ui.contactModelTypecomboBox.findText(QtCore.QString(item.contact_model_type.title()), flags=QtCore.Qt.MatchCaseSensitive)#MatchCaseSensitive
        print "_index =", _index
        self.ui.contactModelTypecomboBox.setCurrentIndex(_index)
        
        #    friction model
        _index = self.ui.frictionModelTypecomboBox.findText(QtCore.QString(item.friction_model_type.title()), flags=QtCore.Qt.MatchContains)#MatchCaseSensitive
        self.ui.frictionModelTypecomboBox.setCurrentIndex(_index)
        
        #    coefficient of restitution
        self.ui.restitutionCoef_lineEdit.setText(str(item.contact_model.coef_of_restitution))
        
        #    coefficient of friction static
        self.ui.frictionCoefStatic_lineEdit.setText(str(item.friction_model.coef_of_friction_static))
        
        #    coefficient of friction dynamic
        self.ui.frictionCoefDynamic_lineEdit.setText(str(item.friction_model.coef_of_friction_dynamic))

        #   z coordinate of position of contact plane
        self.ui.z_dim_lineEdit.setText(str(item.z_dim))
        

        if item._type.capitalize() == "General":
            print "exe???"
            self.ui.uPi_lineEdit.setEnabled(False)
            self.ui.uPj_lineEdit.setEnabled(False)
            self.ui.Ri_lineEdit.setEnabled(False)
            self.ui.Rj_lineEdit.setEnabled(False)
        else:
            self.ui.uPi_lineEdit.setEnabled(True)
            self.ui.uPj_lineEdit.setEnabled(True)
            self.ui.Ri_lineEdit.setEnabled(True)
            self.ui.Rj_lineEdit.setEnabled(True)

            self.ui.uPi_lineEdit.setText(array2string(item.u_iP))
            self.ui.uPj_lineEdit.setText(array2string(item.u_jP))
            self.ui.Ri_lineEdit.setText(str(item.R_i))
            self.ui.Rj_lineEdit.setText(str(item.R_j))


        self.item = item

        self._show()
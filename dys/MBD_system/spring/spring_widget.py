'''
Created on 28. jul. 2015

@author: lskrinjar (email: skrinjar.luka@gmail.com)
'''
import datetime
from pprint import pprint

import numpy as np
from PyQt4 import QtCore, QtGui


from spring_ui import Ui_Form
from MBD_system.spring.spring import Spring
from MBD_system.array2string import array2string
from MBD_system.string2array import string2array


try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s


class SpringWidget(QtGui.QDialog):
    """
    control panel interface
    """
    
    def __init__(self, parent_node=None, parent=None):
        '''
        Constructor
        '''
        super(SpringWidget, self).__init__(parent)
        self._parent = parent


        self.ui = Ui_Form()
        self.ui.setupUi(self)
        
        self.parent_node = parent_node
        
        
        self.ui.save_pushButton.setAutoDefault(True)

        #    set validators to limit user input
        __validator_dbl = QtGui.QDoubleValidator()
        __validator_int = QtGui.QIntValidator()
        
        
        #    validators
        self.ui.springID_lineEdit.setValidator(__validator_int)
        self.ui.bodyIDi_lineEdit.setValidator(__validator_int)
        self.ui.bodyIDj_lineEdit.setValidator(__validator_int)
        self.ui.uPi_lineEdit.setValidator(__validator_dbl)
        self.ui.uPj_lineEdit.setValidator(__validator_dbl)
        self.ui.K_lineEdit.setValidator(__validator_dbl)
        self.ui.D_lineEdit.setValidator(__validator_dbl)

        
        
        #    signals
        self.ui.cancel_pushButton.clicked.connect(self._cancel)
        self.ui.save_pushButton.clicked.connect(self._save)
        self.ui.save_pushButton.setFocus()
        
        
        # pprint(vars(self.group_item))
        self.ui.springID_lineEdit.setText(str(len(self.parent_node._children)))
        
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
        for child in self.parent_node._children:
            print child._name
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
        _type = self.ui.springType_comboBox.currentText()

        try:
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

            _k = float(self.ui.K_lineEdit.text())
            _d = float(self.ui.D_lineEdit.text())
            _l0 = float(self.ui.L0_lineEdit.text())


            #   create new object
            if self.item is not None:
                self.item._name = _name
                self.item.body_id_i = _body_id_i
                self.item.body_id_j = _body_id_j

                self.item.k = _k
                self.item.c = _d
                self.item.l_0 = _l0

            #   update attributes of existing object
            else:
                _type = str(_type).lower()
                self.item = Spring(_name, _type, _body_id_i, _body_id_j, _uPi, _uPj, parent=self.parent_node)

                pos = len(self.parent_node._children)

                self.parent_node._parent.springs.append(self.item)
                self._parent.ui.treeView.model().insertRow(pos, self.item, self.parent_node)

            print "self.item =", self.item
            pprint(vars(self.item))
            self.item = None
            self.close()

        except:
            print self
            QtGui.QMessageBox.warning(self, "Warning!", "Input not correct!",QtGui.QMessageBox.Cancel, QtGui.QMessageBox.NoButton,QtGui.QMessageBox.NoButton)


    def _edit(self, item=None):
        """

        :return:
        """
        self.ui.name_lineEdit.setText(item._name)

        self.ui.springID_lineEdit.setText(str(item.spring_id))

        #   set type line edit
        _index = self.ui.springType_comboBox.findText(QtCore.QString(item.spring_type.capitalize()))
        self.ui.springType_comboBox.setCurrentIndex(_index)

        #   body id i
        self.ui.bodyIDi_lineEdit.setText(str(item.body_id_i))
        #   body id j
        self.ui.bodyIDj_lineEdit.setText(str(item.body_id_j))

        #   uPi
        self.ui.uPi_lineEdit.setText(array2string(item.u_iP))
        #   uPj
        self.ui.uPj_lineEdit.setText(array2string(item.u_jP))

        #   stiffness coef
        self.ui.K_lineEdit.setText(str(item.k))
        #   dumping coef
        self.ui.D_lineEdit.setText(str(item.c))
        #   initial length
        self.ui.L0_lineEdit.setText(str(item.l_0))

        self.item = item

        self._show()
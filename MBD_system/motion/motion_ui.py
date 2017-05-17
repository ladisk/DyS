# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'motion.ui'
#
# Created by: PyQt4 UI code generator 4.11.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName(_fromUtf8("Form"))
        Form.resize(313, 488)
        Form.setMinimumSize(QtCore.QSize(313, 488))
        Form.setMaximumSize(QtCore.QSize(313, 488))
        self.verticalLayout = QtGui.QVBoxLayout(Form)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.groupBox = QtGui.QGroupBox(Form)
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.gridLayout_3 = QtGui.QGridLayout(self.groupBox)
        self.gridLayout_3.setObjectName(_fromUtf8("gridLayout_3"))
        self.name_label = QtGui.QLabel(self.groupBox)
        self.name_label.setObjectName(_fromUtf8("name_label"))
        self.gridLayout_3.addWidget(self.name_label, 0, 0, 1, 1)
        self.name_lineEdit = QtGui.QLineEdit(self.groupBox)
        self.name_lineEdit.setObjectName(_fromUtf8("name_lineEdit"))
        self.gridLayout_3.addWidget(self.name_lineEdit, 0, 1, 1, 1)
        self.motion_id_label = QtGui.QLabel(self.groupBox)
        self.motion_id_label.setObjectName(_fromUtf8("motion_id_label"))
        self.gridLayout_3.addWidget(self.motion_id_label, 1, 0, 1, 1)
        self.motion_id_lineEdit = QtGui.QLineEdit(self.groupBox)
        self.motion_id_lineEdit.setEnabled(False)
        self.motion_id_lineEdit.setObjectName(_fromUtf8("motion_id_lineEdit"))
        self.gridLayout_3.addWidget(self.motion_id_lineEdit, 1, 1, 1, 1)
        self.verticalLayout.addWidget(self.groupBox)
        self.groupBox_simulationParameters = QtGui.QGroupBox(Form)
        self.groupBox_simulationParameters.setObjectName(_fromUtf8("groupBox_simulationParameters"))
        self.gridLayout_4 = QtGui.QGridLayout(self.groupBox_simulationParameters)
        self.gridLayout_4.setObjectName(_fromUtf8("gridLayout_4"))
        self.bodyID_label = QtGui.QLabel(self.groupBox_simulationParameters)
        self.bodyID_label.setObjectName(_fromUtf8("bodyID_label"))
        self.gridLayout_4.addWidget(self.bodyID_label, 0, 0, 1, 1)
        self.bodyID_lineEdit = QtGui.QLineEdit(self.groupBox_simulationParameters)
        self.bodyID_lineEdit.setObjectName(_fromUtf8("bodyID_lineEdit"))
        self.gridLayout_4.addWidget(self.bodyID_lineEdit, 0, 1, 1, 1)
        self.theta_lineEdit = QtGui.QLineEdit(self.groupBox_simulationParameters)
        self.theta_lineEdit.setObjectName(_fromUtf8("theta_lineEdit"))
        self.gridLayout_4.addWidget(self.theta_lineEdit, 4, 1, 1, 1)
        self.dRx_lineEdit = QtGui.QLineEdit(self.groupBox_simulationParameters)
        self.dRx_lineEdit.setObjectName(_fromUtf8("dRx_lineEdit"))
        self.gridLayout_4.addWidget(self.dRx_lineEdit, 5, 1, 1, 1)
        self.dRy_lineEdit = QtGui.QLineEdit(self.groupBox_simulationParameters)
        self.dRy_lineEdit.setObjectName(_fromUtf8("dRy_lineEdit"))
        self.gridLayout_4.addWidget(self.dRy_lineEdit, 6, 1, 1, 1)
        self.Rx_lineEdit = QtGui.QLineEdit(self.groupBox_simulationParameters)
        self.Rx_lineEdit.setObjectName(_fromUtf8("Rx_lineEdit"))
        self.gridLayout_4.addWidget(self.Rx_lineEdit, 2, 1, 1, 1)
        self.line = QtGui.QFrame(self.groupBox_simulationParameters)
        self.line.setFrameShape(QtGui.QFrame.HLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setObjectName(_fromUtf8("line"))
        self.gridLayout_4.addWidget(self.line, 1, 0, 1, 2)
        self.Ry_lineEdit = QtGui.QLineEdit(self.groupBox_simulationParameters)
        self.Ry_lineEdit.setObjectName(_fromUtf8("Ry_lineEdit"))
        self.gridLayout_4.addWidget(self.Ry_lineEdit, 3, 1, 1, 1)
        self.dtheta_lineEdit = QtGui.QLineEdit(self.groupBox_simulationParameters)
        self.dtheta_lineEdit.setObjectName(_fromUtf8("dtheta_lineEdit"))
        self.gridLayout_4.addWidget(self.dtheta_lineEdit, 7, 1, 1, 1)
        self.Rx_label = QtGui.QLabel(self.groupBox_simulationParameters)
        self.Rx_label.setObjectName(_fromUtf8("Rx_label"))
        self.gridLayout_4.addWidget(self.Rx_label, 2, 0, 1, 1)
        self.Ry_label = QtGui.QLabel(self.groupBox_simulationParameters)
        self.Ry_label.setObjectName(_fromUtf8("Ry_label"))
        self.gridLayout_4.addWidget(self.Ry_label, 3, 0, 1, 1)
        self.dRy_label = QtGui.QLabel(self.groupBox_simulationParameters)
        self.dRy_label.setObjectName(_fromUtf8("dRy_label"))
        self.gridLayout_4.addWidget(self.dRy_label, 6, 0, 1, 1)
        self.dRx_label = QtGui.QLabel(self.groupBox_simulationParameters)
        self.dRx_label.setObjectName(_fromUtf8("dRx_label"))
        self.gridLayout_4.addWidget(self.dRx_label, 5, 0, 1, 1)
        self.theta_label = QtGui.QLabel(self.groupBox_simulationParameters)
        self.theta_label.setObjectName(_fromUtf8("theta_label"))
        self.gridLayout_4.addWidget(self.theta_label, 4, 0, 1, 1)
        self.dtheta_label = QtGui.QLabel(self.groupBox_simulationParameters)
        self.dtheta_label.setObjectName(_fromUtf8("dtheta_label"))
        self.gridLayout_4.addWidget(self.dtheta_label, 7, 0, 1, 1)
        self.verticalLayout.addWidget(self.groupBox_simulationParameters)
        self.groupBox_comment = QtGui.QGroupBox(Form)
        self.groupBox_comment.setMinimumSize(QtCore.QSize(0, 40))
        self.groupBox_comment.setObjectName(_fromUtf8("groupBox_comment"))
        self.horizontalLayout = QtGui.QHBoxLayout(self.groupBox_comment)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.comment_textEdit = QtGui.QTextEdit(self.groupBox_comment)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.comment_textEdit.sizePolicy().hasHeightForWidth())
        self.comment_textEdit.setSizePolicy(sizePolicy)
        self.comment_textEdit.setMinimumSize(QtCore.QSize(260, 20))
        self.comment_textEdit.setMaximumSize(QtCore.QSize(260, 30))
        self.comment_textEdit.setObjectName(_fromUtf8("comment_textEdit"))
        self.horizontalLayout.addWidget(self.comment_textEdit)
        self.verticalLayout.addWidget(self.groupBox_comment)
        self.buttonBox = QtGui.QDialogButtonBox(Form)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Save)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.verticalLayout.addWidget(self.buttonBox)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(_translate("Form", "Motion", None))
        self.groupBox.setTitle(_translate("Form", "Basic Info", None))
        self.name_label.setText(_translate("Form", "Name", None))
        self.motion_id_label.setText(_translate("Form", "Motion id", None))
        self.groupBox_simulationParameters.setTitle(_translate("Form", "Simulation Parameters", None))
        self.bodyID_label.setText(_translate("Form", "Body ID i", None))
        self.Rx_label.setText(_translate("Form", "Rx", None))
        self.Ry_label.setText(_translate("Form", "Ry", None))
        self.dRy_label.setText(_translate("Form", "dRy", None))
        self.dRx_label.setText(_translate("Form", "dRx", None))
        self.theta_label.setText(_translate("Form", "theta", None))
        self.dtheta_label.setText(_translate("Form", "dtheta", None))
        self.groupBox_comment.setTitle(_translate("Form", "Comment", None))


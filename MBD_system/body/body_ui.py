# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'body.ui'
#
# Created: Sat Aug 15 14:11:14 2015
#      by: PyQt4 UI code generator 4.10.4
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
        Form.resize(300, 585)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Form.sizePolicy().hasHeightForWidth())
        Form.setSizePolicy(sizePolicy)
        Form.setMinimumSize(QtCore.QSize(300, 585))
        Form.setMaximumSize(QtCore.QSize(300, 585))
        self.verticalLayout = QtGui.QVBoxLayout(Form)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.groupBox_basicInfo = QtGui.QGroupBox(Form)
        self.groupBox_basicInfo.setObjectName(_fromUtf8("groupBox_basicInfo"))
        self.gridLayout = QtGui.QGridLayout(self.groupBox_basicInfo)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.bodyID_label = QtGui.QLabel(self.groupBox_basicInfo)
        self.bodyID_label.setObjectName(_fromUtf8("bodyID_label"))
        self.gridLayout.addWidget(self.bodyID_label, 1, 0, 1, 1)
        self.name_label = QtGui.QLabel(self.groupBox_basicInfo)
        self.name_label.setObjectName(_fromUtf8("name_label"))
        self.gridLayout.addWidget(self.name_label, 0, 0, 1, 1)
        self.bodyID_label_2 = QtGui.QLabel(self.groupBox_basicInfo)
        self.bodyID_label_2.setObjectName(_fromUtf8("bodyID_label_2"))
        self.gridLayout.addWidget(self.bodyID_label_2, 3, 0, 1, 1)
        self.name_lineEdit = QtGui.QLineEdit(self.groupBox_basicInfo)
        self.name_lineEdit.setObjectName(_fromUtf8("name_lineEdit"))
        self.gridLayout.addWidget(self.name_lineEdit, 0, 1, 1, 1)
        self.line = QtGui.QFrame(self.groupBox_basicInfo)
        self.line.setFrameShape(QtGui.QFrame.HLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setObjectName(_fromUtf8("line"))
        self.gridLayout.addWidget(self.line, 2, 0, 1, 2)
        self.bodyID_lineEdit = QtGui.QLineEdit(self.groupBox_basicInfo)
        self.bodyID_lineEdit.setEnabled(False)
        self.bodyID_lineEdit.setObjectName(_fromUtf8("bodyID_lineEdit"))
        self.gridLayout.addWidget(self.bodyID_lineEdit, 1, 1, 1, 1)
        self.load_data_pushButton = QtGui.QPushButton(self.groupBox_basicInfo)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.load_data_pushButton.sizePolicy().hasHeightForWidth())
        self.load_data_pushButton.setSizePolicy(sizePolicy)
        self.load_data_pushButton.setMaximumSize(QtCore.QSize(16777215, 20))
        self.load_data_pushButton.setObjectName(_fromUtf8("load_data_pushButton"))
        self.gridLayout.addWidget(self.load_data_pushButton, 4, 1, 1, 1)
        self.data_file_lineEdit = QtGui.QLineEdit(self.groupBox_basicInfo)
        self.data_file_lineEdit.setObjectName(_fromUtf8("data_file_lineEdit"))
        self.gridLayout.addWidget(self.data_file_lineEdit, 3, 1, 1, 1)
        self.verticalLayout.addWidget(self.groupBox_basicInfo)
        self.groupBox_simulationParameters = QtGui.QGroupBox(Form)
        self.groupBox_simulationParameters.setObjectName(_fromUtf8("groupBox_simulationParameters"))
        self.gridLayout_2 = QtGui.QGridLayout(self.groupBox_simulationParameters)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.J_label = QtGui.QLabel(self.groupBox_simulationParameters)
        self.J_label.setObjectName(_fromUtf8("J_label"))
        self.gridLayout_2.addWidget(self.J_label, 1, 0, 1, 1)
        self.R_label = QtGui.QLabel(self.groupBox_simulationParameters)
        self.R_label.setObjectName(_fromUtf8("R_label"))
        self.gridLayout_2.addWidget(self.R_label, 2, 0, 1, 1)
        self.R_lineEdit = QtGui.QLineEdit(self.groupBox_simulationParameters)
        self.R_lineEdit.setObjectName(_fromUtf8("R_lineEdit"))
        self.gridLayout_2.addWidget(self.R_lineEdit, 2, 1, 1, 1)
        self.dR_label = QtGui.QLabel(self.groupBox_simulationParameters)
        self.dR_label.setObjectName(_fromUtf8("dR_label"))
        self.gridLayout_2.addWidget(self.dR_label, 4, 0, 1, 1)
        self.dR_lineEdit = QtGui.QLineEdit(self.groupBox_simulationParameters)
        self.dR_lineEdit.setObjectName(_fromUtf8("dR_lineEdit"))
        self.gridLayout_2.addWidget(self.dR_lineEdit, 4, 1, 1, 1)
        self.m_lineEdit = QtGui.QLineEdit(self.groupBox_simulationParameters)
        self.m_lineEdit.setObjectName(_fromUtf8("m_lineEdit"))
        self.gridLayout_2.addWidget(self.m_lineEdit, 0, 1, 1, 1)
        self.mass_label = QtGui.QLabel(self.groupBox_simulationParameters)
        self.mass_label.setObjectName(_fromUtf8("mass_label"))
        self.gridLayout_2.addWidget(self.mass_label, 0, 0, 1, 1)
        self.J_lineEdit = QtGui.QLineEdit(self.groupBox_simulationParameters)
        self.J_lineEdit.setObjectName(_fromUtf8("J_lineEdit"))
        self.gridLayout_2.addWidget(self.J_lineEdit, 1, 1, 1, 1)
        self.theta_label = QtGui.QLabel(self.groupBox_simulationParameters)
        self.theta_label.setObjectName(_fromUtf8("theta_label"))
        self.gridLayout_2.addWidget(self.theta_label, 3, 0, 1, 1)
        self.theta_lineEdit = QtGui.QLineEdit(self.groupBox_simulationParameters)
        self.theta_lineEdit.setObjectName(_fromUtf8("theta_lineEdit"))
        self.gridLayout_2.addWidget(self.theta_lineEdit, 3, 1, 1, 1)
        self.dtheta_label = QtGui.QLabel(self.groupBox_simulationParameters)
        self.dtheta_label.setObjectName(_fromUtf8("dtheta_label"))
        self.gridLayout_2.addWidget(self.dtheta_label, 5, 0, 1, 1)
        self.dtheta_lineEdit = QtGui.QLineEdit(self.groupBox_simulationParameters)
        self.dtheta_lineEdit.setObjectName(_fromUtf8("dtheta_lineEdit"))
        self.gridLayout_2.addWidget(self.dtheta_lineEdit, 5, 1, 1, 1)
        self.verticalLayout.addWidget(self.groupBox_simulationParameters)
        self.groupBox_Geometry = QtGui.QGroupBox(Form)
        self.groupBox_Geometry.setObjectName(_fromUtf8("groupBox_Geometry"))
        self.gridLayout_3 = QtGui.QGridLayout(self.groupBox_Geometry)
        self.gridLayout_3.setObjectName(_fromUtf8("gridLayout_3"))
        self.centerOfMass_label = QtGui.QLabel(self.groupBox_Geometry)
        self.centerOfMass_label.setObjectName(_fromUtf8("centerOfMass_label"))
        self.gridLayout_3.addWidget(self.centerOfMass_label, 0, 0, 1, 1)
        self.transparent_label = QtGui.QLabel(self.groupBox_Geometry)
        self.transparent_label.setObjectName(_fromUtf8("transparent_label"))
        self.gridLayout_3.addWidget(self.transparent_label, 7, 0, 1, 1)
        self.geometrySTLfile_lineEdit = QtGui.QLineEdit(self.groupBox_Geometry)
        self.geometrySTLfile_lineEdit.setObjectName(_fromUtf8("geometrySTLfile_lineEdit"))
        self.gridLayout_3.addWidget(self.geometrySTLfile_lineEdit, 2, 1, 1, 1)
        self.display_style_comboBox = QtGui.QComboBox(self.groupBox_Geometry)
        self.display_style_comboBox.setObjectName(_fromUtf8("display_style_comboBox"))
        self.display_style_comboBox.addItem(_fromUtf8(""))
        self.display_style_comboBox.addItem(_fromUtf8(""))
        self.display_style_comboBox.addItem(_fromUtf8(""))
        self.gridLayout_3.addWidget(self.display_style_comboBox, 8, 1, 1, 1)
        self.geometrySTLfile_label = QtGui.QLabel(self.groupBox_Geometry)
        self.geometrySTLfile_label.setObjectName(_fromUtf8("geometrySTLfile_label"))
        self.gridLayout_3.addWidget(self.geometrySTLfile_label, 2, 0, 1, 1)
        self.color_lineEdit = QtGui.QLineEdit(self.groupBox_Geometry)
        self.color_lineEdit.setObjectName(_fromUtf8("color_lineEdit"))
        self.gridLayout_3.addWidget(self.color_lineEdit, 6, 1, 1, 1)
        self.display_style_label = QtGui.QLabel(self.groupBox_Geometry)
        self.display_style_label.setObjectName(_fromUtf8("display_style_label"))
        self.gridLayout_3.addWidget(self.display_style_label, 8, 0, 1, 1)
        self.color_label = QtGui.QLabel(self.groupBox_Geometry)
        self.color_label.setObjectName(_fromUtf8("color_label"))
        self.gridLayout_3.addWidget(self.color_label, 6, 0, 1, 1)
        self.centerOfMass_lineEdit = QtGui.QLineEdit(self.groupBox_Geometry)
        self.centerOfMass_lineEdit.setObjectName(_fromUtf8("centerOfMass_lineEdit"))
        self.gridLayout_3.addWidget(self.centerOfMass_lineEdit, 0, 1, 1, 1)
        self.transparent_lineEdit = QtGui.QLineEdit(self.groupBox_Geometry)
        self.transparent_lineEdit.setObjectName(_fromUtf8("transparent_lineEdit"))
        self.gridLayout_3.addWidget(self.transparent_lineEdit, 7, 1, 1, 1)
        self.load_geometry_pushButton = QtGui.QPushButton(self.groupBox_Geometry)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.load_geometry_pushButton.sizePolicy().hasHeightForWidth())
        self.load_geometry_pushButton.setSizePolicy(sizePolicy)
        self.load_geometry_pushButton.setMinimumSize(QtCore.QSize(0, 20))
        self.load_geometry_pushButton.setObjectName(_fromUtf8("load_geometry_pushButton"))
        self.gridLayout_3.addWidget(self.load_geometry_pushButton, 4, 1, 1, 1)
        self.line_2 = QtGui.QFrame(self.groupBox_Geometry)
        self.line_2.setFrameShape(QtGui.QFrame.HLine)
        self.line_2.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_2.setObjectName(_fromUtf8("line_2"))
        self.gridLayout_3.addWidget(self.line_2, 5, 0, 1, 2)
        self.verticalLayout.addWidget(self.groupBox_Geometry)
        self.bottom_horizontalLayout = QtGui.QHBoxLayout()
        self.bottom_horizontalLayout.setSpacing(6)
        self.bottom_horizontalLayout.setObjectName(_fromUtf8("bottom_horizontalLayout"))
        self.info_lineEdit = QtGui.QLineEdit(Form)
        self.info_lineEdit.setEnabled(False)
        self.info_lineEdit.setObjectName(_fromUtf8("info_lineEdit"))
        self.bottom_horizontalLayout.addWidget(self.info_lineEdit)
        self.save_pushButton = QtGui.QPushButton(Form)
        self.save_pushButton.setObjectName(_fromUtf8("save_pushButton"))
        self.bottom_horizontalLayout.addWidget(self.save_pushButton)
        self.cancel_pushButton = QtGui.QPushButton(Form)
        self.cancel_pushButton.setObjectName(_fromUtf8("cancel_pushButton"))
        self.bottom_horizontalLayout.addWidget(self.cancel_pushButton)
        self.verticalLayout.addLayout(self.bottom_horizontalLayout)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)
        Form.setTabOrder(self.name_lineEdit, self.bodyID_lineEdit)
        Form.setTabOrder(self.bodyID_lineEdit, self.m_lineEdit)
        Form.setTabOrder(self.m_lineEdit, self.J_lineEdit)
        Form.setTabOrder(self.J_lineEdit, self.R_lineEdit)
        Form.setTabOrder(self.R_lineEdit, self.theta_lineEdit)
        Form.setTabOrder(self.theta_lineEdit, self.dR_lineEdit)
        Form.setTabOrder(self.dR_lineEdit, self.dtheta_lineEdit)
        Form.setTabOrder(self.dtheta_lineEdit, self.centerOfMass_lineEdit)
        Form.setTabOrder(self.centerOfMass_lineEdit, self.geometrySTLfile_lineEdit)
        Form.setTabOrder(self.geometrySTLfile_lineEdit, self.color_lineEdit)
        Form.setTabOrder(self.color_lineEdit, self.transparent_lineEdit)
        Form.setTabOrder(self.transparent_lineEdit, self.display_style_comboBox)
        Form.setTabOrder(self.display_style_comboBox, self.info_lineEdit)
        Form.setTabOrder(self.info_lineEdit, self.save_pushButton)
        Form.setTabOrder(self.save_pushButton, self.cancel_pushButton)

    def retranslateUi(self, Form):
        Form.setWindowTitle(_translate("Form", "Body", None))
        self.groupBox_basicInfo.setTitle(_translate("Form", "Basic Info", None))
        self.bodyID_label.setText(_translate("Form", "Body id", None))
        self.name_label.setText(_translate("Form", "Name", None))
        self.bodyID_label_2.setText(_translate("Form", "Load from .dat file", None))
        self.load_data_pushButton.setText(_translate("Form", "Load", None))
        self.groupBox_simulationParameters.setTitle(_translate("Form", "Simulation Parameters", None))
        self.J_label.setText(_translate("Form", "J, kgm^2", None))
        self.R_label.setText(_translate("Form", "R, m", None))
        self.dR_label.setText(_translate("Form", "dR, m/s", None))
        self.mass_label.setText(_translate("Form", "m, kg", None))
        self.theta_label.setText(_translate("Form", "theta, deg", None))
        self.dtheta_label.setText(_translate("Form", "dtheta, deg/s", None))
        self.groupBox_Geometry.setTitle(_translate("Form", "Visualization (Geometry)", None))
        self.centerOfMass_label.setText(_translate("Form", "CM (in CAD), m", None))
        self.transparent_label.setText(_translate("Form", "Transparent", None))
        self.display_style_comboBox.setItemText(0, _translate("Form", "Filled", None))
        self.display_style_comboBox.setItemText(1, _translate("Form", "Wireframe", None))
        self.display_style_comboBox.setItemText(2, _translate("Form", "Points", None))
        self.geometrySTLfile_label.setText(_translate("Form", "Geometry STL file", None))
        self.display_style_label.setText(_translate("Form", "Display Style", None))
        self.color_label.setText(_translate("Form", "Color, RGB", None))
        self.load_geometry_pushButton.setText(_translate("Form", "Load", None))
        self.save_pushButton.setText(_translate("Form", "Save", None))
        self.cancel_pushButton.setText(_translate("Form", "Cancel", None))


# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MBD_system.ui'
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
        Form.resize(698, 425)
        Form.setMinimumSize(QtCore.QSize(698, 425))
        Form.setMaximumSize(QtCore.QSize(698, 425))
        self.gridLayout_4 = QtGui.QGridLayout(Form)
        self.gridLayout_4.setObjectName(_fromUtf8("gridLayout_4"))
        self.gridLayout_2 = QtGui.QGridLayout()
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.General_groupBox = QtGui.QGroupBox(Form)
        self.General_groupBox.setObjectName(_fromUtf8("General_groupBox"))
        self.gridLayout = QtGui.QGridLayout(self.General_groupBox)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.project_folder_label = QtGui.QLabel(self.General_groupBox)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.project_folder_label.sizePolicy().hasHeightForWidth())
        self.project_folder_label.setSizePolicy(sizePolicy)
        self.project_folder_label.setMinimumSize(QtCore.QSize(200, 0))
        self.project_folder_label.setObjectName(_fromUtf8("project_folder_label"))
        self.gridLayout.addWidget(self.project_folder_label, 1, 0, 1, 1)
        self.solutionFiletype_comboBox = QtGui.QComboBox(self.General_groupBox)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.solutionFiletype_comboBox.sizePolicy().hasHeightForWidth())
        self.solutionFiletype_comboBox.setSizePolicy(sizePolicy)
        self.solutionFiletype_comboBox.setMinimumSize(QtCore.QSize(200, 0))
        self.solutionFiletype_comboBox.setObjectName(_fromUtf8("solutionFiletype_comboBox"))
        self.solutionFiletype_comboBox.addItem(_fromUtf8(""))
        self.solutionFiletype_comboBox.addItem(_fromUtf8(""))
        self.solutionFiletype_comboBox.addItem(_fromUtf8(""))
        self.gridLayout.addWidget(self.solutionFiletype_comboBox, 3, 2, 1, 1)
        self.saveSolutionOptions_label = QtGui.QLabel(self.General_groupBox)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.saveSolutionOptions_label.sizePolicy().hasHeightForWidth())
        self.saveSolutionOptions_label.setSizePolicy(sizePolicy)
        self.saveSolutionOptions_label.setMinimumSize(QtCore.QSize(100, 0))
        self.saveSolutionOptions_label.setObjectName(_fromUtf8("saveSolutionOptions_label"))
        self.gridLayout.addWidget(self.saveSolutionOptions_label, 5, 0, 1, 1)
        self.filename_lineEdit = QtGui.QLineEdit(self.General_groupBox)
        self.filename_lineEdit.setObjectName(_fromUtf8("filename_lineEdit"))
        self.gridLayout.addWidget(self.filename_lineEdit, 0, 2, 1, 1)
        self.solutionFiletype_label = QtGui.QLabel(self.General_groupBox)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.solutionFiletype_label.sizePolicy().hasHeightForWidth())
        self.solutionFiletype_label.setSizePolicy(sizePolicy)
        self.solutionFiletype_label.setMinimumSize(QtCore.QSize(200, 0))
        self.solutionFiletype_label.setObjectName(_fromUtf8("solutionFiletype_label"))
        self.gridLayout.addWidget(self.solutionFiletype_label, 3, 0, 1, 1)
        self.filename_label = QtGui.QLabel(self.General_groupBox)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.filename_label.sizePolicy().hasHeightForWidth())
        self.filename_label.setSizePolicy(sizePolicy)
        self.filename_label.setMinimumSize(QtCore.QSize(200, 0))
        self.filename_label.setObjectName(_fromUtf8("filename_label"))
        self.gridLayout.addWidget(self.filename_label, 0, 0, 1, 1)
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem, 0, 3, 1, 1)
        spacerItem1 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem1, 3, 3, 1, 1)
        spacerItem2 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem2, 5, 3, 1, 1)
        self.saveDataToFile_comboBox = QtGui.QComboBox(self.General_groupBox)
        self.saveDataToFile_comboBox.setObjectName(_fromUtf8("saveDataToFile_comboBox"))
        self.saveDataToFile_comboBox.addItem(_fromUtf8(""))
        self.saveDataToFile_comboBox.addItem(_fromUtf8(""))
        self.saveDataToFile_comboBox.addItem(_fromUtf8(""))
        self.gridLayout.addWidget(self.saveDataToFile_comboBox, 5, 2, 1, 1)
        self.project_folder_lineEdit = QtGui.QLineEdit(self.General_groupBox)
        self.project_folder_lineEdit.setObjectName(_fromUtf8("project_folder_lineEdit"))
        self.gridLayout.addWidget(self.project_folder_lineEdit, 2, 0, 1, 4)
        self.set_folder_pushButton = QtGui.QPushButton(self.General_groupBox)
        self.set_folder_pushButton.setEnabled(True)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.set_folder_pushButton.sizePolicy().hasHeightForWidth())
        self.set_folder_pushButton.setSizePolicy(sizePolicy)
        self.set_folder_pushButton.setObjectName(_fromUtf8("set_folder_pushButton"))
        self.gridLayout.addWidget(self.set_folder_pushButton, 1, 2, 1, 1)
        self.gridLayout_2.addWidget(self.General_groupBox, 0, 0, 1, 1)
        self.Simulation_grouoBox = QtGui.QGroupBox(Form)
        self.Simulation_grouoBox.setObjectName(_fromUtf8("Simulation_grouoBox"))
        self.gridLayout_3 = QtGui.QGridLayout(self.Simulation_grouoBox)
        self.gridLayout_3.setObjectName(_fromUtf8("gridLayout_3"))
        self.gravity_vector_label = QtGui.QLabel(self.Simulation_grouoBox)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.gravity_vector_label.sizePolicy().hasHeightForWidth())
        self.gravity_vector_label.setSizePolicy(sizePolicy)
        self.gravity_vector_label.setMinimumSize(QtCore.QSize(200, 0))
        self.gravity_vector_label.setObjectName(_fromUtf8("gravity_vector_label"))
        self.gridLayout_3.addWidget(self.gravity_vector_label, 0, 0, 1, 1)
        self.gravity_value_lineEdit = QtGui.QLineEdit(self.Simulation_grouoBox)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.gravity_value_lineEdit.sizePolicy().hasHeightForWidth())
        self.gravity_value_lineEdit.setSizePolicy(sizePolicy)
        self.gravity_value_lineEdit.setMinimumSize(QtCore.QSize(100, 0))
        self.gravity_value_lineEdit.setObjectName(_fromUtf8("gravity_value_lineEdit"))
        self.gridLayout_3.addWidget(self.gravity_value_lineEdit, 1, 1, 1, 1)
        self.gravity_value_label = QtGui.QLabel(self.Simulation_grouoBox)
        self.gravity_value_label.setObjectName(_fromUtf8("gravity_value_label"))
        self.gridLayout_3.addWidget(self.gravity_value_label, 1, 0, 1, 1)
        self.gravity_vector_lineEdit = QtGui.QLineEdit(self.Simulation_grouoBox)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.gravity_vector_lineEdit.sizePolicy().hasHeightForWidth())
        self.gravity_vector_lineEdit.setSizePolicy(sizePolicy)
        self.gravity_vector_lineEdit.setMinimumSize(QtCore.QSize(100, 0))
        self.gravity_vector_lineEdit.setObjectName(_fromUtf8("gravity_vector_lineEdit"))
        self.gridLayout_3.addWidget(self.gravity_vector_lineEdit, 0, 1, 1, 1)
        spacerItem3 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout_3.addItem(spacerItem3, 0, 2, 1, 1)
        self.gridLayout_2.addWidget(self.Simulation_grouoBox, 1, 0, 1, 1)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        spacerItem4 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem4)
        self.save_pushButton = QtGui.QPushButton(Form)
        self.save_pushButton.setObjectName(_fromUtf8("save_pushButton"))
        self.horizontalLayout.addWidget(self.save_pushButton)
        self.cancel_pushButton = QtGui.QPushButton(Form)
        self.cancel_pushButton.setObjectName(_fromUtf8("cancel_pushButton"))
        self.horizontalLayout.addWidget(self.cancel_pushButton)
        self.gridLayout_2.addLayout(self.horizontalLayout, 4, 0, 1, 1)
        spacerItem5 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout_2.addItem(spacerItem5, 3, 0, 1, 1)
        self.Visualisation_groupBox = QtGui.QGroupBox(Form)
        self.Visualisation_groupBox.setObjectName(_fromUtf8("Visualisation_groupBox"))
        self.gridLayout_5 = QtGui.QGridLayout(self.Visualisation_groupBox)
        self.gridLayout_5.setObjectName(_fromUtf8("gridLayout_5"))
        self.scale_factor_lineEdit = QtGui.QLineEdit(self.Visualisation_groupBox)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.scale_factor_lineEdit.sizePolicy().hasHeightForWidth())
        self.scale_factor_lineEdit.setSizePolicy(sizePolicy)
        self.scale_factor_lineEdit.setMinimumSize(QtCore.QSize(100, 0))
        self.scale_factor_lineEdit.setObjectName(_fromUtf8("scale_factor_lineEdit"))
        self.gridLayout_5.addWidget(self.scale_factor_lineEdit, 0, 1, 1, 1)
        spacerItem6 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout_5.addItem(spacerItem6, 0, 4, 1, 1)
        self.scale_factor_label = QtGui.QLabel(self.Visualisation_groupBox)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.scale_factor_label.sizePolicy().hasHeightForWidth())
        self.scale_factor_label.setSizePolicy(sizePolicy)
        self.scale_factor_label.setMinimumSize(QtCore.QSize(200, 0))
        self.scale_factor_label.setObjectName(_fromUtf8("scale_factor_label"))
        self.gridLayout_5.addWidget(self.scale_factor_label, 0, 0, 1, 1)
        self.line = QtGui.QFrame(self.Visualisation_groupBox)
        self.line.setFrameShape(QtGui.QFrame.VLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setObjectName(_fromUtf8("line"))
        self.gridLayout_5.addWidget(self.line, 0, 2, 1, 1)
        self.showGCS_checkBox = QtGui.QCheckBox(self.Visualisation_groupBox)
        self.showGCS_checkBox.setTristate(False)
        self.showGCS_checkBox.setObjectName(_fromUtf8("showGCS_checkBox"))
        self.gridLayout_5.addWidget(self.showGCS_checkBox, 0, 3, 1, 1)
        self.gridLayout_2.addWidget(self.Visualisation_groupBox, 2, 0, 1, 1)
        self.gridLayout_4.addLayout(self.gridLayout_2, 0, 0, 1, 1)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(_translate("Form", "Properties", None))
        self.General_groupBox.setTitle(_translate("Form", "General", None))
        self.project_folder_label.setText(_translate("Form", "Project Folder (abs. path)", None))
        self.solutionFiletype_comboBox.setItemText(0, _translate("Form", ".sol", None))
        self.solutionFiletype_comboBox.setItemText(1, _translate("Form", ".dat", None))
        self.solutionFiletype_comboBox.setItemText(2, _translate("Form", ".xlsx", None))
        self.saveSolutionOptions_label.setText(_translate("Form", "Save Solution options", None))
        self.solutionFiletype_label.setText(_translate("Form", "Solution filetype", None))
        self.filename_label.setText(_translate("Form", "Filename", None))
        self.saveDataToFile_comboBox.setItemText(0, _translate("Form", "Discard", None))
        self.saveDataToFile_comboBox.setItemText(1, _translate("Form", "Overwrite (existing file)", None))
        self.saveDataToFile_comboBox.setItemText(2, _translate("Form", "Save to new (next available) file", None))
        self.set_folder_pushButton.setText(_translate("Form", "Set Folder", None))
        self.Simulation_grouoBox.setTitle(_translate("Form", "Simulation", None))
        self.gravity_vector_label.setText(_translate("Form", "Gravity - unit vector, [1, 1, 1]", None))
        self.gravity_value_label.setText(_translate("Form", "Gravity - magnitude, [m/s^2]", None))
        self.save_pushButton.setText(_translate("Form", "Save", None))
        self.cancel_pushButton.setText(_translate("Form", "Cancel", None))
        self.Visualisation_groupBox.setTitle(_translate("Form", "Visualization", None))
        self.scale_factor_label.setText(_translate("Form", "Scale Factor", None))
        self.showGCS_checkBox.setText(_translate("Form", "Show GCS", None))


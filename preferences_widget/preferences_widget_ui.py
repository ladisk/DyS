# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'preferences_widget.ui'
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
        Form.resize(769, 499)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Form.sizePolicy().hasHeightForWidth())
        Form.setSizePolicy(sizePolicy)
        Form.setMinimumSize(QtCore.QSize(769, 499))
        self.gridLayout = QtGui.QGridLayout(Form)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.buttonBox = QtGui.QDialogButtonBox(Form)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Save)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.gridLayout.addWidget(self.buttonBox, 1, 0, 1, 1)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.listWidget = QtGui.QListWidget(Form)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.listWidget.sizePolicy().hasHeightForWidth())
        self.listWidget.setSizePolicy(sizePolicy)
        self.listWidget.setObjectName(_fromUtf8("listWidget"))
        item = QtGui.QListWidgetItem()
        self.listWidget.addItem(item)
        item = QtGui.QListWidgetItem()
        self.listWidget.addItem(item)
        item = QtGui.QListWidgetItem()
        self.listWidget.addItem(item)
        item = QtGui.QListWidgetItem()
        self.listWidget.addItem(item)
        item = QtGui.QListWidgetItem()
        self.listWidget.addItem(item)
        self.horizontalLayout.addWidget(self.listWidget)
        self.stackedWidget = QtGui.QStackedWidget(Form)
        self.stackedWidget.setObjectName(_fromUtf8("stackedWidget"))
        self.pageMain = QtGui.QWidget()
        self.pageMain.setObjectName(_fromUtf8("pageMain"))
        self.gridLayout_3 = QtGui.QGridLayout(self.pageMain)
        self.gridLayout_3.setObjectName(_fromUtf8("gridLayout_3"))
        self.mainSettings_groupBox = QtGui.QGroupBox(self.pageMain)
        self.mainSettings_groupBox.setObjectName(_fromUtf8("mainSettings_groupBox"))
        self.gridLayout_9 = QtGui.QGridLayout(self.mainSettings_groupBox)
        self.gridLayout_9.setObjectName(_fromUtf8("gridLayout_9"))
        self.integrationMethod_label = QtGui.QLabel(self.mainSettings_groupBox)
        self.integrationMethod_label.setObjectName(_fromUtf8("integrationMethod_label"))
        self.gridLayout_9.addWidget(self.integrationMethod_label, 1, 0, 1, 1)
        self.endTime_label = QtGui.QLabel(self.mainSettings_groupBox)
        self.endTime_label.setObjectName(_fromUtf8("endTime_label"))
        self.gridLayout_9.addWidget(self.endTime_label, 2, 0, 1, 1)
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout_9.addItem(spacerItem, 0, 1, 1, 1)
        self.endTime = QtGui.QLineEdit(self.mainSettings_groupBox)
        self.endTime.setInputMethodHints(QtCore.Qt.ImhDigitsOnly|QtCore.Qt.ImhPreferNumbers)
        self.endTime.setObjectName(_fromUtf8("endTime"))
        self.gridLayout_9.addWidget(self.endTime, 2, 2, 1, 1)
        self.analysisTypeComboBox = QtGui.QComboBox(self.mainSettings_groupBox)
        self.analysisTypeComboBox.setObjectName(_fromUtf8("analysisTypeComboBox"))
        self.analysisTypeComboBox.addItem(_fromUtf8(""))
        self.analysisTypeComboBox.addItem(_fromUtf8(""))
        self.gridLayout_9.addWidget(self.analysisTypeComboBox, 0, 2, 1, 1)
        self.analysisType_label = QtGui.QLabel(self.mainSettings_groupBox)
        self.analysisType_label.setObjectName(_fromUtf8("analysisType_label"))
        self.gridLayout_9.addWidget(self.analysisType_label, 0, 0, 1, 1)
        self.integrationMethodComboBox = QtGui.QComboBox(self.mainSettings_groupBox)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.integrationMethodComboBox.sizePolicy().hasHeightForWidth())
        self.integrationMethodComboBox.setSizePolicy(sizePolicy)
        self.integrationMethodComboBox.setMinimumSize(QtCore.QSize(0, 0))
        self.integrationMethodComboBox.setObjectName(_fromUtf8("integrationMethodComboBox"))
        self.integrationMethodComboBox.addItem(_fromUtf8(""))
        self.integrationMethodComboBox.addItem(_fromUtf8(""))
        self.gridLayout_9.addWidget(self.integrationMethodComboBox, 1, 2, 1, 1)
        self.gridLayout_3.addWidget(self.mainSettings_groupBox, 1, 0, 1, 1)
        spacerItem1 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout_3.addItem(spacerItem1, 2, 0, 1, 1)
        self.stackedWidget.addWidget(self.pageMain)
        self.pageSimulation = QtGui.QWidget()
        self.pageSimulation.setObjectName(_fromUtf8("pageSimulation"))
        self.gridLayout_10 = QtGui.QGridLayout(self.pageSimulation)
        self.gridLayout_10.setObjectName(_fromUtf8("gridLayout_10"))
        self.simulationSettings_groupBox = QtGui.QGroupBox(self.pageSimulation)
        self.simulationSettings_groupBox.setObjectName(_fromUtf8("simulationSettings_groupBox"))
        self.gridLayout_11 = QtGui.QGridLayout(self.simulationSettings_groupBox)
        self.gridLayout_11.setObjectName(_fromUtf8("gridLayout_11"))
        self.relTol = QtGui.QLineEdit(self.simulationSettings_groupBox)
        self.relTol.setEnabled(True)
        self.relTol.setInputMethodHints(QtCore.Qt.ImhDigitsOnly|QtCore.Qt.ImhPreferNumbers)
        self.relTol.setObjectName(_fromUtf8("relTol"))
        self.gridLayout_11.addWidget(self.relTol, 3, 2, 1, 1)
        self.TOL_dq_i = QtGui.QLineEdit(self.simulationSettings_groupBox)
        self.TOL_dq_i.setEnabled(True)
        self.TOL_dq_i.setInputMethodHints(QtCore.Qt.ImhDigitsOnly|QtCore.Qt.ImhPreferNumbers)
        self.TOL_dq_i.setObjectName(_fromUtf8("TOL_dq_i"))
        self.gridLayout_11.addWidget(self.TOL_dq_i, 4, 2, 1, 1)
        self.TOL_C_label = QtGui.QLabel(self.simulationSettings_groupBox)
        self.TOL_C_label.setObjectName(_fromUtf8("TOL_C_label"))
        self.gridLayout_11.addWidget(self.TOL_C_label, 5, 0, 1, 1)
        self.Hmin = QtGui.QLineEdit(self.simulationSettings_groupBox)
        self.Hmin.setEnabled(True)
        self.Hmin.setInputMethodHints(QtCore.Qt.ImhDigitsOnly|QtCore.Qt.ImhPreferNumbers)
        self.Hmin.setObjectName(_fromUtf8("Hmin"))
        self.gridLayout_11.addWidget(self.Hmin, 1, 2, 1, 1)
        self.absTol_label = QtGui.QLabel(self.simulationSettings_groupBox)
        self.absTol_label.setObjectName(_fromUtf8("absTol_label"))
        self.gridLayout_11.addWidget(self.absTol_label, 2, 0, 1, 1)
        self.relTolLabel = QtGui.QLabel(self.simulationSettings_groupBox)
        self.relTolLabel.setObjectName(_fromUtf8("relTolLabel"))
        self.gridLayout_11.addWidget(self.relTolLabel, 3, 0, 1, 1)
        self.absTol = QtGui.QLineEdit(self.simulationSettings_groupBox)
        self.absTol.setEnabled(True)
        self.absTol.setInputMethodHints(QtCore.Qt.ImhDigitsOnly|QtCore.Qt.ImhPreferNumbers)
        self.absTol.setObjectName(_fromUtf8("absTol"))
        self.gridLayout_11.addWidget(self.absTol, 2, 2, 1, 1)
        self.TOL_dq_i_label = QtGui.QLabel(self.simulationSettings_groupBox)
        self.TOL_dq_i_label.setObjectName(_fromUtf8("TOL_dq_i_label"))
        self.gridLayout_11.addWidget(self.TOL_dq_i_label, 4, 0, 1, 1)
        self.Hmin_label = QtGui.QLabel(self.simulationSettings_groupBox)
        self.Hmin_label.setObjectName(_fromUtf8("Hmin_label"))
        self.gridLayout_11.addWidget(self.Hmin_label, 1, 0, 1, 1)
        self.TOL_C = QtGui.QLineEdit(self.simulationSettings_groupBox)
        self.TOL_C.setEnabled(True)
        self.TOL_C.setInputMethodHints(QtCore.Qt.ImhDigitsOnly|QtCore.Qt.ImhPreferNumbers)
        self.TOL_C.setObjectName(_fromUtf8("TOL_C"))
        self.gridLayout_11.addWidget(self.TOL_C, 5, 2, 1, 1)
        spacerItem2 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout_11.addItem(spacerItem2, 0, 1, 1, 1)
        self.Hmax_label = QtGui.QLabel(self.simulationSettings_groupBox)
        self.Hmax_label.setObjectName(_fromUtf8("Hmax_label"))
        self.gridLayout_11.addWidget(self.Hmax_label, 0, 0, 1, 1)
        self.Hmax = QtGui.QLineEdit(self.simulationSettings_groupBox)
        self.Hmax.setEnabled(True)
        self.Hmax.setMinimumSize(QtCore.QSize(0, 0))
        self.Hmax.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.Hmax.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.Hmax.setInputMethodHints(QtCore.Qt.ImhDigitsOnly|QtCore.Qt.ImhPreferNumbers)
        self.Hmax.setObjectName(_fromUtf8("Hmax"))
        self.gridLayout_11.addWidget(self.Hmax, 0, 2, 1, 1)
        self.gridLayout_10.addWidget(self.simulationSettings_groupBox, 1, 0, 1, 1)
        self.simulation_label = QtGui.QLabel(self.pageSimulation)
        self.simulation_label.setObjectName(_fromUtf8("simulation_label"))
        self.gridLayout_10.addWidget(self.simulation_label, 0, 0, 1, 1)
        spacerItem3 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout_10.addItem(spacerItem3, 2, 0, 1, 1)
        self.stackedWidget.addWidget(self.pageSimulation)
        self.pageAnimation = QtGui.QWidget()
        self.pageAnimation.setObjectName(_fromUtf8("pageAnimation"))
        self.gridLayout_6 = QtGui.QGridLayout(self.pageAnimation)
        self.gridLayout_6.setObjectName(_fromUtf8("gridLayout_6"))
        spacerItem4 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout_6.addItem(spacerItem4, 3, 0, 1, 1)
        self.animationSettings_groupBox = QtGui.QGroupBox(self.pageAnimation)
        self.animationSettings_groupBox.setObjectName(_fromUtf8("animationSettings_groupBox"))
        self.gridLayout_8 = QtGui.QGridLayout(self.animationSettings_groupBox)
        self.gridLayout_8.setObjectName(_fromUtf8("gridLayout_8"))
        self.dt_lineEdit = QtGui.QLineEdit(self.animationSettings_groupBox)
        self.dt_lineEdit.setObjectName(_fromUtf8("dt_lineEdit"))
        self.gridLayout_8.addWidget(self.dt_lineEdit, 3, 2, 1, 1)
        self.updateDisplay_label = QtGui.QLabel(self.animationSettings_groupBox)
        self.updateDisplay_label.setObjectName(_fromUtf8("updateDisplay_label"))
        self.gridLayout_8.addWidget(self.updateDisplay_label, 0, 0, 1, 1)
        self.playbackSpeed_label = QtGui.QLabel(self.animationSettings_groupBox)
        self.playbackSpeed_label.setObjectName(_fromUtf8("playbackSpeed_label"))
        self.gridLayout_8.addWidget(self.playbackSpeed_label, 4, 0, 1, 1)
        self.steps_lineEdit = QtGui.QLineEdit(self.animationSettings_groupBox)
        self.steps_lineEdit.setObjectName(_fromUtf8("steps_lineEdit"))
        self.gridLayout_8.addWidget(self.steps_lineEdit, 2, 2, 1, 1)
        self.steps_label = QtGui.QLabel(self.animationSettings_groupBox)
        self.steps_label.setObjectName(_fromUtf8("steps_label"))
        self.gridLayout_8.addWidget(self.steps_label, 3, 0, 1, 1)
        self.updateDisplay_comboBox = QtGui.QComboBox(self.animationSettings_groupBox)
        self.updateDisplay_comboBox.setObjectName(_fromUtf8("updateDisplay_comboBox"))
        self.updateDisplay_comboBox.addItem(_fromUtf8(""))
        self.updateDisplay_comboBox.addItem(_fromUtf8(""))
        self.gridLayout_8.addWidget(self.updateDisplay_comboBox, 0, 2, 1, 1)
        self.dt_label = QtGui.QLabel(self.animationSettings_groupBox)
        self.dt_label.setObjectName(_fromUtf8("dt_label"))
        self.gridLayout_8.addWidget(self.dt_label, 2, 0, 1, 1)
        self.playbackSpeed_doubleSpinBox = QtGui.QDoubleSpinBox(self.animationSettings_groupBox)
        self.playbackSpeed_doubleSpinBox.setButtonSymbols(QtGui.QAbstractSpinBox.NoButtons)
        self.playbackSpeed_doubleSpinBox.setKeyboardTracking(True)
        self.playbackSpeed_doubleSpinBox.setProperty("value", 1.0)
        self.playbackSpeed_doubleSpinBox.setObjectName(_fromUtf8("playbackSpeed_doubleSpinBox"))
        self.gridLayout_8.addWidget(self.playbackSpeed_doubleSpinBox, 4, 2, 1, 1)
        self.fps_label = QtGui.QLabel(self.animationSettings_groupBox)
        self.fps_label.setObjectName(_fromUtf8("fps_label"))
        self.gridLayout_8.addWidget(self.fps_label, 5, 0, 1, 1)
        self.fps_spinBox = QtGui.QSpinBox(self.animationSettings_groupBox)
        self.fps_spinBox.setButtonSymbols(QtGui.QAbstractSpinBox.NoButtons)
        self.fps_spinBox.setProperty("value", 24)
        self.fps_spinBox.setObjectName(_fromUtf8("fps_spinBox"))
        self.gridLayout_8.addWidget(self.fps_spinBox, 5, 2, 1, 1)
        self.gridLayout_6.addWidget(self.animationSettings_groupBox, 1, 0, 1, 1)
        self.label_5 = QtGui.QLabel(self.pageAnimation)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.gridLayout_6.addWidget(self.label_5, 0, 0, 1, 1)
        spacerItem5 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout_6.addItem(spacerItem5, 1, 1, 1, 1)
        self.stackedWidget.addWidget(self.pageAnimation)
        self.pageVisualization = QtGui.QWidget()
        self.pageVisualization.setObjectName(_fromUtf8("pageVisualization"))
        self.gridLayout_7 = QtGui.QGridLayout(self.pageVisualization)
        self.gridLayout_7.setObjectName(_fromUtf8("gridLayout_7"))
        self.visualization_label = QtGui.QLabel(self.pageVisualization)
        self.visualization_label.setObjectName(_fromUtf8("visualization_label"))
        self.gridLayout_7.addWidget(self.visualization_label, 0, 0, 1, 1)
        spacerItem6 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout_7.addItem(spacerItem6, 2, 2, 1, 1)
        spacerItem7 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout_7.addItem(spacerItem7, 4, 0, 1, 1)
        self.visualization_groupBox = QtGui.QGroupBox(self.pageVisualization)
        self.visualization_groupBox.setObjectName(_fromUtf8("visualization_groupBox"))
        self.gridLayout_5 = QtGui.QGridLayout(self.visualization_groupBox)
        self.gridLayout_5.setObjectName(_fromUtf8("gridLayout_5"))
        self.editBackgroundColor_pushButton = QtGui.QPushButton(self.visualization_groupBox)
        self.editBackgroundColor_pushButton.setObjectName(_fromUtf8("editBackgroundColor_pushButton"))
        self.gridLayout_5.addWidget(self.editBackgroundColor_pushButton, 0, 2, 1, 1)
        self.backgroundColor_lineEdit = QtGui.QLineEdit(self.visualization_groupBox)
        self.backgroundColor_lineEdit.setInputMethodHints(QtCore.Qt.ImhDigitsOnly|QtCore.Qt.ImhPreferNumbers)
        self.backgroundColor_lineEdit.setObjectName(_fromUtf8("backgroundColor_lineEdit"))
        self.gridLayout_5.addWidget(self.backgroundColor_lineEdit, 0, 1, 1, 1)
        self.updateDisplayStep = QtGui.QLineEdit(self.visualization_groupBox)
        self.updateDisplayStep.setInputMethodHints(QtCore.Qt.ImhDigitsOnly|QtCore.Qt.ImhPreferNumbers)
        self.updateDisplayStep.setObjectName(_fromUtf8("updateDisplayStep"))
        self.gridLayout_5.addWidget(self.updateDisplayStep, 2, 1, 1, 1)
        self.label_3 = QtGui.QLabel(self.visualization_groupBox)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.gridLayout_5.addWidget(self.label_3, 0, 0, 1, 1)
        self.updateDisplayStepLabel = QtGui.QLabel(self.visualization_groupBox)
        self.updateDisplayStepLabel.setObjectName(_fromUtf8("updateDisplayStepLabel"))
        self.gridLayout_5.addWidget(self.updateDisplayStepLabel, 2, 0, 1, 1)
        spacerItem8 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout_5.addItem(spacerItem8, 2, 3, 1, 1)
        self.gridLayout_7.addWidget(self.visualization_groupBox, 1, 0, 1, 3)
        self.infoDisplay_groupBox = QtGui.QGroupBox(self.pageVisualization)
        self.infoDisplay_groupBox.setObjectName(_fromUtf8("infoDisplay_groupBox"))
        self.verticalLayout = QtGui.QVBoxLayout(self.infoDisplay_groupBox)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.filename_checkBox = QtGui.QCheckBox(self.infoDisplay_groupBox)
        self.filename_checkBox.setObjectName(_fromUtf8("filename_checkBox"))
        self.verticalLayout.addWidget(self.filename_checkBox)
        self.simulationTime_checkBox = QtGui.QCheckBox(self.infoDisplay_groupBox)
        self.simulationTime_checkBox.setObjectName(_fromUtf8("simulationTime_checkBox"))
        self.verticalLayout.addWidget(self.simulationTime_checkBox)
        self.simulationStepNumber_checkBox = QtGui.QCheckBox(self.infoDisplay_groupBox)
        self.simulationStepNumber_checkBox.setObjectName(_fromUtf8("simulationStepNumber_checkBox"))
        self.verticalLayout.addWidget(self.simulationStepNumber_checkBox)
        self.timeAndDate_checkBox = QtGui.QCheckBox(self.infoDisplay_groupBox)
        self.timeAndDate_checkBox.setObjectName(_fromUtf8("timeAndDate_checkBox"))
        self.verticalLayout.addWidget(self.timeAndDate_checkBox)
        self.gridLayout_7.addWidget(self.infoDisplay_groupBox, 2, 0, 1, 1)
        self.GCS_groupBox = QtGui.QGroupBox(self.pageVisualization)
        self.GCS_groupBox.setObjectName(_fromUtf8("GCS_groupBox"))
        self.gridLayout_2 = QtGui.QGridLayout(self.GCS_groupBox)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.showGCS_checkBox = QtGui.QCheckBox(self.GCS_groupBox)
        self.showGCS_checkBox.setObjectName(_fromUtf8("showGCS_checkBox"))
        self.gridLayout_2.addWidget(self.showGCS_checkBox, 0, 0, 1, 1)
        self.showGCS_checkBox_2 = QtGui.QCheckBox(self.GCS_groupBox)
        self.showGCS_checkBox_2.setObjectName(_fromUtf8("showGCS_checkBox_2"))
        self.gridLayout_2.addWidget(self.showGCS_checkBox_2, 1, 0, 1, 1)
        self.gridLayout_7.addWidget(self.GCS_groupBox, 3, 0, 1, 1)
        self.stackedWidget.addWidget(self.pageVisualization)
        self.solutionSettings_page = QtGui.QWidget()
        self.solutionSettings_page.setObjectName(_fromUtf8("solutionSettings_page"))
        self.gridLayout_12 = QtGui.QGridLayout(self.solutionSettings_page)
        self.gridLayout_12.setObjectName(_fromUtf8("gridLayout_12"))
        self.solution_label = QtGui.QLabel(self.solutionSettings_page)
        self.solution_label.setObjectName(_fromUtf8("solution_label"))
        self.gridLayout_12.addWidget(self.solution_label, 0, 0, 1, 1)
        self.whenSimulationIsFinished_groupBox = QtGui.QGroupBox(self.solutionSettings_page)
        self.whenSimulationIsFinished_groupBox.setObjectName(_fromUtf8("whenSimulationIsFinished_groupBox"))
        self.gridLayout_4 = QtGui.QGridLayout(self.whenSimulationIsFinished_groupBox)
        self.gridLayout_4.setObjectName(_fromUtf8("gridLayout_4"))
        self.restoreInitialConditionsStatus = QtGui.QCheckBox(self.whenSimulationIsFinished_groupBox)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.restoreInitialConditionsStatus.sizePolicy().hasHeightForWidth())
        self.restoreInitialConditionsStatus.setSizePolicy(sizePolicy)
        self.restoreInitialConditionsStatus.setMinimumSize(QtCore.QSize(200, 0))
        self.restoreInitialConditionsStatus.setObjectName(_fromUtf8("restoreInitialConditionsStatus"))
        self.gridLayout_4.addWidget(self.restoreInitialConditionsStatus, 0, 0, 1, 1)
        self.loadSolutionFileStatus = QtGui.QCheckBox(self.whenSimulationIsFinished_groupBox)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.loadSolutionFileStatus.sizePolicy().hasHeightForWidth())
        self.loadSolutionFileStatus.setSizePolicy(sizePolicy)
        self.loadSolutionFileStatus.setMinimumSize(QtCore.QSize(200, 0))
        self.loadSolutionFileStatus.setObjectName(_fromUtf8("loadSolutionFileStatus"))
        self.gridLayout_4.addWidget(self.loadSolutionFileStatus, 1, 0, 1, 1)
        self.solutionFiletypes_label = QtGui.QLabel(self.whenSimulationIsFinished_groupBox)
        self.solutionFiletypes_label.setObjectName(_fromUtf8("solutionFiletypes_label"))
        self.gridLayout_4.addWidget(self.solutionFiletypes_label, 2, 0, 1, 1)
        self.solutuionFiletypes_comboBox = QtGui.QComboBox(self.whenSimulationIsFinished_groupBox)
        self.solutuionFiletypes_comboBox.setObjectName(_fromUtf8("solutuionFiletypes_comboBox"))
        self.solutuionFiletypes_comboBox.addItem(_fromUtf8(""))
        self.solutuionFiletypes_comboBox.addItem(_fromUtf8(""))
        self.solutuionFiletypes_comboBox.addItem(_fromUtf8(""))
        self.gridLayout_4.addWidget(self.solutuionFiletypes_comboBox, 2, 1, 1, 1)
        self.gridLayout_12.addWidget(self.whenSimulationIsFinished_groupBox, 1, 0, 1, 1)
        spacerItem9 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout_12.addItem(spacerItem9, 2, 0, 1, 1)
        spacerItem10 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout_12.addItem(spacerItem10, 1, 1, 1, 1)
        self.stackedWidget.addWidget(self.solutionSettings_page)
        self.horizontalLayout.addWidget(self.stackedWidget)
        self.gridLayout.addLayout(self.horizontalLayout, 0, 0, 1, 1)

        self.retranslateUi(Form)
        self.stackedWidget.setCurrentIndex(2)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(_translate("Form", "Preferences", None))
        __sortingEnabled = self.listWidget.isSortingEnabled()
        self.listWidget.setSortingEnabled(False)
        item = self.listWidget.item(0)
        item.setText(_translate("Form", "Main", None))
        item = self.listWidget.item(1)
        item.setText(_translate("Form", "Simulation", None))
        item = self.listWidget.item(2)
        item.setText(_translate("Form", "Animation", None))
        item = self.listWidget.item(3)
        item.setText(_translate("Form", "Visualization", None))
        item = self.listWidget.item(4)
        item.setText(_translate("Form", "Solution", None))
        self.listWidget.setSortingEnabled(__sortingEnabled)
        self.mainSettings_groupBox.setTitle(_translate("Form", "Main settings", None))
        self.integrationMethod_label.setText(_translate("Form", "Integration method", None))
        self.endTime_label.setText(_translate("Form", "End time, s", None))
        self.analysisTypeComboBox.setItemText(0, _translate("Form", "Kinematic", None))
        self.analysisTypeComboBox.setItemText(1, _translate("Form", "Dynamic", None))
        self.analysisType_label.setText(_translate("Form", "Analysis type", None))
        self.integrationMethodComboBox.setItemText(0, _translate("Form", "Euler", None))
        self.integrationMethodComboBox.setItemText(1, _translate("Form", "Runge-Kutta", None))
        self.simulationSettings_groupBox.setTitle(_translate("Form", "Simulation settings", None))
        self.TOL_C_label.setText(_translate("Form", "TOL C(q, t)", None))
        self.absTol_label.setText(_translate("Form", "AbsTol", None))
        self.relTolLabel.setText(_translate("Form", "RelTol", None))
        self.TOL_dq_i_label.setText(_translate("Form", "TOL dq_i", None))
        self.Hmin_label.setText(_translate("Form", "Hmin", None))
        self.Hmax_label.setText(_translate("Form", "Hmax", None))
        self.simulation_label.setText(_translate("Form", "Simulation", None))
        self.animationSettings_groupBox.setTitle(_translate("Form", "Animation settings", None))
        self.updateDisplay_label.setText(_translate("Form", "During simulation update display on every", None))
        self.playbackSpeed_label.setText(_translate("Form", "Playback speed", None))
        self.steps_label.setText(_translate("Form", "step", None))
        self.updateDisplay_comboBox.setItemText(0, _translate("Form", "dt", None))
        self.updateDisplay_comboBox.setItemText(1, _translate("Form", "step", None))
        self.dt_label.setText(_translate("Form", "dt", None))
        self.fps_label.setText(_translate("Form", "Frames per second (for video file)", None))
        self.label_5.setText(_translate("Form", "Animation", None))
        self.visualization_label.setText(_translate("Form", "Visualization", None))
        self.visualization_groupBox.setTitle(_translate("Form", "Open GL", None))
        self.editBackgroundColor_pushButton.setText(_translate("Form", "Edit", None))
        self.label_3.setText(_translate("Form", "Background color", None))
        self.updateDisplayStepLabel.setText(_translate("Form", "Update display on i-th stp", None))
        self.infoDisplay_groupBox.setTitle(_translate("Form", "Info display", None))
        self.filename_checkBox.setText(_translate("Form", "Filename", None))
        self.simulationTime_checkBox.setText(_translate("Form", "Simulation time", None))
        self.simulationStepNumber_checkBox.setText(_translate("Form", "Step number", None))
        self.timeAndDate_checkBox.setText(_translate("Form", "Time and date", None))
        self.GCS_groupBox.setTitle(_translate("Form", "GCS ", None))
        self.showGCS_checkBox.setText(_translate("Form", "Show GCS", None))
        self.showGCS_checkBox_2.setText(_translate("Form", "Show GCS View", None))
        self.solution_label.setText(_translate("Form", "Solution settings", None))
        self.whenSimulationIsFinished_groupBox.setTitle(_translate("Form", "When simulation is finished", None))
        self.restoreInitialConditionsStatus.setText(_translate("Form", "Restore initial conditions", None))
        self.loadSolutionFileStatus.setText(_translate("Form", "Load solution", None))
        self.solutionFiletypes_label.setText(_translate("Form", "Solution filetype", None))
        self.solutuionFiletypes_comboBox.setItemText(0, _translate("Form", ".sol", None))
        self.solutuionFiletypes_comboBox.setItemText(1, _translate("Form", ".csv", None))
        self.solutuionFiletypes_comboBox.setItemText(2, _translate("Form", ".xlsx", None))


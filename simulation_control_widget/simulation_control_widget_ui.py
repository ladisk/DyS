# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'simulation_control_widget.ui'
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
        Form.resize(300, 685)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Form.sizePolicy().hasHeightForWidth())
        Form.setSizePolicy(sizePolicy)
        Form.setMinimumSize(QtCore.QSize(300, 685))
        Form.setMaximumSize(QtCore.QSize(360, 685))
        self.gridLayout_2 = QtGui.QGridLayout(Form)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.animationSettings = QtGui.QGroupBox(Form)
        self.animationSettings.setObjectName(_fromUtf8("animationSettings"))
        self.gridLayout_3 = QtGui.QGridLayout(self.animationSettings)
        self.gridLayout_3.setObjectName(_fromUtf8("gridLayout_3"))
        self.updateDisplayStepLabel = QtGui.QLabel(self.animationSettings)
        self.updateDisplayStepLabel.setObjectName(_fromUtf8("updateDisplayStepLabel"))
        self.gridLayout_3.addWidget(self.updateDisplayStepLabel, 0, 0, 1, 1)
        self.numberOfSteps_label_2 = QtGui.QLabel(self.animationSettings)
        self.numberOfSteps_label_2.setObjectName(_fromUtf8("numberOfSteps_label_2"))
        self.gridLayout_3.addWidget(self.numberOfSteps_label_2, 3, 0, 1, 1)
        self.numberOfSteps = QtGui.QLineEdit(self.animationSettings)
        self.numberOfSteps.setEnabled(False)
        self.numberOfSteps.setInputMethodHints(QtCore.Qt.ImhDigitsOnly|QtCore.Qt.ImhPreferNumbers)
        self.numberOfSteps.setObjectName(_fromUtf8("numberOfSteps"))
        self.gridLayout_3.addWidget(self.numberOfSteps, 2, 1, 1, 1)
        self.numberOfSteps_label = QtGui.QLabel(self.animationSettings)
        self.numberOfSteps_label.setObjectName(_fromUtf8("numberOfSteps_label"))
        self.gridLayout_3.addWidget(self.numberOfSteps_label, 2, 0, 1, 1)
        self.currentStepLabel = QtGui.QLabel(self.animationSettings)
        self.currentStepLabel.setObjectName(_fromUtf8("currentStepLabel"))
        self.gridLayout_3.addWidget(self.currentStepLabel, 1, 0, 1, 1)
        self.updateDisplayStep = QtGui.QLineEdit(self.animationSettings)
        self.updateDisplayStep.setInputMethodHints(QtCore.Qt.ImhDigitsOnly|QtCore.Qt.ImhPreferNumbers)
        self.updateDisplayStep.setObjectName(_fromUtf8("updateDisplayStep"))
        self.gridLayout_3.addWidget(self.updateDisplayStep, 0, 1, 1, 1)
        self.playBackSpeed_doubleSpinBox = QtGui.QDoubleSpinBox(self.animationSettings)
        self.playBackSpeed_doubleSpinBox.setProperty("value", 1.0)
        self.playBackSpeed_doubleSpinBox.setObjectName(_fromUtf8("playBackSpeed_doubleSpinBox"))
        self.gridLayout_3.addWidget(self.playBackSpeed_doubleSpinBox, 3, 1, 1, 1)
        self.currentStep = QtGui.QLineEdit(self.animationSettings)
        self.currentStep.setEnabled(True)
        self.currentStep.setInputMethodHints(QtCore.Qt.ImhDigitsOnly|QtCore.Qt.ImhPreferNumbers)
        self.currentStep.setObjectName(_fromUtf8("currentStep"))
        self.gridLayout_3.addWidget(self.currentStep, 1, 1, 1, 1)
        self.gridLayout_2.addWidget(self.animationSettings, 2, 0, 1, 1)
        self.animationControl = QtGui.QHBoxLayout()
        self.animationControl.setObjectName(_fromUtf8("animationControl"))
        self.backwardButton = QtGui.QToolButton(Form)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.backwardButton.sizePolicy().hasHeightForWidth())
        self.backwardButton.setSizePolicy(sizePolicy)
        self.backwardButton.setAccessibleName(_fromUtf8(""))
        self.backwardButton.setText(_fromUtf8(""))
        self.backwardButton.setArrowType(QtCore.Qt.LeftArrow)
        self.backwardButton.setObjectName(_fromUtf8("backwardButton"))
        self.animationControl.addWidget(self.backwardButton)
        self.playButton = QtGui.QPushButton(Form)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.playButton.sizePolicy().hasHeightForWidth())
        self.playButton.setSizePolicy(sizePolicy)
        self.playButton.setObjectName(_fromUtf8("playButton"))
        self.animationControl.addWidget(self.playButton)
        self.forwardButton = QtGui.QToolButton(Form)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.forwardButton.sizePolicy().hasHeightForWidth())
        self.forwardButton.setSizePolicy(sizePolicy)
        self.forwardButton.setText(_fromUtf8(""))
        self.forwardButton.setArrowType(QtCore.Qt.RightArrow)
        self.forwardButton.setObjectName(_fromUtf8("forwardButton"))
        self.animationControl.addWidget(self.forwardButton)
        self.gridLayout_2.addLayout(self.animationControl, 4, 0, 1, 1)
        self.info = QtGui.QGroupBox(Form)
        self.info.setObjectName(_fromUtf8("info"))
        self.gridLayout_4 = QtGui.QGridLayout(self.info)
        self.gridLayout_4.setObjectName(_fromUtf8("gridLayout_4"))
        self.solutionFileLoaded_display = QtGui.QLineEdit(self.info)
        self.solutionFileLoaded_display.setEnabled(False)
        self.solutionFileLoaded_display.setObjectName(_fromUtf8("solutionFileLoaded_display"))
        self.gridLayout_4.addWidget(self.solutionFileLoaded_display, 1, 0, 1, 1)
        self.loadedSolution_label = QtGui.QLabel(self.info)
        self.loadedSolution_label.setObjectName(_fromUtf8("loadedSolution_label"))
        self.gridLayout_4.addWidget(self.loadedSolution_label, 0, 0, 1, 1)
        self.gridLayout_2.addWidget(self.info, 3, 0, 1, 1)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.simulationStartButton = QtGui.QPushButton(Form)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.simulationStartButton.sizePolicy().hasHeightForWidth())
        self.simulationStartButton.setSizePolicy(sizePolicy)
        self.simulationStartButton.setObjectName(_fromUtf8("simulationStartButton"))
        self.horizontalLayout.addWidget(self.simulationStartButton)
        self.simulationStopButton = QtGui.QPushButton(Form)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.simulationStopButton.sizePolicy().hasHeightForWidth())
        self.simulationStopButton.setSizePolicy(sizePolicy)
        self.simulationStopButton.setObjectName(_fromUtf8("simulationStopButton"))
        self.horizontalLayout.addWidget(self.simulationStopButton)
        self.simulationResetButton = QtGui.QPushButton(Form)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.simulationResetButton.sizePolicy().hasHeightForWidth())
        self.simulationResetButton.setSizePolicy(sizePolicy)
        self.simulationResetButton.setObjectName(_fromUtf8("simulationResetButton"))
        self.horizontalLayout.addWidget(self.simulationResetButton)
        self.gridLayout_2.addLayout(self.horizontalLayout, 0, 0, 1, 1)
        self.integrationSettings = QtGui.QGroupBox(Form)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.integrationSettings.sizePolicy().hasHeightForWidth())
        self.integrationSettings.setSizePolicy(sizePolicy)
        self.integrationSettings.setObjectName(_fromUtf8("integrationSettings"))
        self.gridLayout = QtGui.QGridLayout(self.integrationSettings)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.line_3 = QtGui.QFrame(self.integrationSettings)
        self.line_3.setFrameShape(QtGui.QFrame.HLine)
        self.line_3.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_3.setObjectName(_fromUtf8("line_3"))
        self.gridLayout.addWidget(self.line_3, 3, 2, 1, 3)
        self.absTol = QtGui.QLineEdit(self.integrationSettings)
        self.absTol.setEnabled(True)
        self.absTol.setInputMethodHints(QtCore.Qt.ImhDigitsOnly|QtCore.Qt.ImhPreferNumbers)
        self.absTol.setObjectName(_fromUtf8("absTol"))
        self.gridLayout.addWidget(self.absTol, 6, 4, 1, 1)
        self.line_4 = QtGui.QFrame(self.integrationSettings)
        self.line_4.setFrameShape(QtGui.QFrame.HLine)
        self.line_4.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_4.setObjectName(_fromUtf8("line_4"))
        self.gridLayout.addWidget(self.line_4, 10, 2, 1, 3)
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem, 1, 3, 1, 1)
        self.integrationMethodComboBox = QtGui.QComboBox(self.integrationSettings)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.integrationMethodComboBox.sizePolicy().hasHeightForWidth())
        self.integrationMethodComboBox.setSizePolicy(sizePolicy)
        self.integrationMethodComboBox.setMinimumSize(QtCore.QSize(100, 0))
        self.integrationMethodComboBox.setObjectName(_fromUtf8("integrationMethodComboBox"))
        self.integrationMethodComboBox.addItem(_fromUtf8(""))
        self.integrationMethodComboBox.addItem(_fromUtf8(""))
        self.gridLayout.addWidget(self.integrationMethodComboBox, 1, 4, 1, 1)
        self.HminLabel = QtGui.QLabel(self.integrationSettings)
        self.HminLabel.setObjectName(_fromUtf8("HminLabel"))
        self.gridLayout.addWidget(self.HminLabel, 5, 2, 1, 1)
        self.relTol = QtGui.QLineEdit(self.integrationSettings)
        self.relTol.setEnabled(True)
        self.relTol.setInputMethodHints(QtCore.Qt.ImhDigitsOnly|QtCore.Qt.ImhPreferNumbers)
        self.relTol.setObjectName(_fromUtf8("relTol"))
        self.gridLayout.addWidget(self.relTol, 7, 4, 1, 1)
        self.label = QtGui.QLabel(self.integrationSettings)
        self.label.setObjectName(_fromUtf8("label"))
        self.gridLayout.addWidget(self.label, 0, 2, 1, 1)
        self.absTolLabel = QtGui.QLabel(self.integrationSettings)
        self.absTolLabel.setObjectName(_fromUtf8("absTolLabel"))
        self.gridLayout.addWidget(self.absTolLabel, 6, 2, 1, 1)
        self.Hmin = QtGui.QLineEdit(self.integrationSettings)
        self.Hmin.setEnabled(True)
        self.Hmin.setInputMethodHints(QtCore.Qt.ImhDigitsOnly|QtCore.Qt.ImhPreferNumbers)
        self.Hmin.setObjectName(_fromUtf8("Hmin"))
        self.gridLayout.addWidget(self.Hmin, 5, 4, 1, 1)
        self.analysisTypeComboBox = QtGui.QComboBox(self.integrationSettings)
        self.analysisTypeComboBox.setObjectName(_fromUtf8("analysisTypeComboBox"))
        self.analysisTypeComboBox.addItem(_fromUtf8(""))
        self.analysisTypeComboBox.addItem(_fromUtf8(""))
        self.gridLayout.addWidget(self.analysisTypeComboBox, 0, 4, 1, 1)
        self.endTime = QtGui.QLineEdit(self.integrationSettings)
        self.endTime.setInputMethodHints(QtCore.Qt.ImhDigitsOnly|QtCore.Qt.ImhPreferNumbers)
        self.endTime.setObjectName(_fromUtf8("endTime"))
        self.gridLayout.addWidget(self.endTime, 2, 4, 1, 1)
        self.relTolLabel = QtGui.QLabel(self.integrationSettings)
        self.relTolLabel.setObjectName(_fromUtf8("relTolLabel"))
        self.gridLayout.addWidget(self.relTolLabel, 7, 2, 1, 1)
        self.endTimeLabel = QtGui.QLabel(self.integrationSettings)
        self.endTimeLabel.setObjectName(_fromUtf8("endTimeLabel"))
        self.gridLayout.addWidget(self.endTimeLabel, 2, 2, 1, 1)
        self.endTimeLabel_3 = QtGui.QLabel(self.integrationSettings)
        self.endTimeLabel_3.setObjectName(_fromUtf8("endTimeLabel_3"))
        self.gridLayout.addWidget(self.endTimeLabel_3, 1, 2, 1, 1)
        self.loadSolutionFileStatus = QtGui.QCheckBox(self.integrationSettings)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.loadSolutionFileStatus.sizePolicy().hasHeightForWidth())
        self.loadSolutionFileStatus.setSizePolicy(sizePolicy)
        self.loadSolutionFileStatus.setMinimumSize(QtCore.QSize(200, 0))
        self.loadSolutionFileStatus.setObjectName(_fromUtf8("loadSolutionFileStatus"))
        self.gridLayout.addWidget(self.loadSolutionFileStatus, 11, 2, 1, 1)
        self.HmaxLabel = QtGui.QLabel(self.integrationSettings)
        self.HmaxLabel.setObjectName(_fromUtf8("HmaxLabel"))
        self.gridLayout.addWidget(self.HmaxLabel, 4, 2, 1, 1)
        self.Hmax = QtGui.QLineEdit(self.integrationSettings)
        self.Hmax.setEnabled(True)
        self.Hmax.setInputMethodHints(QtCore.Qt.ImhDigitsOnly|QtCore.Qt.ImhPreferNumbers)
        self.Hmax.setObjectName(_fromUtf8("Hmax"))
        self.gridLayout.addWidget(self.Hmax, 4, 4, 1, 1)
        self.TOL_dq_i = QtGui.QLineEdit(self.integrationSettings)
        self.TOL_dq_i.setEnabled(True)
        self.TOL_dq_i.setInputMethodHints(QtCore.Qt.ImhDigitsOnly|QtCore.Qt.ImhPreferNumbers)
        self.TOL_dq_i.setObjectName(_fromUtf8("TOL_dq_i"))
        self.gridLayout.addWidget(self.TOL_dq_i, 8, 4, 1, 1)
        self.TOL_dq_i_Label = QtGui.QLabel(self.integrationSettings)
        self.TOL_dq_i_Label.setObjectName(_fromUtf8("TOL_dq_i_Label"))
        self.gridLayout.addWidget(self.TOL_dq_i_Label, 8, 2, 1, 1)
        self.TOL_C_Label = QtGui.QLabel(self.integrationSettings)
        self.TOL_C_Label.setObjectName(_fromUtf8("TOL_C_Label"))
        self.gridLayout.addWidget(self.TOL_C_Label, 9, 2, 1, 1)
        self.TOL_C = QtGui.QLineEdit(self.integrationSettings)
        self.TOL_C.setEnabled(True)
        self.TOL_C.setInputMethodHints(QtCore.Qt.ImhDigitsOnly|QtCore.Qt.ImhPreferNumbers)
        self.TOL_C.setObjectName(_fromUtf8("TOL_C"))
        self.gridLayout.addWidget(self.TOL_C, 9, 4, 1, 1)
        self.gridLayout_2.addWidget(self.integrationSettings, 1, 0, 1, 1)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(_translate("Form", "Simulation Control", None))
        self.animationSettings.setTitle(_translate("Form", "Animation Settings", None))
        self.updateDisplayStepLabel.setText(_translate("Form", "Update disp. on i-th stp", None))
        self.numberOfSteps_label_2.setText(_translate("Form", "Playback speed", None))
        self.numberOfSteps_label.setText(_translate("Form", "Number of steps", None))
        self.currentStepLabel.setText(_translate("Form", "Current step ", None))
        self.playButton.setText(_translate("Form", "Play", None))
        self.info.setTitle(_translate("Form", "Info", None))
        self.loadedSolution_label.setText(_translate("Form", "Loaded solution", None))
        self.simulationStartButton.setText(_translate("Form", "Start", None))
        self.simulationStopButton.setText(_translate("Form", "Stop", None))
        self.simulationResetButton.setText(_translate("Form", "Reset", None))
        self.integrationSettings.setTitle(_translate("Form", "Simulation Settings", None))
        self.integrationMethodComboBox.setItemText(0, _translate("Form", "Euler", None))
        self.integrationMethodComboBox.setItemText(1, _translate("Form", "Runge-Kutta", None))
        self.HminLabel.setText(_translate("Form", "Hmin", None))
        self.label.setText(_translate("Form", "Analysis type", None))
        self.absTolLabel.setText(_translate("Form", "AbsTol", None))
        self.analysisTypeComboBox.setItemText(0, _translate("Form", "Kinematic", None))
        self.analysisTypeComboBox.setItemText(1, _translate("Form", "Dynamic", None))
        self.relTolLabel.setText(_translate("Form", "RelTol", None))
        self.endTimeLabel.setText(_translate("Form", "End time, s", None))
        self.endTimeLabel_3.setText(_translate("Form", "Integration method", None))
        self.loadSolutionFileStatus.setText(_translate("Form", "Load solution when finished", None))
        self.HmaxLabel.setText(_translate("Form", "Hmax", None))
        self.TOL_dq_i_Label.setText(_translate("Form", "TOL dq_i", None))
        self.TOL_C_Label.setText(_translate("Form", "TOL C(q, t)", None))


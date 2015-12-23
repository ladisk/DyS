# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'analysis.ui'
#
# Created: Wed Sep 30 20:26:27 2015
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
        Form.setWindowModality(QtCore.Qt.NonModal)
        Form.resize(1000, 496)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Form.sizePolicy().hasHeightForWidth())
        Form.setSizePolicy(sizePolicy)
        Form.setMinimumSize(QtCore.QSize(600, 400))
        Form.setMaximumSize(QtCore.QSize(16777215, 16777215))
        Form.setMouseTracking(True)
        Form.setFocusPolicy(QtCore.Qt.NoFocus)
        Form.setAcceptDrops(False)
        self.gridLayout_2 = QtGui.QGridLayout(Form)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setSizeConstraint(QtGui.QLayout.SetNoConstraint)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.splitter = QtGui.QSplitter(Form)
        self.splitter.setOrientation(QtCore.Qt.Horizontal)
        self.splitter.setObjectName(_fromUtf8("splitter"))
        self.treeView = QtGui.QTreeView(self.splitter)
        self.treeView.setMinimumSize(QtCore.QSize(100, 200))
        self.treeView.setAutoFillBackground(False)
        self.treeView.setObjectName(_fromUtf8("treeView"))
        self.treeView.header().setSortIndicatorShown(False)
        self.tableView = QtGui.QTableView(self.splitter)
        self.tableView.setMinimumSize(QtCore.QSize(200, 0))
        self.tableView.setObjectName(_fromUtf8("tableView"))
        self.verticalLayout.addWidget(self.splitter)
        self.groupBox = QtGui.QGroupBox(Form)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox.sizePolicy().hasHeightForWidth())
        self.groupBox.setSizePolicy(sizePolicy)
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.gridLayout = QtGui.QGridLayout(self.groupBox)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.progressBar = QtGui.QProgressBar(self.groupBox)
        self.progressBar.setProperty("value", 24)
        self.progressBar.setObjectName(_fromUtf8("progressBar"))
        self.gridLayout.addWidget(self.progressBar, 2, 0, 1, 7)
        self.number_of_simulations = QtGui.QLineEdit(self.groupBox)
        self.number_of_simulations.setObjectName(_fromUtf8("number_of_simulations"))
        self.gridLayout.addWidget(self.number_of_simulations, 1, 1, 1, 1)
        self.startButton = QtGui.QPushButton(self.groupBox)
        self.startButton.setObjectName(_fromUtf8("startButton"))
        self.gridLayout.addWidget(self.startButton, 1, 4, 1, 1)
        self.label = QtGui.QLabel(self.groupBox)
        self.label.setObjectName(_fromUtf8("label"))
        self.gridLayout.addWidget(self.label, 1, 0, 1, 1)
        self.loadButton = QtGui.QPushButton(self.groupBox)
        self.loadButton.setObjectName(_fromUtf8("loadButton"))
        self.gridLayout.addWidget(self.loadButton, 0, 4, 1, 1)
        self.closeButton = QtGui.QPushButton(self.groupBox)
        self.closeButton.setObjectName(_fromUtf8("closeButton"))
        self.gridLayout.addWidget(self.closeButton, 1, 5, 1, 1)
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem, 1, 2, 1, 1)
        self.label_2 = QtGui.QLabel(self.groupBox)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.gridLayout.addWidget(self.label_2, 0, 0, 1, 1)
        spacerItem1 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem1, 1, 3, 1, 1)
        self.verticalLayout.addWidget(self.groupBox)
        self.gridLayout_2.addLayout(self.verticalLayout, 0, 0, 1, 1)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(_translate("Form", "Analysis", None))
        self.startButton.setText(_translate("Form", "Start", None))
        self.label.setText(_translate("Form", "Number of simulations ", None))
        self.loadButton.setText(_translate("Form", "Load", None))
        self.closeButton.setText(_translate("Form", "Close", None))
        self.label_2.setText(_translate("Form", "Load parameters from file", None))


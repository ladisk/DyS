# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'output_widget.ui'
#
# Created: Thu Jul 03 17:21:41 2014
#      by: PyQt4 UI code generator 4.9.6
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
        Form.resize(500, 217)
        Form.setMinimumSize(QtCore.QSize(500, 200))
        Form.setCursor(QtGui.QCursor(QtCore.Qt.IBeamCursor))
        Form.setAutoFillBackground(True)
        self.gridLayout = QtGui.QGridLayout(Form)
        self.gridLayout.setSizeConstraint(QtGui.QLayout.SetNoConstraint)
        self.gridLayout.setMargin(0)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.textOutputWidget = QtGui.QTextBrowser(Form)
        self.textOutputWidget.setMinimumSize(QtCore.QSize(500, 200))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Consolas"))
        font.setPointSize(9)
        font.setKerning(False)
        self.textOutputWidget.setFont(font)
        self.textOutputWidget.viewport().setProperty("cursor", QtGui.QCursor(QtCore.Qt.IBeamCursor))
        self.textOutputWidget.setAutoFillBackground(True)
        self.textOutputWidget.setObjectName(_fromUtf8("textOutputWidget"))
        self.gridLayout.addWidget(self.textOutputWidget, 0, 0, 1, 1)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(_translate("Form", "Output window", None))


__author__ = 'lskrinjar'
import numpy as np
from pprint import pprint
from PyQt4 import QtCore, QtGui
from tree_view_widget.tree_model import TreeModel
from MBD_system.MBD_system_items import AnalysisTreeItem


class AnalysisListModel(QtCore.QAbstractListModel):#QtCore.QAbstractTableModel, TreeModel, QtCore.QAbstractListModel

    """INPUTS: Node, QObject"""
    def __init__(self, data, parent=None):
        super(AnalysisListModel, self).__init__(parent)
        self._parent = parent
        self._data = data
        print "ListModel()"
        print "self._data =", self._data
        # #   list through all parameters
        # for i in range(0, 10):
        #     AnalysisTreeItem("mass", 1, parent=self._rootNode)

        self.headers = ["Parameter", "Simulation status", "Item", "Mean value", "Random value"]

    def headerData(self, column, orientation, role):
        """
        INPUTS: int, Qt::Orientation, int
        OUTPUT: QVariant, strings are cast to QString which is a QVariant
        """
        self.headers = ["Parameter", "Simulation status", "Item", "Mean value", "Random value"]
        # print "column =", column, "orientation =", orientation, "role =", role
        if role == QtCore.Qt.DisplayRole:
            #    column headers
            if orientation == QtCore.Qt.Horizontal:
                return self.headers[column]
            #    row headers
            if orientation == QtCore.Qt.Vertical:
                # self.setRowHeight(column, 10)
                return column
            # # if section == 0:
            # #     return "Selected"
            # if section == 0:
            #     return "Parameter"
            # elif section == 1:
            #     return "Item"
            # elif section == 2:
            #     return "Mean value"
            # elif section == 3:
            #     return "Min"
            # elif section == 4:
            #     return "Max"
            # else:
            #     return "Typeinfo"

    def flags(self, index):

        item = index.internalPointer()
        #   non editable paramater name and object reference
        if index.column() == 0 or index.column() == 1:
            return QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsUserCheckable
        #   editable parameter mean, min, max values
        else:
            return QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEditable


    """INPUTS: QModelIndex"""
    """OUTPUT: int"""
    def columnCount(self, parent):
        """
        Number of columns
        :param parent:
        :return:
        """
        return len(self.headers)

    def rowCount(self, parent):
        """
        
        """
        return len(self._data)

    def data(self, index, role):
        """
        INPUTS: QModelIndex, int
        QVariant, strings are cast to QString which is a QVariant
        :param index:
        :param role:
        :return:
        """
        if not index.isValid():
            return None

        item = index.internalPointer()
        print "item =", item
        # pprint(vars(item))

        # brush = QtGui.QBrush(QtGui.QColor(255, 0, 0))
        # brush.setStyle(QtCore.Qt.SolidPattern)
        # item.setBackground(brush)

        #   background color
        if role == QtCore.Qt.BackgroundRole:
            if index.row() % 2 == 0:
                return QtGui.QBrush(QtGui.QColor(240, 240, 240))
            return QtGui.QBrush(QtCore.Qt.transparent)

        # if role == QtCore.Qt.CheckStateRole:
        #     # if self.dataObj[index.row()] == 0:
        #     return QtCore.QVariant(QtCore.Qt.Unchecked)
        # else:
        #     return QtCore.QVariant(QtCore.Qt.Checked)
        
        if role == QtCore.Qt.DisplayRole or role == QtCore.Qt.EditRole:
            #   row number
            if index.column() == 0:
                return "Test"#, item.name()
            #   object
            if index.column() == 1:
                try:
                    return item._item.name()
                except:
                    return None
            #   mean(nominal) value
            if index.column() == 2:
                return str(item.mean_value())
            #   min limit value
            if index.column() == 3:
                return str(item.min_value())
            #   max limit value
            if index.column() == 4:
                return str(item.max_value())

        # if role == QtCore.Qt.CheckStateRole:
        #     if index.column() == 0:
        #         # print "item._selected =", item._selected
        #         if item._selected:
        #             return QtCore.Qt.Checked
        #         else:
        #             return QtCore.Qt.Unchecked

        # if role == Qt.CheckStateRole:
        #     row = index.row()
        #     print self.args[row].checked
        #     if self.args[row].checked == False:
        #         return QVariant(Qt.Unchecked)
        #     else:
        #         return QVariant(Qt.Checked)

    def setData(self, index, value, role):

        item = index.internalPointer()

        if role == QtCore.Qt.CheckStateRole:
            row = index.row()
            item._selected = not item._selected
        return True
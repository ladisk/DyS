__author__ = 'lskrinjar'
from pprint import pprint
from PyQt4 import QtCore, QtGui
from tree_view_widget.tree_model import TreeModel
from MBD_system.MBD_system_items import AnalysisTreeItem


class AnalysisTreeModel(TreeModel):
    """
    INPUTS: Node, QObject
    """
    def __init__(self, MBD_system, parent=None):
        super(AnalysisTreeModel, self).__init__(parent)
        self._parent = parent
        self._rootNode = MBD_system
        print "TreeModel()"
        print "self._rootNode =", self._rootNode
        #   column headers
        self.headers = ["Parameter", "Item (Name)", "Mean value", "Min", "Max"]
        # #   list through all parameters
        # for i in range(0, 10):
        #     AnalysisTreeItem("mass", 1, parent=self._rootNode)

    def headerData(self, section, orientation, role):
        """
        INPUTS: int, Qt::Orientation, int
        OUTPUT: QVariant, strings are cast to QString which is a QVariant
        """
        if role == QtCore.Qt.DisplayRole:
            return self.headers[section]

        if role == QtCore.Qt.SizeHintRole:
            if section == 0:
                pass
                # self.setColumnWidth()
                # self.setColumnWidth(20)
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
        """
        
        """
        item = index.internalPointer()
        #   non editable parameter name and object reference
        if index.column() == 0 or index.column() == 1 or index.column() == 2:
            return QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsUserCheckable
        #   editable parameter mean, min, max values
        else:
            return QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEditable

    def rowCount(self, parent):
        """
        INPUTS: QModelIndex
        OUTPUT: int
        """
        if not parent.isValid():
            parentNode = self._rootNode
        else:
            parentNode = parent.internalPointer()

        return parentNode.childCount()

    def columnCount(self, parent):
        """
        Number of columns
        :param parent:
        :return:
        INPUTS: QModelIndex
        OUTPUT: int
        """
        return len(self.headers)

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

        # brush = QtGui.QBrush(QtGui.QColor(255, 0, 0))
        # brush.setStyle(QtCore.Qt.SolidPattern)
        # item.setBackground(brush)

        #   backgrund color
        # if role == QtCore.Qt.BackgroundRole:
        #     if index.row() % 2 == 0:
        #         return QtGui.QBrush(QtGui.QColor(240, 240, 240))
        #     return QtGui.QBrush(QtCore.Qt.transparent)

        # if role == QtCore.Qt.CheckStateRole:
        #     # if self.dataObj[index.row()] == 0:
        #     return QtCore.QVariant(QtCore.Qt.Unchecked)
        # else:
        #     return QtCore.QVariant(QtCore.Qt.Checked)
        
        if role == QtCore.Qt.DisplayRole or role == QtCore.Qt.EditRole:
            #   parameter name
            if index.column() == 0:
                return item.name()
            #   object
            if index.column() == 1:
#                 print item
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

        if role == QtCore.Qt.CheckStateRole:
            if index.column() == 0:
                # print "item._selected =", item._selected
                if item._selected:
                    print "self._parent.MBD_system._children[0].parameters_selected =", self._parent.MBD_system._children[0].parameters_selected
                    self._parent.MBD_system._children[0].parameters_selected.append(item)
                    return QtCore.Qt.Checked
                else:
                    return QtCore.Qt.Unchecked

        # if role == Qt.CheckStateRole:
        #     row = index.row()
        #     print self.args[row].checked
        #     if self.args[row].checked == False:
        #         return QVariant(Qt.Unchecked)
        #     else:
        #         return QVariant(Qt.Checked)

    def setData(self, index, value, role):
        """
        INPUTS: QModelIndex, QVariant, int (flag)
        """
        item = index.internalPointer()
        if index.isValid():
            if role == QtCore.Qt.CheckStateRole:
                row = index.row()
                item._selected = not item._selected
                return True
            if role == QtCore.Qt.DisplayRole or role == QtCore.Qt.EditRole:
                #    set min val
                if index.column() == 3:
                    item._min_value = value.toFloat()[0]
                    return True
                #    set max val
                if index.column() == 4:
                    item._max_value = value.toFloat()[0]
                    return True
        return False

    def index(self, row, column, parent):
        """INPUTS: int, int, QModelIndex
        OUTPUT: QModelIndex
        Should return a QModelIndex that corresponds to the given row, column and parent node"""
        parentNode = self.getNode(parent)

        childItem = parentNode.child(row)

        if childItem:
            return self.createIndex(row, column, childItem)
        else:
            print QtCore.QModelIndex()
            return QtCore.QModelIndex()

    def parent(self, index):
        """
        Should return the parent of the node with the given QModelIndex
        :param QModelIndex:
        :return QModelIndex:
        """
        node = self.getNode(index)
        parentNode = node.parent()

        if parentNode == self._rootNode:
            return QtCore.QModelIndex()

        return self.createIndex(parentNode.row(), 0, parentNode)

    def getNode(self, index):
        """
        CUSTOM
        INPUTS: QModelIndex
        """
        if index.isValid():
            node = index.internalPointer()
            if node:
                return node
        return self._rootNode

    def insertRows(self, position, rows, parent=QtCore.QModelIndex()):
        """
        INPUTS: int, int, QModelIndex
        """
        parentNode = self.getNode(parent)

        self.beginInsertRows(parent, position, position + rows - 1)

        for row in range(rows):

            childCount = parentNode.childCount()
            childNode = TreeItem("untitled" + str(childCount))
            success = parentNode.insertChild(position, childNode)

        self.endInsertRows()
        return success

    def insertRow(self, position, child, parent=QtCore.QModelIndex()):
        """

        """
        self.beginInsertRows(QtCore.QModelIndex(), parent.row() + len(parent._children), parent.row() + len(parent._children) + 1)

        success = parent.insertChild(QtCore.QModelIndex().row() + len(parent._children) + 1, child)

        self.endInsertRows()

        return success

    def removeRows(self, position, rows, parent=QtCore.QModelIndex()):
        """
        INPUTS: int, int, QModelIndex
        """
        parentNode = self.getNode(parent)
        self.beginRemoveRows(parent, position, position + rows - 1)

        for row in range(rows):
            success = parentNode.removeChild(position)

        self.endRemoveRows()

        return success


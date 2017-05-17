__author__ = 'lskrinjar'


from PyQt4 import QtCore, QtGui


class TreeModel(QtCore.QAbstractItemModel):
    """INPUTS: Node, QObject"""
    def __init__(self, root, parent=None):
        super(TreeModel, self).__init__(parent)
        self._rootNode = root
        
        self._editable_types = ["body", "force", "joint", "contact", "spring", "motion", "measure", "variable"]

    def headerData(self, section, orientation, role):
        """
        INPUTS: int, Qt::Orientation, int
        OUTPUT: QVariant, strings are cast to QString which is a QVariant
        """
        if role == QtCore.Qt.DisplayRole:
            if section == 0:
                return self._rootNode.name()
            else:
                return "Typeinfo"

    def flags(self, index):
        """

        :param index: QModelIndex
        :return flag: int
        """
        item = index.internalPointer()

        if item.typeInfo() in self._editable_types:# or item._name != "Ground":
            return QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsEditable
        #     if item._children == []:
        #         return QtCore.Qt.ItemIsSelectable

        return QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsEditable

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
        return 1

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
        # print "item =", item
        #   background color
        # if role == QtCore.Qt.BackgroundRole:
        #         return QtGui.QBrush(QtGui.QColor(240, 240, 240))

        if role == QtCore.Qt.DisplayRole or role == QtCore.Qt.EditRole:
            if index.column() == 0:
                return item.name()
            # if index.column() == 1:
            #     return node.typeInfo()

        #     add icons here
        if role == QtCore.Qt.DecorationRole:

            if index.column() == 0:
                typeInfo = item.typeInfo()

                if typeInfo == "LIGHT":
                    return QtGui.QIcon(QtGui.QPixmap(":/Light.png"))

                if typeInfo == "TRANSFORM":
                    return QtGui.QIcon(QtGui.QPixmap(":/Transform.png"))

                if typeInfo == "CAMERA":
                    return QtGui.QIcon(QtGui.QPixmap(":/Camera.png"))

                if typeInfo == "group":
                    return QtGui.qApp.style().standardIcon(QtGui.QStyle.SP_DirClosedIcon)
#                    return QtGui.QIcon(":/folder.png")

                if typeInfo == "settings":
                    return QtGui.qApp.style().standardIcon(QtGui.QStyle.SP_DirClosedIcon)

                if typeInfo == "solution":
                    return QtGui.qApp.style().standardIcon(QtGui.QStyle.SP_DirClosedIcon)
#                    print type(QtGui.QStyle.SP_DirIcon)
#                    __icon = QtGui.QIcon()
#                    __icon.addPixmap()
#                    style = QtGui.QStyle("")
#                    __icon.addPixmap(QtGui.QStyle.standardPixmap(QtGui.QStyle.SP_DirClosedIcon),QtGui.QIcon.Normal, QtGui.QIcon.Off)
#                    __icon.standardIcon("SP_DirIcon")
#                    return __icon
#                    return QtGui.QStyle.SP_DirIcon
#                    return QtGui.QStandardItem(QtGui.QStyle.SP_DirIcon)
#                    return QtGui.QIcon.fromTheme("folder-new")
#                    return QtGui.QStyle.standardIcon(QtGui.QStyle.SP_DirClosedIcon)
#                    return QtGui.QIcon(QtGui.QPixmap("c:/Users/reti-luka/Dropbox/DyS/src/V55/icons/folder.svg"))
#                    return QtGui.QIcon(QtGui.QStyle.SP_DirClosedIcon)#(QtGui.QPixmap("c:/Users/reti-luka/Dropbox/DyS/src/V55/icons/folder.svg"))
#                    QtGui.QStyle.SP_DirClosedIcon

    def setData(self, index, value, role=QtCore.Qt.EditRole):
        """
        INPUTS: QModelIndex, QVariant, int (flag)
        """
        if index.isValid():
            if role == QtCore.Qt.EditRole:
                item = index.internalPointer()
                item.setName(value)
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
            return QtCore.QModelIndex()

    def parent(self, index):
        """
        Should return the parent of the node with the given QModelIndex
        :param QModelIndex:
        :return QModelIndex:
        """
        node = self.getNode(index)
        parentNode = node.parent()

        if parentNode == self._rootNode or parentNode.row() is None:
            return QtCore.QModelIndex()

        # if not hasattr(parentNode, "row"):

        if parentNode.row() is None:
            print "-----------------"
            print node._name
            print parentNode

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
#         print "treeView.model().insertRow()"
# 
#         print "child =", child, id(child)
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
#     def createGroups(self):
#         print "EXE?"
#         print self._rootNode.typeInfo()
#         if self._rootNode.typeInfo().lower() == "mbdsystem":
#             None
# #             for group in self._rootNode.list_of_object_groups:
#                 print "group =", group
#                 group_obj = TreeItem(group, rootNode)
#                 print "------------------------"
#                 for obj in getattr(self._rootNode, group.lower()):
#                     print obj
#                     group_obj.addChild(obj)
#                 group_obj
#                 if getattr(self._rootNode, group.lower()) != []:
#                     for obj in getattr(self._rootNode, group.lower()):
#                         if group
#                         print obj.typeInfo(), obj._name

#                     print getattr(self._rootNode, group.lower())
#                     for item in getattr(self._rootNode, group.lower()):
#                         print "item =", item
#                         TreeItem(item, group)





if __name__ == "__main__":
#     print os.listdir("..//dynamic_systems//dynamic_system_0")


#     pprint(vars(a))
#     a.addBody("body_1")
#     pprint(vars(a))
#     for joint in a.joints:
#         pprint(vars(joint))
#     print a.create_q0()
#     print a.bodies[1].CM[0:2]
#     print a.forces
    app = QtGui.QApplication(sys.argv)
    app.setStyle("plastique")

    rootNode = TreeItem("rootNode")
    MBDNode = MBDsystem(MBD_folder_abs_path="..//dynamic_systems//dynamic_system_0", parent=rootNode)
#     a = MBDsystem(MBD_folder_abs_path="..//dynamic_systems//dynamic_system_0")
#     rootNode.addChild(a)
#     print "rootNode ="
#     pprint(vars(a))
#     childNode0 = TransformNode("RightPirateLeg", rootNode)
#     childNode1 = TreeItem("RightPirateLeg_END", rootNode)
#     childNode2 = CameraNode("LeftFemur", rootNode)
#     childNode3 = Node("LeftTibia", childNode2)
#     childNode4 = Node("LeftFoot", rootNode)
#     childNode5 = LightNode("LeftFoot_END", childNode4)

#     print "rootNode =", rootNode
#     print "childNode1 =", childNode1

    model = TreeModel(rootNode)
#     model.createGroups()
#     root_index = model.setRootPath()

    treeView = QtGui.QTreeView()
    treeView.show()
    treeView.setModel(model)
    treeView.setRootIsDecorated(False)
    treeView.expandAll()
#     treeView.setRootIndex(rootNode.index)
#     rightPirateLeg = model.index(0, 0, QtCore.QModelIndex())

    sys.exit(app.exec_())
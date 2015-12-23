'''
Created on 27. jan. 2014

@author: lskrinjar
'''
import os
import sys
import subprocess
from pprint import pprint
import numpy as np

from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from tree_view_widget_ui import Ui_Form
from MBD_system.solution_data.solution_data import SolutionData
from MBD_system.body.body_widget import BodyWidget
from MBD_system.force.force_widget import ForceWidget
from MBD_system.joint.joint_widget import JointWidget
from MBD_system.contact.contact_widget import ContactWidget
from MBD_system.spring.spring_widget import SpringWidget
from MBD_system.item_widget import ItemWidget
from tree_model import TreeModel
from MBD_system.MBD_system_items import SolutionDataItem
from analysis.analysis_widget import AnalysisWidget
from MBD_system.MBD_system_widget import MBDSystemWidget

class solutionFilenameSignal(QtCore.QObject):
    signal_filename = QtCore.pyqtSignal(str, name='')

class loadSolutionFile(QtCore.QObject):
    signal_loadSolutionFile = QtCore.pyqtSignal()

class TreeViewWidget(QWidget):  # QMainWindow#, QAbstractItemView 
    """
    classdocs
    """
    def __init__(self, MBD_system_=None, parent=None, flags=0):
        """
        Constructor
        """
        super(TreeViewWidget, self).__init__(parent)
        self.ui = Ui_Form()
        self.ui.setupUi(self)
        
        self.filename_signal = solutionFilenameSignal()
        self.load_loadSolutionFile = loadSolutionFile()
        
        self._parent = parent

        self.MBD_system = MBD_system_
        # _model = MBD_system.TreeModel(self.MBD_system)
        _model = TreeModel(self.MBD_system)

        self.ui.treeView.setModel(_model)
        
        self.ui.treeView.setRootIsDecorated(True)
        self.ui.treeView.setHeaderHidden(True)
        self.ui.treeView.expandAll()

        # we want our listview to have a context self.menu taken from the actions on this widget
        # those actions will be to delete an item :)
        self.ui.treeView.customContextMenuRequested.connect(self.contextMenuEvent)
        QtCore.QObject.connect(self.ui.treeView, QtCore.SIGNAL("clicked (QModelIndex)"),self.row_clicked)
 
        self._widget = ItemWidget(parent=self)

    def row_clicked(self, index):
        """
        when a row is clicked... show the name
        """

    def setWindowFlags(self, flags):
        """

        """
        super(TreeViewWidget, self).setWindowFlags(flags)

    @pyqtSlot()
    def onTriggered(self, event):
        """
        
        """
        # tell our model to remove the selected row.
#         self.ui.treeView.model()
  
#         self.model().removeRows(self.currentIndex().row(), 1)
    def create_action(self, event):
        """

        """
        print event.pos()

    def setSelection(self, current, old):
        """

        """
        current = self._proxyModel.mapToSource(current)

    def contextMenuEvent(self, event):
        """
        Context menu
        """
        self.menu = QtGui.QMenu(self)

        indexes = self.ui.treeView.selectedIndexes()

        self._item = self.ui.treeView.model().getNode(indexes[0])

        if len(indexes) > 0:
            level = 0
            index = indexes[0]
            while index.parent().isValid():
                index = index.parent()
                level += 1

        if self._item._typeInfo == "group":
            createAction = self.menu.addAction("Create")
            createAction.triggered.connect(self._create)
            
            self.menu.addSeparator()
            
            loadGroupItemsAction = self.menu.addAction("Load items from file")
            loadGroupItemsAction.triggered.connect(self._load_group_items)
        
            self.menu.addSeparator()
            
        if self._item._typeInfo != "group" or self._item._typeInfo != "solution":# or self._item._name.lower() != "solution":
            editAction = self.menu.addAction("Edit")
            editAction.triggered.connect(self._edit)

            deleteAction = self.menu.addAction("Delete")
            deleteAction.setEnabled(True)
            deleteAction.triggered.connect(self._delete)

            self.menu.addSeparator()

        if self._item._typeInfo == "MBDsystem" or self._item._typeInfo == "settings":
            self.menu.addSeparator()
            show_GCS_Action = QtGui.QAction("Show GCS ", self, checkable=True, checked=self._parent.SimulationControlWidget.OpenGLWidget.GCS._visible)
            show_GCS_Action.triggered.connect(self._parent.SimulationControlWidget.OpenGLWidget.GCS._show)
            self.menu.addAction(show_GCS_Action)

            propertiesAction = self.menu.addAction("Properties")
            propertiesAction.triggered.connect(self._properties)
            
            self.menu.addSeparator()
            
            getDOFsAction = self.menu.addAction("DOF")
            getDOFsAction.triggered.connect(self.__getDOF)

            if self._item._typeInfo == "MBDsystem":
                evaluate_C_Action = self.menu.addAction("C(q, t)")
                evaluate_C_Action.triggered.connect(self._evaluate_C)

                get_q_Action = self.menu.addAction("Get q")
                get_q_Action.triggered.connect(self._get_q)
                
                list_parameters_Action = self.menu.addAction("List parameters")
                list_parameters_Action.triggered.connect(self._list_parameters)

                analysis_Action = self.menu.addAction("Analysis")
                analysis_Action.triggered.connect(self._analysis)
                analysis_Action.setShortcut(QtGui.QKeySequence("Alt+A"))

                self.menu.addSeparator()

        elif self._item._typeInfo == "solution":
            loadAction = self.menu.addAction("Load from File")
            loadAction.triggered.connect(self._load_solution_file_to_project)

        elif self._item._typeInfo == "solutiondata":
            loadAction = self.menu.addAction("Load solution")
            loadAction.triggered.connect(self._load_solution_data)

            propertiesAction = self.menu.addAction("Properties")
            propertiesAction.triggered.connect(self.solution_properties)

        elif self._item._typeInfo == "group":

            propertiesAction = self.menu.addAction("Properties")
            propertiesAction.triggered.connect(self._properties)

#         elif self._item._typeInfo == "solutiondata":
#             
#             openAction = self.menu.addAction("Open")
#             openAction.triggered.connect(self._open)
#             
#             self.menu.addSeparator()
# 
#             deleteAction = self.menu.addAction("Delete")
#             deleteAction.triggered.connect(self._delete)
            
        elif self._item._typeInfo == "force":
            show_force_Action = QtGui.QAction("Show force ", self, checkable=True, checked=self._item._visible)
            show_force_Action.triggered.connect(self._item._show_force)
            self.menu.addAction(show_force_Action)

            self.menu.addSeparator()
            for marker in self._item.markers:
                show_marker_Action = QtGui.QAction("Show u_P on Body "+str(self._item.body_id), self, checkable=True, checked=marker._visible)
                show_marker_Action.triggered.connect(marker._show)
                self.menu.addAction(show_marker_Action)

        elif self._item._typeInfo == "joint" or self._item._typeInfo == "spring":
            for force, ID in zip(self._item.force_list, ["i", "j"]):
                show_force_Action = QtGui.QAction("Show force on Body "+ID, self, checkable=True, checked=force._visible)
                show_force_Action.triggered.connect(force._show_force)
                self.menu.addAction(show_force_Action)

            self.menu.addSeparator()
            for marker, ID in  zip(self._item.markers, ["i", "j"]):
                show_marker_Action = QtGui.QAction("Show u_P on body "+ID, self, checkable=True, checked=marker._visible)
                show_marker_Action.triggered.connect(marker._show)
                self.menu.addAction(show_marker_Action)

        elif self._item._typeInfo == "contact":
            for force, ID in zip(self._item._Fn_list, ["i", "j"]):
                show_force_Action = QtGui.QAction("Show Fn on body "+ID, self, checkable=True, checked=force._visible)
                show_force_Action.triggered.connect(force._show_force)
                self.menu.addAction(show_force_Action)

            self.menu.addSeparator()
            for force, ID in zip(self._item._Ft_list, ["i", "j"]):
                show_force_Action = QtGui.QAction("Show Ft on body "+ID, self, checkable=True, checked=force._visible)
                show_force_Action.triggered.connect(force._show_force)
                self.menu.addAction(show_force_Action)

            self.menu.addSeparator()
            if self._item._type == "general":
                # try:
                for body_id, ID in zip(self._item.body_id_list, ["i", "j"]):
                    _body = self.MBD_system._children[0].bodies[body_id]
                    _AABB = _body.AABBtree

                    show_AABB_action = QtGui.QAction("Show AABB on body "+ID, self, checkable=True, checked=_AABB._visible)
                    show_AABB_action.triggered.connect(_AABB._show_AABB)
                    self.menu.addAction(show_AABB_action)

                self.menu.addSeparator()

                for body_id, ID in zip(self._item.body_id_list, ["i", "j"]):
                    _body = self.MBD_system._children[0].bodies[body_id]
                    _AABB = _body.AABBtree

                    properties_of_AABB_action = self.menu.addAction("Properties of AABB on body "+ID)
                    properties_of_AABB_action.triggered.connect(_AABB._get_properties)

                self.menu.addSeparator()
                # except:
                #     pass
            try:
                for contact_geometry, ID in zip(self._item.contact_geometry_list, ["i", "j"]):
                    # get_nodes_action = QtGui.QAction("Get contact nodes "+ID, self, checkable=True, checked=_AABB._visible)
                    get_nodes_action = self.menu.addAction("Get contact nodes "+ID)
                    get_nodes_action.triggered.connect(contact_geometry.get_nodes)
                    # self.menu.addAction(get_nodes_action)
            except:
                pass

            contactForceAction = self.menu.addAction("Contact force")
            contactForceAction.triggered.connect(self._item.get_contact_force)

            self.menu.addSeparator()

            if self._item._type == "Revolute Clearance Joint":
                for marker in self._item.markers:
                    show_marker_Action = QtGui.QAction("Show Body "+ID+" joint center", self, checkable=True, checked=marker._visible)
                    show_marker_Action.triggered.connect(marker._show)
                    self.menu.addAction(show_marker_Action)

        elif self._item._typeInfo == "body":
            if self._item._name == "Ground" or self._item.body_id == -1:
                pass
            else:
                show_body_Action = QtGui.QAction("Show Body", self, checkable=True, checked=self._item._visible)
                show_body_Action.triggered.connect(self._item._show)
                self.menu.addAction(show_body_Action)

                show_LCS_Action = QtGui.QAction("Show LCS", self, checkable=True, checked=self._item.LCS._visible)
                show_LCS_Action.triggered.connect(self._item.LCS._show)
                self.menu.addAction(show_LCS_Action)

                show_CAD_CS_Action = QtGui.QAction("Show CAD CS", self, checkable=True, checked=self._item.CAD_CS._visible)
                show_CAD_CS_Action.triggered.connect(self._item.CAD_CS._show)
                self.menu.addAction(show_CAD_CS_Action)
                self.menu.addSeparator()

                getBody_q_Action = self.menu.addAction("Get body q")
                getBody_q_Action.triggered.connect(self._item.get_q)

                getBody_dq_Action = self.menu.addAction("Get body dq")
                getBody_dq_Action.triggered.connect(self._item.get_dq)

                getBody_qdq_Action = self.menu.addAction("Get body [q dq]")
                getBody_qdq_Action.triggered.connect(self._item.get_qdq)
                self.menu.addSeparator()

        self.menu.exec_(event.globalPos())
        self._parent.SimulationControlWidget.OpenGLWidget._repaintGL()

    def closeEvent(self, event):
        """

        """
        self._parent.SimulationControlWidget.OpenGLWidget._repaintGL()

    def _list_parameters(self):
        """
        
        """
        self._item._list_parameters(item=self.MBD_system)

    def _analysis(self):
        """
        
        """
        self._widget = AnalysisWidget(self.MBD_system)
        self._widget.show()

    def _load_group_items(self):
        """
        
        """
        print "Under construction!", os.path.realpath(__file__)

    def _load_solution_data(self):
        """

        """
        #   this has to be put in a thread
        data = self._item._load_dat_file()
        self._parent.SimulationControlWidget.load_solution_file(filename = self._item._name, solution_data=data)

    def _load_solution_file_to_project(self):
        """
        Load solution file to opened project if there is any solution file.
        TODO - check if right solution data is loaded to project
        """
        #   filetype filter
        _filter = "Data File (*.dat);;Excel (*.xlsx);;CSV (*.csv)"

        #   open dialog
        _open_file = QtGui.QFileDialog()
        _open_file.setFileMode(QFileDialog.ExistingFiles)
        
        #   directory to open
        _directory = self._item._parent._parent._children[0].MBD_folder_abs_path_
        #    set directory from where to open
        _open_file.setDirectory(_directory)

        #    file abs path and filetype(extension)
        file, filetype = _open_file.getOpenFileNameAndFilter(self, "Load Solution File", _directory, _filter)
        
        #    get filename from abspath
        filename = os.path.basename(str(file))
        
        #   filename without extension
        name, _extension = filename.split(".")

        #   create solution data object
        _solution = SolutionData(name, file=file)

        #   add solution data object item to treeview
        self.ui.treeView.model().insertRow(len(self._item._children), _solution, self._item)

    def _open(self):
        """
        
        """
        subprocess.call(['notepad.exe', str(self.__filename)])

    def add_solution_data(self, solution_data_object_ID=None, filename=None):
        """
        Function adds created solution object item in tree view in folder solutions
        """
        print "solution_data_object_ID =", solution_data_object_ID
        print "filename =", filename
        #   root index
        root_index = self.ui.treeView.rootIndex()

        #   get root node with root index
        self._item = self.ui.treeView.model().getNode(root_index)
        MBD_system_item = self._item._children[0]
        # pprint(vars(MBD_system_item))
        #    solutions (group) item
        _solution_group_item = MBD_system_item._children[6]

        #    add child (solution data) to solution item
        pos = len(_solution_group_item._children)
        _childNode = MBD_system_item.solutions[-1]
        self.ui.treeView.model().insertRow(pos, _childNode, _solution_group_item)

    def solution_properties(self):
        pprint(vars(self._item))

    def _create(self, index):
        """
        
        """
        self._parentNode = self._item
        self._childNode = None
        if self._item._typeInfo == "group":
            if self._item._name.lower() == "bodies":
                self._widget = BodyWidget(self._parentNode, parent=self)

            if self._item._name.lower() == "forces":
                self._widget = ForceWidget(self._parentNode, parent=self)

            if self._item._name.lower() == "joints":
                self._widget = JointWidget(self._parentNode, parent=self)

            if self._item._name.lower() == "contacts":
                self._widget = ContactWidget(self._parentNode, parent=self)

            if self._item._name.lower() == "springs":
                self._widget = SpringWidget(self._parentNode, parent=self)
                # self._childNode = MBD_system_items.SpringItem("Spring_" + str(len(self._parentNode._children)))

            # if self._childNode is not None:
            #     print "created =", self._childNode
            #     if len(self._parentNode._children) == 0:
            #         pos = 1
            #     else:
            #         pos = len(self._parentNode._children)
            #     self.ui.treeView.model().insertRow(pos, self._childNode, self._parentNode)
    
    def _insertTreeItem(self):
        """
        
        """
        pos = len(self._parentNode._children)
        self._childNode = self._widget.item
        self.ui.treeView.model().insertRow(pos, self._childNode, self._parentNode)

    def __getDOF(self):
        """
        
        """
        C_q, C_qT = self._parent.SimulationControlWidget.solver.solveODE.ode_fun.getDOF()

    def _get_q(self):
        """

        """
        q = self.MBD_system._children[0].get_q()
        print "q ="
        print q

    def _evaluate_C(self):
        """

        """
        t = 0
        q_ = self.MBD_system._children[0].get_q()
        print "C(q, t)"
        print self._parent.SimulationControlWidget.solver.solveODE.ode_fun.create_C(t, q_)

    def _edit(self, index):
        """

        """
        self._parentNode = self._item._parent

        if self._item._typeInfo.lower() == "body":
            self._widget = BodyWidget(self._parentNode, parent=self)

        if self._item._typeInfo.lower() == "contact":
            self._widget = ContactWidget(self._parentNode, parent=self)

        if self._item._typeInfo.lower() == "force":
            self._widget = ForceWidget(self._parentNode, parent=self)

        if self._item._typeInfo.lower() == "joint":
            self._widget = JointWidget(self._parentNode, parent=self)

        if self._item._typeInfo.lower() == "spring":
            self._widget = SpringWidget(self._parentNode, parent=self)

        # if self._item._typeInfo.lower() == "MBDsystem":
        #     self._widget = MBDSystemWidget(self._item)

        self._widget._edit(self._item)
        # except:
        #     raise ValueError, "edit() not completed"

    def _delete(self):
        """

        """
        pos = self._item._parent._children.index(self._item)
        self._item._parent.removeChild(pos)
        del(self._item)

    def _properties(self, index):
        """
        Open properties window
        """
        _widget = MBDSystemWidget(self._item, parent=self)
        _widget._edit(self._item)

    def _save_contact_solution(self, index):
        """
        
        """
        self._item.save_solution_data()
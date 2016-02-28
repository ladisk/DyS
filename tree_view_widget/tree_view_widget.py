"""
Created on 27. jan. 2014

@author: lskrinjar
"""
import inspect
import os
import subprocess
from pprint import pprint

import numpy as np
from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from MBD_system import read_and_write
from MBD_system.MBD_system_widget import MBDSystemWidget
from MBD_system.body.body_widget import BodyWidget
from MBD_system.contact.contact_widget import ContactWidget
from MBD_system.force.force_widget import ForceWidget
from MBD_system.item_widget import ItemWidget
from MBD_system.joint.joint_widget import JointWidget
from MBD_system.motion.motion_widget import MotionWidget
from MBD_system.q2R_i import q2R_i
from MBD_system.q2theta_i import q2theta_i
from MBD_system.solution_data.solution_data import SolutionData
from MBD_system.spring.spring_widget import SpringWidget
from analysis.analysis_widget import AnalysisWidget
from tree_model import TreeModel
from tree_view_widget_ui import Ui_Form


class solutionFilenameSignal(QtCore.QObject):
    signal_filename = QtCore.pyqtSignal(str, name='')

class loadSolutionFile(QtCore.QObject):
    signal_loadSolutionFile = QtCore.pyqtSignal()

class CreateAnimationFile(QtCore.QObject):
    signal_createAnimationFile = QtCore.pyqtSignal()

class TreeViewWidget(QWidget):  # QMainWindow#, QAbstractItemView 
    """
    classdocs
    """
    def __init__(self, project_item, MBD_system, parent=None, flags=0):
        """
        Constructor
        """
        super(TreeViewWidget, self).__init__(parent)
        self.ui = Ui_Form()
        self.ui.setupUi(self)

        #   signals
        self.filename_signal = solutionFilenameSignal()
        self.load_loadSolutionFile = loadSolutionFile()
        self.create_animation_file = CreateAnimationFile()

        #   parent
        self._parent = parent

        self.project_item = project_item
        # _model = MBD_system.TreeModel(self.MBD_system)
        _model = TreeModel(self.project_item)

        self.MBD_system = MBD_system

        self.ui.treeView.setModel(_model)
        
        self.ui.treeView.setRootIsDecorated(True)
        self.ui.treeView.setHeaderHidden(True)
        self.ui.treeView.expandAll()

        # we want our listview to have a context self.menu taken from the actions on this widget
        # those actions will be to delete an item :)
        self.ui.treeView.customContextMenuRequested.connect(self.contextMenuEvent)
        QtCore.QObject.connect(self.ui.treeView, QtCore.SIGNAL("clicked (QModelIndex)"),self.row_clicked)

        # self._create_animation_file.connect(self._parent.simulation_control_widget._create_animation_file)
 
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

            if self._item._name == "Motions":
                pprint(vars(self._item))
            
        if self._item._typeInfo != "group" or self._item._typeInfo != "solution":# or self._item._name.lower() != "solution":
            editAction = self.menu.addAction("Edit")
            editAction.triggered.connect(self._edit)

            deleteAction = self.menu.addAction("Delete")
            deleteAction.setEnabled(True)
            deleteAction.triggered.connect(self._delete)

            saveAction = self.menu.addAction("Save")
            saveAction.triggered.connect(self._save)

            self.menu.addSeparator()

        if self._item._typeInfo == "MBDsystem" or self._item._typeInfo == "settings":
            self.menu.addSeparator()
            show_GCS_Action = QtGui.QAction("Show GCS ", self, checkable=True, checked=self._parent.simulation_control_widget.opengl_widget.GCS._visible)
            show_GCS_Action.triggered.connect(self._parent.simulation_control_widget.opengl_widget.GCS._show)
            self.menu.addAction(show_GCS_Action)

            propertiesAction = self.menu.addAction("Properties")
            propertiesAction.triggered.connect(self._properties)
            
            self.menu.addSeparator()

            profileAction = self.menu.addAction("Profile")
            profileAction.triggered.connect(self._parent.simulation_control_widget.profile_functions)
            
            getDOFsAction = self.menu.addAction("DOF")
            getDOFsAction.triggered.connect(self.__getDOF)

            if self._item._typeInfo == "MBDsystem":
                get_q_Action = self.menu.addAction("q")
                get_q_Action.triggered.connect(self._get_q)
                
                save_q_Action = self.menu.addAction("Save q")
                save_q_Action.triggered.connect(self._save_q)
                
                load_q_Action = self.menu.addAction("Load q")
                load_q_Action.triggered.connect(self._load_q)
                
                self.menu.addSeparator()

                evaluate_C_Action = self.menu.addAction("C(q, t)")
                evaluate_C_Action.triggered.connect(self._evaluate_C)

                evaluate_C_q_Action = self.menu.addAction("C_q(q)")
                evaluate_C_q_Action.triggered.connect(self._evaluate_C_q)

                evaluate_C_t_Action = self.menu.addAction("C_t(t)")
                evaluate_C_t_Action.triggered.connect(self._evaluate_C_t)
                
                evaluate_M_Action = self.menu.addAction("M")
                evaluate_M_Action.triggered.connect(self.evaluate_M)

                self.menu.addSeparator()
                list_parameters_Action = self.menu.addAction("List parameters")
                list_parameters_Action.triggered.connect(self._list_parameters)

                analysis_Action = self.menu.addAction("Analysis")
                analysis_Action.triggered.connect(self._analysis)
                analysis_Action.setShortcut(QtGui.QKeySequence("Alt+A"))

                self.menu.addSeparator()

        elif self._item._typeInfo == "solution":
            loadAction = self.menu.addAction("Load from File")
            loadAction.triggered.connect(self._load_solution_file_to_project)

        #   solution data item context menu
        elif self._item._typeInfo == "solutiondata":
            loadAction = self.menu.addAction("Load solution")
            loadAction.triggered.connect(self._load_solution_data)

            propertiesAction = self.menu.addAction("Properties")
            propertiesAction.triggered.connect(self.solution_properties)

            if self._item.loaded:
                saveAnimationFileAction = self.menu.addAction("Save animation to file")
                saveAnimationFileAction.triggered.connect(self._item._create_animation_file)
                saveAnimationFileAction.triggered.connect(self._create_animation_file)

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
            for force, ID in zip(self._item._Fn_list, ["i", "j"]):
                show_force_Action = QtGui.QAction("Show force on Body "+ID, self, checkable=True, checked=force._visible)
                show_force_Action.triggered.connect(force._show_force)
                self.menu.addAction(show_force_Action)

            self.menu.addSeparator()
            for marker, ID, body_id in  zip(self._item.markers, ["i", "j"], self._item.body_id_list):
                if isinstance(body_id, int):
                    name = self.MBD_system.bodies[body_id]._name
                else:
                    name = "ground"
                
                show_marker_Action = QtGui.QAction("Show u_P on body "+ID+" id = "+str(body_id)+" ("+name+")", self, checkable=True, checked=marker._visible)
                show_marker_Action.triggered.connect(marker._show)
                self.menu.addAction(show_marker_Action)

            if self._item._typeInfo == "joint":
                if self._item.joint_type == "prismatic":
                    marker = self._item.markers[-1]
                    show_marker_uQ_Action = QtGui.QAction("Show u_Q on body "+"i", self, checkable=True, checked=marker._visible)
                    show_marker_uQ_Action.triggered.connect(marker._show)
                    self.menu.addAction(show_marker_uQ_Action)

            self.menu.addSeparator()
            if self._item._typeInfo == "joint":
                evaluate_C_Action = QtGui.QAction("Evaluate C(q, t)", self)
                evaluate_C_Action.triggered.connect(self.evaluate_C)
                self.menu.addAction(evaluate_C_Action)

                evaluate_C_q_Action = QtGui.QAction("Evaluate C_q(q)", self)
                evaluate_C_q_Action.triggered.connect(self.evaluate_C_q)
                self.menu.addAction(evaluate_C_q_Action)

                evaluate_Q_d_Action = QtGui.QAction("Evaluate Q_d", self)
                evaluate_Q_d_Action.triggered.connect(self.evaluate_Q_d)
                self.menu.addAction(evaluate_Q_d_Action)

            evaluate_rijP_Action = QtGui.QAction("Evaluate d", self)
            evaluate_rijP_Action.triggered.connect(self.evaluate_d)
            self.menu.addAction(evaluate_rijP_Action)

        elif self._item._typeInfo == "contact":
            for force, ID, body_id in zip(self._item._Fn_list, ["i", "j"], self._item.body_id_list):
                if isinstance(body_id, int):
                    name = self.MBD_system.bodies[body_id]._name
                else:
                    name = "ground"
                show_force_Action = QtGui.QAction("Show Fn on body "+ID+" id = "+str(body_id)+" ("+name+")", self, checkable=True, checked=force._visible)
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

            contactPointAction = QtGui.QAction("Contact point in GCS", self)
            # contactPointAction.setEnabled(self._item._contact_point_found)
            contactPointAction.triggered.connect(self._contact_point)
            self.menu.addAction(contactPointAction)

            self.menu.addSeparator()

            if self._item._type.lower() == "revolute clearance joint":
                for marker, ID in zip(self._item.markers, ["i", "j"]):
                    show_marker_Action = QtGui.QAction("Show Body "+ID+" joint center", self, checkable=True, checked=marker._visible)
                    show_marker_Action.triggered.connect(marker._show)
                    self.menu.addAction(show_marker_Action)

            if self._item._type.lower() == "pin-slot clearance joint linear":
                for marker, ID in zip(self._item.markers, ["i", "j", "j"]):
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
                if self._item.VBO_created:
                    show_body_Action.setEnabled(True)
                else:
                    show_body_Action.setEnabled(False)

                show_LCS_Action = QtGui.QAction("Show LCS", self, checkable=True, checked=self._item.LCS._visible)
                show_LCS_Action.triggered.connect(self._item.LCS._show)
                self.menu.addAction(show_LCS_Action)
                if self._item.LCS._VBO_created:
                    show_LCS_Action.setEnabled(True)
                else:
                    show_LCS_Action.setEnabled(False)

                show_CAD_CS_Action = QtGui.QAction("Show CAD CS", self, checkable=True, checked=self._item.CAD_CS._visible)
                show_CAD_CS_Action.triggered.connect(self._item.CAD_CS._show)
                self.menu.addAction(show_CAD_CS_Action)
                if self._item.CAD_CS._VBO_created:
                    show_CAD_CS_Action.setEnabled(True)
                else:
                    show_CAD_CS_Action.setEnabled(False)

                show_ID_Action = QtGui.QAction("Show ID", self, checkable=True, checked=self._item._idVisible)
                show_ID_Action.triggered.connect(self._item._showID)
                self.menu.addAction(show_ID_Action)

                self.menu.addSeparator()

                getBody_q_Action = self.menu.addAction("Get body q")
                getBody_q_Action.triggered.connect(self._item.get_q)

                getBody_dq_Action = self.menu.addAction("Get body dq")
                getBody_dq_Action.triggered.connect(self._item.get_dq)

                getBody_qdq_Action = self.menu.addAction("Get body [q dq]")
                getBody_qdq_Action.triggered.connect(self._item.get_qdq)
                self.menu.addSeparator()

        elif self._item._typeInfo == "motion":
            evaluate_C_Action = QtGui.QAction("Evaluate C(q, t)", self)
            evaluate_C_Action.triggered.connect(self.evaluate_C)
            self.menu.addAction(evaluate_C_Action)

            evaluate_C_q_Action = QtGui.QAction("Evaluate C_q(q)", self)
            evaluate_C_q_Action.triggered.connect(self.evaluate_C_q)
            self.menu.addAction(evaluate_C_q_Action)

            evaluate_C_t_Action = QtGui.QAction("Evaluate C_t(t)", self)
            evaluate_C_t_Action.triggered.connect(self.evaluate_C_t)
            self.menu.addAction(evaluate_C_t_Action)
        else:
            pass

        self.menu.exec_(event.globalPos())
        self._parent.simulation_control_widget.opengl_widget.updateGL()

    def closeEvent(self, event):
        """

        """
        self._parent.simulation_control_widget.opengl_widget.updateGL()

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
        # data = self._item._load_dat_file()
        self._parent.simulation_control_widget.load_solution_file(filename=self._item._name)

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
        curframe = inspect.currentframe()
        calframe = inspect.getouterframes(curframe, 2)
#         print 'caller name:', calframe
        # pprint(vars(calframe))

#         print "add_solution_data()"
        #   root index
        root_index = self.ui.treeView.rootIndex()

        #   get root node with root index
        self._item = self.ui.treeView.model().getNode(root_index)
        MBD_system_item = self._item._children[0]
        # pprint(vars(MBD_system_item))
        #    solutions (group) item
        # _solution_group_item = MBD_system_item._children[6]


        _solution_group_item = [child for child in MBD_system_item._children if child._name == "Solution"][0]

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
        C_q, C_qT = self._parent.simulation_control_widget.solver.solveODE.ode_fun.getDOF()

    def _get_q(self):
        """

        """
        q = self.MBD_system.get_q()
        print "q ="
        print q
    
    def _load_q(self):
        """
        Function loads vector q from solution file to MBD system object and its bodies
        """
        _path = self.MBD_system.MBD_folder_abs_path_
        _load_file = QtGui.QFileDialog(self)
        _load_file.setDirectory(_path)

        _supported_types = self.__supported_types_solution()

        #   get file path from open file dialog
        file_path = _load_file.getOpenFileName(parent=self._parent, caption=QString("Load q"), directory=QString(_path), filter=QString(_supported_types))
        file_path.replace("/", "\\")

        #   load file data if file is selected
        if file_path:
            #   read data file and store data to solution data object
            solution_data = SolutionData(_file=file_path, MBD_system=self.MBD_system)

            q = self.MBD_system.q0 = solution_data._q_sol_container

            for body in self.MBD_system.bodies:
                #   q to body
                body.R[0:2] = q2R_i(q, body.body_id)
                #   dq to body
                body.theta[2] = q2theta_i(q, body.body_id) 

    def _save_q(self):
        """
        Function saves curent vector q to solution file
        """
        _path = self.MBD_system.MBD_folder_abs_path_
        _save_file = QtGui.QFileDialog(self)
        _save_file.setDirectory(_path)
        
        #    get supported types as string
        _supported_types = self.__supported_types_solution()
        
        #    filepath and filetype
        file_path, filetype = _save_file.getSaveFileNameAndFilter(parent=self._parent, caption=QString("Save q"), directory=QString(_path), filter=QString(_supported_types))
        
        #    get current vector q of MBD system
        q = self.MBD_system.get_q()
        
        #    create solution data object to store vector q to file
        solution_data = SolutionData(name=None, _file=str(file_path), MBD_system=self.MBD_system)
       
        _data = np.array([np.concatenate((np.zeros(5), q), axis=0)])
        solution_data.add_data(_data)
        
        solution_data.write_to_file(_file_abs_path=file_path)
        
    def __supported_types_solution(self):
        """
        
        """
        #   predefine empty string
        _supported_types = ""
        for solution_type, type_name in zip(self.MBD_system._solution_filetypes, self.MBD_system._soltuion_filetype_software):
            _supported_type = type_name + " (*" + solution_type + ")" +";;"

            _supported_types = _supported_types + _supported_type

        #   remove last ;;
        _supported_types = _supported_types[0:-2]
        
        return _supported_types

    def _evaluate_C(self):
        """

        """
        t = 0
        q = self.MBD_system.get_q()
        print "C(q, t):"
        print self._parent.simulation_control_widget.solver.analysis.DAE_fun.evaluate_C(q, t)

    def _evaluate_C_q(self):
        """

        :return:
        """
        q = self.MBD_system.get_q()
        print "C_q(q, t):"
        print self._parent.simulation_control_widget.solver.analysis.DAE_fun.evaluate_C_q(q, 0)

    def _evaluate_C_t(self):
        """

        :return:
        """
        t = 0
        q = self.MBD_system.get_q()
        print "C_t(t):"
        print self._parent.simulation_control_widget.solver.analysis.DAE_fun.evaluate_C_t(q, 0)

    def _edit(self, index):
        """

        """
        self._parentNode = self._item._parent

        if self._item._typeInfo.lower() == "body":
            self._widget = BodyWidget(self._parentNode, parent=self)

        if self._item._typeInfo.lower() == "contact":
            self._widget = ContactWidget(self._parentNode, parent=self)
            pprint(vars(self._item))

        if self._item._typeInfo.lower() == "force":
            self._widget = ForceWidget(self._parentNode, parent=self)

        if self._item._typeInfo.lower() == "joint":
            pprint(vars(self._item))
            self._widget = JointWidget(self._parentNode, parent=self)

        if self._item._typeInfo.lower() == "spring":
            self._widget = SpringWidget(self._parentNode, parent=self)

        if self._item._typeInfo.lower() == "motion":
            self._widget = MotionWidget(self._parentNode, parent=self)

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

    def _save(self):
        """

        :return:
        """
        obj = self._item
        name = self._item._name+".dprj"
        read_and_write.write(obj, name)

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

    def _create_animation_file(self):
        """

        :return:
        """
        print "exe in tree view widget"
        self.create_animation_file.signal_createAnimationFile.emit()
        print "signal emited!!!"

    def evaluate_C(self):
        """

        :return:
        """
        #   get current q vector of MBD system
        q = self.MBD_system.evaluate_q0()

        print self._item.evaluate_C(q)

    def evaluate_C_q(self):
        """

        :return:
        """
        #   get current q vector of MBD system
        q = self.MBD_system.evaluate_q0()

        C_q = self._item.evaluate_C_q(q)
        print "self._item =", self._item
        if self._item._typeInfo == "joint":
            for C_q_i, id in zip(C_q, ["i", "j"]):
                print "C_q %s="%id
                print C_q_i
        else:
            print "C_q ="
            print C_q

    def evaluate_C_t(self):
        """

        :return:
        """
        #   get current q vector of MBD system
        q = self.MBD_system._children[0].evaluate_q0()

        print self._item.evaluate_C_t(q=q, t=0)

    def evaluate_Q_d(self):
        """

        :return:
        """
        #   get current q vector of MBD system
        q = self.MBD_system._children[0].evaluate_q0()

        print self._item.evaluate_Q_d(q)

    def evaluate_d(self):
        """

        :return:
        """
        #   get current q vector of MBD system
        q = self.MBD_system._children[0].evaluate_q0()
        print self._item.evaluate_d(q)
    
    def evaluate_M(self):
        """
        
        :return:
        """
        self._parent.simulation_control_widget.solver.analysis.DAE_fun.evaluate_M()
        M = self._parent.simulation_control_widget.solver.analysis.DAE_fun.M
        print "M ="
        print M

    def _contact_point(self):
        """
        Method gets information of point of contact
        :return:
        """
        if self._parent.simulation_control_widget._status == "animation":
            self._item.get_contact_point(step=self._parent.simulation_control_widget._step)



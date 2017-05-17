"""
Created on 28. jun. 2015

@author: lskrinjar (email: skrinjar.luka@gmail.com)
"""
from pprint import pprint

from PyQt4 import QtCore, QtGui

from MBD_system.array2string import array2string
from MBD_system.contact.contact import Contact
from MBD_system.fix_string import fix_string
from MBD_system.force.force_widget import ForceWidget
from MBD_system.string2array import string2array
from contact_ui import Ui_Form
from MBD_system.checkable_combobox import CheckableComboBox
from MBD_system.contact_model.contact_model import ContactModel


try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s


class ContactWidget(QtGui.QDialog):
    """
    control panel interface
    """
    
    def __init__(self, parent_node=None, parent=None):
        """
        Constructor
        """
        super(ContactWidget, self).__init__(parent)
        self._parent = parent

        self.ui = Ui_Form()
        self.ui.setupUi(self)

        self.ui.contactModelTypecomboBox = CheckableComboBox(self.ui.groupBox_basicInfo)
        self.ui.contactModelTypecomboBox.setMouseTracking(True)
        self.ui.contactModelTypecomboBox.setObjectName(_fromUtf8("contactModelTypecomboBox"))
        self.ui.gridLayout.addWidget(self.ui.contactModelTypecomboBox, 3, 2, 1, 1)

        #   set focus on first widget (index = 0)
        self.ui.tabWidget.setCurrentIndex(0)
        
        self.parent_node = parent_node

        self.ui.save_pushButton.setAutoDefault(True)

        #    set validators to limit user input
        __validator_dbl = QtGui.QDoubleValidator()
        __validator_int = QtGui.QIntValidator()
        __validator_scientific = QtGui.QDoubleValidator()
        __validator_scientific.setNotation(1)

        #    validators
        self.ui.contactID_lineEdit.setValidator(__validator_int)
        self.ui.contact_stiffness_lineEdit.setValidator(__validator_scientific)

        #   body id list at custom context menu
        # self.ui.bodyIDi_lineEdit.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        # self.ui.bodyIDi_lineEdit.customContextMenuRequested.connect(self._contextMenuEvent)
        # self.ui.uPi_lineEdit.setValidator(__validator_dbl)
        # self.ui.Fx_lineEdit.setValidator(__validator_dbl)
        # self.ui.Fy_lineEdit.setValidator(__validator_dbl)
        # self.ui.Mz_lineEdit.setValidator(__validator_dbl)

        #   validators for line edit user input
        self.ui.z_dim_lineEdit.setValidator(__validator_dbl)

        #    signals
        self.ui.cancel_pushButton.clicked.connect(self._cancel)
        self.ui.save_pushButton.clicked.connect(self._save)
        self.ui.save_pushButton.setFocus()
        self.ui.contactTypecomboBox.currentIndexChanged.connect(self._type_changed)

        # pprint(vars(self.group_item))
        self.ui.contactID_lineEdit.setText(str(len(self.parent_node._children)))
        
        #   move to center of screen
        self._move()

        #    show widget if new item is being created
        if parent_node is not None:
            self.show()
        #    signals
#         self.setAttribute(QtCore.Qt.WA_DeleteOnClose)

    def contextMenuEvent(self, event):
        """

        :return:
        """
        menu = QtGui.QMenu(self)

        indexes = self.ui.forces_listWidget.selectedIndexes()

        self._item = self.ui.forces_listWidget.itemFromIndex(indexes[0])

        print self._item
        print self._item.data()
        pprint(vars(self._item))
        # if self._item.typeInfo().lower() == "force":
        #
        viewForceInfo = menu.addAction("View Force Info")
        #
        viewForceInfo.triggered.connect(self._edit)
        menu.addAction(viewForceInfo)
        menu.exec_(event.globalPos())

    def _edit(self, item):
        """

        :return:
        """
        self._widget = ForceWidget(self._item, parent=self)
        self._widget._edit(self._item)

    def _type_changed(self, _index):
        """

        :param type:
        :return:
        """
        _type = self.ui.contactTypecomboBox.itemText(_index)

        if _type == "General":
            self.ui.uPi_lineEdit.setEnabled(False)
            self.ui.uPj_lineEdit.setEnabled(False)
            self.ui.Ri_lineEdit.setEnabled(False)
            self.ui.Rj_lineEdit.setEnabled(False)

        if _type == "Revolute Clearance Joint":
            self.ui.uPi_lineEdit.setEnabled(True)
            self.ui.uPj_lineEdit.setEnabled(True)
            self.ui.Ri_lineEdit.setEnabled(True)
            self.ui.Rj_lineEdit.setEnabled(True)

    def _move(self):
        frameGm = self.frameGeometry()
        screen = QtGui.QApplication.desktop().screenNumber(QtGui.QApplication.desktop().cursor().pos())
        centerPoint = QtGui.QApplication.desktop().screenGeometry(screen).center()
        frameGm.moveCenter(centerPoint)
        self.move(frameGm.topLeft())

    def _cancel(self):
        """
        
        """
        self.close()

#         except:
#             QtGui.QMessageBox.warning(self, "Warning!",
#                 "Input not correct!",
#                 QtGui.QMessageBox.Cancel, QtGui.QMessageBox.NoButton,
#                 QtGui.QMessageBox.NoButton)

    def edit(self, item):
        """

        :return:
        """
        self.item = item

        #   name of contact
        self.ui.name_lineEdit.setText(item._name)

        #   contact id
        self.ui.contactID_lineEdit.setText(str(item.contact_id))

        #   type of contact
        #   add types programamtically
        self.ui.contactTypecomboBox.addItems(self.item._contact_types)
        _index = self.ui.contactTypecomboBox.findText(QtCore.QString(self.item._contact_type))
        self.ui.contactTypecomboBox.setCurrentIndex(_index)

        #   contact model type
        if self.item.contact_model is not None:
            self.ui.contactModelTypecomboBox.addItems(self.item.contact_model._types)
            _index = self.ui.contactModelTypecomboBox.findText(self.item.contact_model._type)
            self.ui.contactModelTypecomboBox.setCurrentIndex(_index)
        else:
            self.ui.contactModelTypecomboBox.addItems(ContactModel._supported_types())
            for i in range(self.ui.contactModelTypecomboBox.count()):
                item = self.ui.contactModelTypecomboBox.model().item(i, 0)
                item.setCheckState(QtCore.Qt.Unchecked)

            for contact_model in self.item.contact_models:
                print "contact_model._type =", contact_model._type
                _index = self.ui.contactModelTypecomboBox.findText(contact_model._type)
                print "_index =" , _index
                if _index != -1:
                    item = self.ui.contactModelTypecomboBox.model().item(_index, 0)
                    item.setCheckState(QtCore.Qt.Checked)

            # for contact_model in ContactModel._get_types():

        #   friction model type
        self.ui.frictionModelTypecomboBox.addItems(self.item.friction_model._types)
        print "GUI friction model"
        print self.item.friction_model._type
        _index = self.ui.frictionModelTypecomboBox.findText(self.item.friction_model._type.lower())
        self.ui.frictionModelTypecomboBox.setCurrentIndex(_index)

        #   coefficient of restitution
        if hasattr(item.contact_model, "_type"):
            if item.contact_model._type.title() == "Hertz":
                self.ui.restitution_coefficient_lineEdit.setEnabled(False)
            else:
                self.ui.restitution_coefficient_lineEdit.setText(str(item.contact_model.c_r))

        #   contact stiffness
        if hasattr(item.contact_model, "K"):
            self.ui.contact_stiffness_lineEdit.setText(str(item.contact_model.K))

        #   list of contact models where n is not used
        # _n_unused = ["Hertz", "Kelvin-Voigt", "Hunt-Crossley", "Herbert-McWhannell", "Flores et al", "Zhiying-Qishao", "Gonthier et al"]

        #   power exponent - n
        # print "contact_model_type =", contact_model_type
        if hasattr(item.contact_model, "n"):
            if item.contact_model != None:
                self.ui.power_exponent_lineEdit.setText(str(item.contact_model.n))
        else:
            self.ui.power_exponent_lineEdit.setEnabled(False)

        #   body id i
        self.ui.bodyIDi_lineEdit.setText(str(item.body_id_i))
        #   body id j
        self.ui.bodyIDj_lineEdit.setText(str(item.body_id_j))
        #   save options
        _index = self.ui.solutionOptions_comboBox.findText(self.item._solution_save_options[0].upper() + self.item._solution_save_options[1:])
        self.ui.solutionOptions_comboBox.setCurrentIndex(_index)
        #   color
        self.ui.color_lineEdit.setText(array2string(item.color))
        #   scale
        self.ui.scale_factor_lineEdit.setText(str(item.scale))

        #   solution options
        self.ui.solutionFiletype_comboBox.addItems(self.item._solution_filetypes)
        _index = self.ui.solutionFiletype_comboBox.findText(self.item._solution_filetype)
        self.ui.solutionFiletype_comboBox.setCurrentIndex(_index)

        #    coefficient of restitution
        # self.ui.restitutionCoef_lineEdit.setText(str(item.contact_model.coef_of_restitution))
        
        #    coefficient of friction static
        # self.ui.frictionCoefStatic_lineEdit.setText(str(item.friction_model.coef_of_friction_static))
        
        #    coefficient of friction dynamic
        # self.ui.frictionCoefDynamic_lineEdit.setText(str(item.friction_model.coef_of_friction_dynamic))

        #   z coordinate of position of contact plane
        # self.ui.z_dim_lineEdit.setText(str(item.z_dim))

        #   TOL of distance
        self.ui.TOL_distance_lineEdit.setText(str(item.distance_TOL))

        #   default dissable all input and enable later only selected one
        self.ui.GeneralContactAddParams.setEnabled(False)
        self.ui.RevoluteClearanceJointAddParams.setEnabled(False)
        self.ui.SphereSphereContactAddParams.setEnabled(False)
        self.ui.PlaneSphereContactAddParams.setEnabled(False)

        #   list of GUI editable properties/parameters for each type of contact
        #   general contact
        self._GeneralContactAddParams = [self.ui.z_dim_lineEdit,
                                        self.ui.skin_thickness_i_lineEdit,
                                        self.ui.skin_thickness_j_lineEdit]

        #   revolute clearance joint contact
        self._RevoluteClearanceJointAddParams = [self.ui.uPi_lineEdit,
                                            self.ui.uPj_lineEdit,
                                            self.ui.Ri_lineEdit,
                                            self.ui.Rj_lineEdit]

        #   contact sphere-sphere
        self._ContactSphereSphere = [self.ui.R0i_lineEdit,
                                    self.ui.R0j_lineEdit]

        #   contact plane-sphere
        self._ContactPlaneSphere = [self.ui.normal_i_lineEdit,
                                    self.ui.R0j_sphere_lineEdit]

        _allAddParams = self._GeneralContactAddParams + self._RevoluteClearanceJointAddParams + self._ContactSphereSphere

        for _param in _allAddParams:
            _param.setEnabled(False)

        if item._type == "general":
            self.ui.GeneralContactAddParams.setEnabled(True)

        elif item._contact_type == "revolute clearance joint":
            self.ui.RevoluteClearanceJointAddParams.setEnabled(True)

            self.ui.uPi_lineEdit.setText(array2string(item.u_iP))
            self.ui.uPj_lineEdit.setText(array2string(item.u_jP))
            self.ui.Ri_lineEdit.setText(str(item.R0_i))
            self.ui.Rj_lineEdit.setText(str(item.R0_j))

        elif item._type.title() == "Contact Sphere-Sphere":
            self.ui.SphereSphereContactAddParams.setEnabled(True)
            self._updateAddParams(self._ContactSphereSphere)

            #   set parameters
            for R0_lineEdit, R_ in zip([self.ui.R0i_lineEdit, self.ui.R0j_lineEdit], item.R0_list):
                R0_lineEdit.setText(str(R_))

        elif item._type.title() == "Contact Plane-Sphere":
            self.ui.PlaneSphereContactAddParams.setEnabled(True)

            self.ui.normal_i_lineEdit.setText(array2string(item.n_i_LCS))
            self.ui.R0j_sphere_lineEdit.setText(str(item.R0_j))

        elif item._contact_type.lower() == "roughness profile":
            for profile, lineEdit in zip(item.roughness_profile_list, [self.ui.roughness_profile_i_lineEdit, self.ui.roughness_profile_j_lineEdit]):
                lineEdit.setText(profile.filename)

        #   add markest to list to display in widget
        # self.ui.markers_listWidget.setFont(self._parent._parent.font)
        # for marker in self.item.markers:
        #     self.ui.markers_listWidget.addItem(marker._name)

        #   add contact forces to list to display in widget
        self.ui.forces_listWidget.setFont(self._parent._parent.font)
        for Fn, Ft in zip(self.item._Fn_list, self.item._Ft_list):
            self.ui.forces_listWidget.addItem(Fn._name + " id = " + str(Fn.force_id))
            self.ui.forces_listWidget.addItem(Ft._name + " id = " + str(Ft.force_id))

        self.show()

    def _save(self, item=None):
        """

        """
        #   name
        _name = self.ui.name_lineEdit.text()

        #   contact type
        _type = self.ui.contactTypecomboBox.currentText()

        #   contact model type
        _type = str(self.ui.contactModelTypecomboBox.currentText()).lower()
        print "_type =", _type
        self.item.contact_model._type = _type
        print "self.item.contact_model._type =", self.item.contact_model._type
        #   distance tolerance
        self.item.distance_TOL = float(self.ui.TOL_distance_lineEdit.text())

#         try:
        _body_id_i = int(self.ui.bodyIDi_lineEdit.text())
        _body_id_j = int(self.ui.bodyIDj_lineEdit.text())

        _coef_of_friction_static = self.ui.frictionCoefStatic_lineEdit.text().toDouble()
        _coef_of_friction_dynamic = self.ui.frictionCoefDynamic_lineEdit.text().toDouble()
        _coef_of_restitution = self.ui.restitution_coefficient_lineEdit.text().toDouble()

        if self.item._type == "Revolute Clearance Joint":
            _uPi = string2array(self.ui.uPi_lineEdit.text())
            _uPj = string2array(self.ui.uPj_lineEdit.text())

            _Ri = float(self.ui.Ri_lineEdit.text())
            _Rj = float(self.ui.Rj_lineEdit.text())

        if self.item._type == "General":
            pass

        if self.item is not None:
            self.item.body_id_i = _body_id_i
            self.item.body_id_j = _body_id_j

            self.item.friction_model.coef_of_friction_static = _coef_of_friction_static
            self.item.friction_model.coef_of_friction_dynamic = _coef_of_friction_dynamic
            self.item.contact_model.coef_of_restitution = _coef_of_restitution

            if self.item._type == "Revolute Clearance Joint":
                self.item.u_iP = _uPi
                self.item.u_jP = _uPj
                self.item.R_i = _Ri
                self.item.R_j = _Rj

        else:
            _item = Contact(body_id_i=_body_id_i, body_id_j=_body_id_j)

            pos = len(self.parent_node._children)
            self.parent_node._parent.forces.append(_item)
            self._parent.ui.treeView.model().insertRow(pos, _item, self.parent_node)

        self.close()

    def _updateAddParams(self, params_list):

        for _param in params_list:
            if _param.isEnabled():
                _param.setEnabled(False)
            else:
                _param.setEnabled(True)
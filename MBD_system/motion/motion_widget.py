"""

created by: lskrinjar
date of creation: 02/02/2016
time of creation: 14:28
"""

from PyQt4 import QtGui

from motion_ui import Ui_Form


class MotionWidget(QtGui.QDialog):
    """
    control panel interface
    """

    def __init__(self, parent_node=None, parent=None):
        """
        Constructor
        """
        super(MotionWidget, self).__init__(parent)
        self._parent = parent


        self.ui = Ui_Form()
        self.ui.setupUi(self)

        self.item = None

        #    signals
        self.ui.buttonBox.button(QtGui.QDialogButtonBox.Cancel).clicked.connect(self._cancel)
        self.ui.buttonBox.button(QtGui.QDialogButtonBox.Save).clicked.connect(self._save)
        self.ui.buttonBox.button(QtGui.QDialogButtonBox.Save).setFocus()


        self._move()

        #   q_lineEdit
        self.q_lineEdit = [self.ui.Rx_lineEdit,
                           self.ui.Ry_lineEdit,
                           self.ui.theta_lineEdit,
                           self.ui.dRx_lineEdit,
                           self.ui.dRy_lineEdit,
                           self.ui.dtheta_lineEdit]

    def _move(self):
        """
        Move widget (dialog window) to screen center
        """
        frameGm = self.frameGeometry()
        screen = QtGui.QApplication.desktop().screenNumber(QtGui.QApplication.desktop().cursor().pos())
        centerPoint = QtGui.QApplication.desktop().screenGeometry(screen).center()
        frameGm.moveCenter(centerPoint)
        self.move(frameGm.topLeft())

    def _cancel(self):
        """

        """
        self.close()

    def _show(self):
        """

        :return:
        """
        self.show()

    def _edit(self, item=None):
        """

        :return:
        """
        self.item = item

        #   name
        self.ui.name_lineEdit.setText(self.item._name)
        #   motion id
        self.ui.motion_id_lineEdit.setText(str(self.item._id))
        #   body id
        self.ui.bodyID_lineEdit.setText(str(self.item.body_id))

        for lineEdit, q in zip(self.q_lineEdit, self.item.q):
            if q is None:
                pass
            else:
                lineEdit.setText(str(q))

        self.show()

    def _save(self):
        """

        :return:
        """
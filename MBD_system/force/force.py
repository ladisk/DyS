'''
Created on 21. feb. 2014

@author: lskrinjar
'''
import numpy as np

import itertools
from OpenGL.GL import *
import numpy as np
import xlsxwriter
from sympy import Symbol
from sympy.parsing.sympy_parser import parse_expr

try:
    from ..MBD_system import *
    from ..cad2cm_lcs import cad2cm_lcs
    from ..force_matrix import Force_Q_e_matrix
except:
    None

try:
    from ..MBD_system_items import ForceItem
except:
    from MBD_system_items import ForceItem

from MBD_system.q2theta_i import q2theta_i
from simulation_control_widget.opengl_widget.marker.marker import Marker
from MBD_system.transform_cs import u_P_cad2cm_lcs

class Force(ForceItem):
    '''
    classdocs
    '''
    __id = itertools.count(-1)

    def __init__(self, body_id, Fx=0, Fy=0, Mz=0, u_iP_f=np.array([0, 0, 0]), force_name=None, data=None, parent=None):
        super(Force, self).__init__(force_name, parent)
        """
        Constructor of force class
        in:
            force_name - string
            body_id - number of body on which the force is applied
            F_x - x component of force
            F_y - y component of force
            M_z - z component of moment
            u_iP_f - vector of acting force in CAD LCS of a body
            
        """
        self._parent = parent
        #    number
        self.force_id = self.__id.next()
        self._comments = ""
        #    name as string
        if force_name is None:
            self._name = "Force_" + str(self.force_id)
        else:
            self._name = force_name
        
        self.body_id = body_id

        self._body_assigned = False

        #   force components Fx, Fy, Mz
        self.Fx = Fx
        self.Fy = Fy
        self.Mz = Mz

        #   function of time
        self._Fx_t = None
        self._Fy_t = None
        self._Mz_t = None
        
        #    position of acting force in CAD LCS of a body
        self.u_iP_f = u_iP_f[0:2]
        self.z_dim_lcs = u_iP_f[2]

        #   force data from file - matrix F(t)
        self.data = data

        #   opengl properties
        self._scale_GL = 1E-1
        self._visible = False

        #   marker parameters
        self.markers = []


        if self._parent is not None:
            if self._parent._name.lower() == "forces":
                #   body handle
                _body = self._parent._parent.bodies[self.body_id]
                #   force position in body LCS (center of mass)
                u_P = u_P_cad2cm_lcs(self.body_id, _body, self.u_iP_f)
                #   added z dim for visualization
                _node = np.array(np.append(u_P, self.z_dim_lcs), dtype='float32')
                #   create marker
                _marker = Marker(_node, visible=True, parent=_body)
                #   add marker to list of markers of a body
                self.markers.append(_marker)
                #   append marker to list
                #   (only one marker as there is only 1 position of applied force on a body)
                _body.markers.append(_marker)
        
        self._solution_containers()


    def load_data_file(self, filename):
        """

        :return:
        """


    def _load_csv(self):
        """

        :return:
        """


    def _show_force(self):
        """

        :return:
        """
        if self._visible:
            self._visible = False
        else:
            self._visible = True


    def _solution_containers(self):
        """
        
        """
        self._u_iP_f_solution_container = []
        self._Fx_solution_container = []
        self._Fy_solution_container = []
        self._Mz_solution_container = []

        self._step_num_solution_container = [0]


    def _evaluate_F(self, t):
        """
        Function evaluates force function if is function of time. String symbol is evaluated using sympy package with
        function evalf()
        :param t:
        :return:
        """
        time = Symbol('time')
        
        
        if isinstance(self.Fx, basestring):
            _exp = parse_expr(self.Fx)
            _Fx = _exp.evalf(subs={time:t})
        else:
            _Fx = self.Fx
            
        
        if isinstance(self.Fy, basestring):
            _exp = parse_expr(self.Fy)
            _Fy = _exp.evalf(subs={time:t})
        else:
            _Fy = self.Fy
        
        
        force = np.array([_Fx, _Fy])
        return force


    def create_force_Q_e_vector(self, t, q):
        """
        calculates force vector from force object properties and vector q (only angle theta is used from q vector)
        in:
            q - MBS system vector
        out:
            vector of force acting on a body (Fx, Fy, Mz)
        """
        #   construct force vector (can be function of time)
        force = self._evaluate_F(t)
        
        
        #   get theta of body from q vector
        _theta = q2theta_i(q, self.body_id)

        #   construct a matrix (function of vector q)
        matrix = Force_Q_e_matrix(u_iP=self.u_iP_f, theta=_theta, body_id=self.body_id).matrix_

        #   form a vector
        vector = np.dot(matrix, force)
        return vector


    def update(self, step, F = np.zeros(3), u_P = np.zeros(2)):
        """

        :param step:
        :param Fx:
        :param Fy:
        :param Mz:
        :return:
        """
        self._step_num_solution_container.append(step)

        self._update_force_vector(F)
        self._update_u_P_vector(u_P)


    def _update_force_vector(self, F):
        """

        :return: None
        """
        self.Fx = F[0]
        self.Fy = F[1]
        self._Fx_solution_container.append(self.Fx)
        self._Fy_solution_container.append(self.Fy)

        
    
    def _update_u_P_vector(self, u_P):
        """
        
        """
        self.u_iP_f = u_P
        self._u_iP_f_solution_container.append(u_P)
        
        
    def _reset_to_initial_state(self):
        """

        :return:
        """
        self.Fx = 0
        self.Fy = 0
        self.Mz = 0

        self.u_iP_f = np.array([0, 0])


    def get_data(self):
        """

        :return:
        """
        print "self._Fx_solution_container ="
        print self._Fx_solution_container
        print len(self._Fx_solution_container)
        print len(self._Fy_solution_container)
        print len(self._Mz_solution_container)


    def save_solution_data(self, filename=None):
        """

        :return:
        """
        self._comments = "Contact force on body id: "+str(self.body_id)
        if filename is None:
            self._solution_filename = "_" + self._name + "_sol"+".xlsx"
        else:
            self._solution_filename = filename

        #   column headers
        _header = ["i-th step", "u_Px", "u_Py", "Fx", "Fy", "Mz"]

        # Create an new Excel file and add a worksheet.
        workbook = xlsxwriter.Workbook(self._solution_filename, {'nan_inf_to_errors': True})
        format_1 = workbook.add_format({'num_format': '#0.000000000000000'})
        format_2 = workbook.add_format({'num_format': '#0.00000000'})
        worksheet = workbook.add_worksheet("data")

        #   write comments
        worksheet.write(0, 0, self._comments)
        #   write header
        worksheet.write_row(1,0,_header)


        worksheet.write_column(2, 0, np.array(self._step_num_solution_container, dtype="float32"))
        worksheet.write_column(2, 1, np.array(self._u_iP_f_solution_container, dtype="float32")[:,0], format_1)
        worksheet.write_column(2, 2, np.array(self._u_iP_f_solution_container, dtype="float32")[:,1], format_1)
        worksheet.write_column(2, 3, np.array(self._Fx_solution_container, dtype="float32"), format_1)
        worksheet.write_column(2, 4, np.array(self._Fy_solution_container, dtype="float32"), format_1)
        worksheet.write_column(2, 5, np.array(self._Mz_solution_container, dtype="float32"), format_1)

        #   column width
        worksheet.set_column('B:E', 20)

        #   close file
        workbook.close()


    def _paint_GL(self, step=None):
        """
        Display force in opengl scene
        :return:
        """
        if self._visible:
            #   paint during integration
            if step is None:
                glColor3f(1.0, 0.0, 0.0)
                glBegin(GL_LINES)
                glVertex3f(self.u_iP_f[0], self.u_iP_f[1], 0.0)
                glVertex3f(self.u_iP_f[0]+self._scale_GL*self.Fx, self.u_iP_f[1]+self._scale_GL*self.Fy, 0)
                glEnd()


            #   paint during animation of results
            else:
                try:
                    glColor3f(1.0, 0.0, 0.0)
                    glBegin(GL_LINES)
                    glVertex3f(self._u_iP_f_solution_container[step][0], self._u_iP_f_solution_container[step][1], 0.)
                    glVertex3f(self._u_iP_f_solution_container[step][0]+self._Fx_solution_container[step]*self._scale_GL, self._u_iP_f_solution_container[step][1]+self._Fy_solution_container[step]*self._scale_GL, 0.)
                    glEnd()
                except:
                    pass
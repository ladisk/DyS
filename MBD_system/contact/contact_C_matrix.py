'''
Created on 20. apr. 2015

@author: lskrinjar
'''
from pprint import pprint
import itertools
import numpy as np
from MBD_system.A import A_matrix
from MBD_system.A_2theta import A_2theta_matrix
from MBD_system.A_theta import A_theta_matrix
from d_ds import d_ds
from d2_ds2 import d2_ds2

class Contact_C_q_matrix(object):
    '''
    classdocs
    '''
    __id = itertools.cycle([0, 1])

    def __init__(self):
        '''
        Constructor
        '''
        self.id = self.__id.next()
        
        if self.id == 0:
            self._sign = +1
        elif self.id == 1:
            self._sign = -1
            
        
        self.matrix = np.zeros([3, 3])
        self.matrix[0:2, 0:2] = self._sign * np.eye(2)
            

class Contact_C_s_matrix(object):
    '''
    classdocs
    '''
    __id = itertools.cycle([0, 1])

    def __init__(self):
        '''
        Constructor
        '''
        self.id = self.__id.next()
        
        self.matrix = np.zeros([3, 1])
        
        if self.id == 0:
            self._sign = +1
        elif self.id == 1:
            self._sign = -1
            

class Contact_C_ss_matrix(object):
    '''
    classdocs
    '''
    __id = itertools.cycle([0, 1])

    def __init__(self, parent=None):
        '''
        Constructor
        '''
        self._parent = parent
        self.id = self.__id.next()
        
        self.matrix = np.zeros([3])

        if self.id == 0:
            self._sign = +1
        elif self.id == 1:
            self._sign = -1
        
    def evaluate(self, theta, u_P, n, t):  #    manjka se tretja komponenta od produkta n.t
        """
        Function calculate values in predefined zero matrix
        """
        dt_ds = d_ds(t)
        self.matrix[0:2] = self._sign * np.dot(A_matrix(theta), dt_ds)

        
        #    body i
        if self.id == 0:
            _theta_i = theta
            _t_i = d2_ds2(t)
            _theta_j = self._parent.body_list[1].theta[2]
            _n_j = self._parent._n_list[1]
            self.matrix[-1] = self._sign * np.dot(np.dot(A_2theta_matrix(_theta_i), _t_i), np.dot(A_matrix(_theta_j), _n_j))
        
        #    body j
        elif self.id == 1:
            _theta_i = self._parent.body_list[0].theta[2]
            _t_i = self._parent._t_list[0]
            _theta_j = theta
            _n_j = d2_ds2(n)
            self.matrix[-1] = self._sign * np.dot(np.dot(A_matrix(_theta_i), _t_i), np.dot(A_2theta_matrix(_theta_j), _n_j))
                    



class Contact_C_qq_matrix(object):
    '''
    classdocs
    '''
    __id = itertools.cycle([0, 1])
    
    def __init__(self, parent=None):
        '''
        Constructor
        '''
        self._parent = parent
        self.id = self.__id.next()
        
        self.matrix = np.zeros([3, 3])

        if self.id == 0:
            self._sign = +1
        elif self.id == 1:
            self._sign = -1
            
            
    def evaluate(self, theta, u_P, n, t):
        """
        
        """
        self.matrix[0:2, -1] = self._sign * np.dot(A_2theta_matrix(theta), u_P)
        
        #    body i
        if self.id == 0:
            _theta_i = theta
            #    second partial derivative of t_i with respect to s_i
            _d2t_i_ds2_i = d2_ds2(t)
            _theta_j = self._parent.body_list[1].theta[2]
            _n_j = self._parent._n_list[1]
            self.matrix[-1, -1] = self._sign * np.dot(np.dot(A_2theta_matrix(_theta_i), _d2t_i_ds2_i), np.dot(A_matrix(_theta_j), _n_j))
        
        #    body j
        elif self.id == 1:
            _theta_i = self._parent.body_list[0].theta[2]
            _t_i = self._parent._t_list[0]
            _theta_j = theta
            #    second partial derivative of n_j with respect to s_j
            _d2n_j_ds2_j = d2_ds2(n)
            self.matrix[-1, -1] = self._sign * np.dot(np.dot(A_matrix(_theta_i), _t_i), np.dot(A_2theta_matrix(_theta_j), _d2n_j_ds2_j))
                    

class Contact_C_qs_matrix(object):
    '''
    classdocs
    '''
    __id = itertools.cycle([0, 1])
    
    def __init__(self, parent=None):
        self._parent = parent
        self.id = self.__id.next()
        
        self.matrix = np.zeros([3, 3])

        if self.id == 0:
            self._sign = +1
        elif self.id == 1:
            self._sign = -1
    
    
    def evaluate(self, theta, u_P, n, t):
        """
        
        """
        self.matrix[-1, 0:2] = self._sign * np.dot(A_theta_matrix(theta), t)

        #    body i
        if self.id == 0:
            _theta_i = theta
            #    second partial derivative of t_i with respect to s_i
            _dt_i_ds_i = d_ds(t)
            _theta_j = self._parent.body_list[1].theta[2]
            _n_j = self._parent._n_list[1]
            self.matrix[-1, -1] = self._sign * np.dot(np.dot(A_2theta_matrix(_theta_i), _dt_i_ds_i), np.dot(A_matrix(_theta_j), _n_j))
        
        #    body j
        elif self.id == 1:
            _theta_i = self._parent.body_list[0].theta[2]
            _t_i = self._parent._t_list[0]
            _theta_j = theta
            #    second partial derivative of n_j with respect to s_j
            _dn_j_ds_j = d_ds(n)
            self.matrix[-1, -1] = self._sign * np.dot(np.dot(A_matrix(_theta_i), _t_i), np.dot(A_2theta_matrix(_theta_j), _dn_j_ds_j))
                   

class Contact_C_matrix(object):
    '''
    classdocs
    '''
    __id = itertools.cycle([0, 1])

    def __init__(self, parent=None):
        '''
        Constructor
        '''
        self.C_qq_matrix = None
        self.C_ss_matrix = None
        self.C_qs_matrix = None
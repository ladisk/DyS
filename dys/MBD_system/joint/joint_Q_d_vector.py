'''
Created on 17. apr. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
'''

import itertools

import numpy as np


class Joint_Q_d_vector(object):
    '''
    classdocs
    '''
    __id = itertools.cycle([0, 1])
    
    
    def __init__(self, joint_type_ = [], body_connected_to_ground = False):
        """
        Constructor of joint Q_d vector class for each joint the constructor is called twice
        in:
        """
            
        if body_connected_to_ground:
            self.id = 0
        else:
            self.id = self.__id.next()

        
        if joint_type_ == "fixed":
            self.joint_Q_d_vector = np.zeros(3)
        elif joint_type_ == "revolute":
            self.joint_Q_d_vector = np.zeros(2)
        elif joint_type_ == "prismatic":
            self.joint_Q_d_vector = np.zeros(2)
        else:
            raise ValueError, "Joint Q_d vector not constructed!"
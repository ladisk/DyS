"""
Created on 8. jun. 2016

@author: skrinjar.luka@gmail.com

"""
import numpy as np
def find_nearest(array, value):
    idx = (np.abs(array-value)).argmin()
    return idx
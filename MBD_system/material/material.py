"""
Created on 12. feb. 2016

@author: luka.skrinjar
"""

import itertools


class Material(object):
    """
    classdocs
    """
    __id = itertools.count(0)

    def __init__(self, name=None, _file=None, MBD_system=None, parent=None):
        """
        Constructor of solution data class
        """
        super(Material, self).__init__(name, parent)
        
        #    name of material
        self._name = name
        
        #    density
        self.density = density
        
        #    module of elasticity
        self.module_of_elasticiy = module_of_elasticiy
        
        #    poison ratio
        self.poison_ratio = poison_ratio
'''
Created on 9. jul. 2014

@author: lskrinjar
'''
import itertools
import xlsxwriter
import numpy as np


try:
    from ..MBD_system import *
    from ..cad2cm_lcs import cad2cm_lcs
    from ..force_matrix import Force_Q_e_matrix
except:
    None

try:
    from ..MBD_system_items import SolutionDataItem
except:
    from MBD_system_items import SolutionDataItem


class SolutionData(SolutionDataItem):
    '''
    classdocs
    '''
    __id = itertools.count(0)

    def __init__(self, name, file=None, parent=None):
        super(SolutionData, self).__init__(name, parent)
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

        #   abs path to file
        self._file = file
        
        #    set file type: .dat, .xlsx, .csv
        self._filetype = ".dat" 


    def set_filetype(self, filetype):
        self._filetype = filetype


    def load_file(self):
        """

        """
        if self._filetype == ".dat":
            return self._load_dat_file()
        elif self._filetype == ".xlsx":
            self._load_excel_file()
        elif self._filetype == ".csv":
            self._load_csv_file()


    def _load_dat_file(self):
        """

        """
        data = np.loadtxt(str(self._file), skiprows=2)
        return data

    def _load_excel_file(self):
        print "Under Construction: ", os.path.realpath(__file__)


    def _load_csv_file(self):
        print "Under Construction: ", os.path.realpath(__file__)
    
    
    def add_data(self):
        """
        
        """
        
    
    def save_to_file(self):
        """
        Function saves data to file
        """
        
        
    def _save_to_dat_file(self):
        """
        Function saves data to .dat filetype
        """
        
    def _save_to_excel_file(self):
        """
        Function saves data to .xlsx filetype
        """
        
    def _save_to_csv_file(self):
        """
        Function saves data to .csv filetype
        """

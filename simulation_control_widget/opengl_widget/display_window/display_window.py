'''
Created on 22. nov. 2013

@author: lskrinjar
'''

class Window(object):
    '''
    classdocs
    '''


    def __init__(self, width, height, left_GL, right_GL, bottom_GL, top_GL, near_GL, far_GL, zoom, scale_factor_Translation, scale_factor_Rotation):
        '''
        Constructor
        '''
        self.width = width
        self.height = height
        self.aspect = float(self.width) / self.height
        
#         print "------------------------------"
#         print "class Window props"
#         print "self.width =", self.width
#         print "self.height =", self.height
#         print "self.aspect =", self.aspect
#         print "------------------------------"
        
        #    opengl properties for glOrtho()
        self.left_GL = left_GL
        self.right_GL = right_GL
        self.bottom_GL = bottom_GL 
        self.top_GL = top_GL
        self.near_GL = near_GL 
        self.far_GL = far_GL
        
        #   scaling factors
        self.scale_factor = zoom
        self.scale_factor_Translation = scale_factor_Translation * self.scale_factor
        self.scale_factor_Rotation = scale_factor_Translation * self.scale_factor
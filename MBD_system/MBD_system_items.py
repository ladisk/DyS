'''
Created on 2. jan. 2015

@author: lskrinjar
'''
import itertools


class TreeItem(object):
    """
    
    """
    __id = itertools.count(0)

    def __init__(self, name, parent=None):
        self._name = name
        self._children = []
        self._parent = parent
        self._type = "node"
        self._comments = ""
        self._id = self.__id.next()
        self._constructed = True
        self._selected = False
        self._save_options = None
        if parent is not None:
            parent.addChild(self)
        
#             if name in parent.list_of_object_groups:
#                 self._type = "group"

    def typeInfo(self):
        return self._type

    def addChild(self, child):
        self._children.append(child)

    def insertChild(self, position, child):
        
        if position < 0 or position > len(self._children):
            return False
        
        self._children.insert(position, child)
        child._parent = self
        return True

    def removeChild(self, position):
        
        if position < 0 or position > len(self._children):
            return False
        
        child = self._children.pop(position)
        child._parent = None

        return True

    def name(self):
        return self._name

    def setName(self, name):
        self._name = name

    def child(self, row):
        return self._children[row]
    
    def childCount(self):
        return len(self._children)

    def parent(self):
        return self._parent
    
    def setParent(self, parent):
        self._parent = parent
    
    def row(self):
        if self._parent is not None:
            return self._parent._children.index(self)

class ProjectItem(TreeItem):
    def __init__(self, name, parent=None):
        super(ProjectItem, self).__init__(name, parent)
        self._type = "projectNode"

class AnalysisTreeItem(TreeItem):
    def __init__(self, name, item=None, mean_value=None, parent=None):
        super(AnalysisTreeItem, self).__init__(name, parent)
        self._typeInfo = None
        self._type = "analysisNode"
        
        self._name = name
        self._item = item
        
        #   attributes to build tree item object of each parameter of a MBD system for monte carlo simulation
        self._mean_value = mean_value
        self._min_value = None
        self._max_value = None
        self._random_value = None

    def mean_value(self):
        return self._mean_value

    def min_value(self):
        return self._min_value

    def max_value(self):
        return self._max_value

class MBDsystemItem(TreeItem):
    def __init__(self, name, parent=None):
        super(MBDsystemItem, self).__init__(name, parent)
        self._typeInfo = "MBDsystem"
        self._type = "projectNode"
    def typeInfo(self):
        return "MBDsystem"

class SettingsGroupItem(TreeItem):
    def __init__(self, name, parent=None):
        super(SettingsGroupItem, self).__init__(name, parent)
        self._typeInfo = "settings"
    
    def typeInfo(self):
        return self._typeInfo
    
class SolutionGroupItem(TreeItem):
    def __init__(self, name, parent=None):
        super(SolutionGroupItem, self).__init__(name, parent)
        self._typeInfo = "solution"
    
    def typeInfo(self):
        return self._typeInfo

class SolutionDataItem(TreeItem):
    def __init__(self, name, parent=None):
        super(SolutionDataItem, self).__init__(name, parent)
        self._typeInfo = "solutiondata"
    
    def typeInfo(self):
        return self._typeInfo

class GroupItem(TreeItem):
    def __init__(self, name, parent=None):
        super(GroupItem, self).__init__(name, parent)
        self._typeInfo = "group"
    
    def typeInfo(self):
        return self._typeInfo
        
class BodyItem(TreeItem):
    def __init__(self, name, parent=None):
        super(BodyItem, self).__init__(name, parent)
        self._typeInfo = "body"
    
    def typeInfo(self):
        return self._typeInfo

class ContactItem(TreeItem):
    def __init__(self, name, parent=None):
        super(ContactItem, self).__init__(name, parent)
        self._typeInfo = "contact"
        
    def typeInfo(self):
        return self._typeInfo
    
class ForceItem(TreeItem):
    def __init__(self, name, parent=None):
        super(ForceItem, self).__init__(name, parent)
        self._typeInfo = "force"
    
    def typeInfo(self):
        return self._typeInfo
    
class JointItem(TreeItem):
    def __init__(self, name, parent=None):
        super(JointItem, self).__init__(name, parent)
        self._typeInfo = "joint"
        
    def typeInfo(self):
        return self._typeInfo

class SpringItem(TreeItem):
    def __init__(self, name, parent=None):
        super(SpringItem, self).__init__(name, parent)
        self._typeInfo = "spring"
        
    def typeInfo(self):
        return self._typeInfo
    
class AABBItem(TreeItem):
    def __init__(self, name, parent=None):
        super(AABBItem, self).__init__(name, parent)
        self._typeInfo = "aabb"
    
    def typeInfo(self):
        return self._typeInfo
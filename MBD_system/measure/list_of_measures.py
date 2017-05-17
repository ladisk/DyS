"""

created by: lskrinjar
date of creation: 23/04/2016
time of creation: 11:06
"""
import os
from pprint import pprint
import re
import numpy as np


def create_list(filename, parent=None):
    """

    :param motions_lst: parent of motion objects
    :param filename:    absolute path to file to read motions
    :param parent:      parent to motion object
    :return:
    """
    from MBD_system.measure.measure import Measure
    from MBD_system.measure.measure_body import MeasureBody
    from MBD_system.measure.measure_contact import MeasureContact
    from MBD_system.measure.measure_force import MeasureForce
    from MBD_system.measure.measure_MBD_system import MeasureMBDsystem
    measures = []
    if os.path.isfile(filename):
        with open(filename, "r") as _file:

            #   define variables
            _type = None
            _body_id = None
            _y_variable = None
            _x_variable = "time"
            _id = None
            _MBD_system = None
            _node_id = None

            for line in _file:
                line = line.strip()

                if line.startswith("MEASURES FILE-START") or line.startswith("#"):
                    pass
                
                elif  line.startswith("MEASURES FILE-END"):
                    _file.close()
                    break
                    
                    
                elif line.startswith("type"):
                    match_equal = re.search(r"=", line)
                    _type = line[match_equal.start() + 1::].strip()

                elif line.startswith("body_id"):
                    match_equal = re.search("=", line)
                    _body_id = np.int(line[match_equal.start() +1::])

                elif line.startswith("x_variable"):
                    match_equal = re.search(r"=", line)
                    _x_variable = line[match_equal.start() + 1::].strip()

                elif line.startswith("y_variable"):
                    match_equal = re.search(r"=", line)
                    _y_variable = line[match_equal.start() + 1::].strip()
                
                elif line.startswith("contact_id"):
                    match_equal = re.search(r"=", line)
                    _id = int(line[match_equal.start() + 1::])
                
                elif line.startswith("force_id"):
                    match_equal = re.search(r"=", line)
                    _id = int(line[match_equal.start() + 1::])

                elif line.startswith("node_id"):
                    match_equal = re.search(r"=", line)
                    _node_id = int(line[match_equal.start() + 1::])
                else:
                    pass

                # print "_type =", _type
                # print "_body_id =", _body_id
                # print "__measure_type =", _measure_type
                # print _type == "body" and _body_id!=None and _measure_type!=None
                if _type == "body" and _body_id is not None and _y_variable is not None:
                    measure = MeasureBody(_y_variable, _body_id, x_variable=_x_variable, node_id=_node_id, parent=parent)

                    measures.append(measure)

                    #   reset values
                    _type = None
                    _body_id = None
                    _y_variable = None
                    _x_variable = "time"
                    _node_id = None
                
                if _type == "contact" and _id is not None and _y_variable is not None:
                    measure = MeasureContact(_y_variable, _id, x_variable=_x_variable, parent=parent)

                    measures.append(measure)
                    
                    #   reset values
                    _type = None
                    _id = None
                    _y_variable = None
                    _x_variable = "time"
                
                if _type == "force" and _id is not None and _y_variable is not None:
                    measure = MeasureForce(_y_variable, _id, x_variable=_x_variable, parent=parent)

                    measures.append(measure)

                    #   reset values
                    _type = None
                    _id = None
                    _y_variable = None
                    _x_variable = "time"

                if _type == "MBDsystem" and _y_variable is not None:
                    measure = MeasureMBDsystem(_y_variable, parent=parent)

                    measures.append(measure)

                    #   reset values
                    _type = None
                    _y_variable = None
                    _x_variable = "time"

    return measures

if __name__ == "__main__":
    measures = create_list("motions.dat")
    print "list of motions:", len(measures)
    for i, object in enumerate(measures):
        print "------------------------"
        print i, object
        pprint(vars(object))
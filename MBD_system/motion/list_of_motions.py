"""

created by: lskrinjar
date of creation: 01/02/2016
time of creation: 10:22
"""
import os
import re
from pprint import pprint

import numpy as np

from MBD_system.motion.motion import Motion


def create_list(filename, parent=None):
    """

    :param motions_lst: parent of motion objects
    :param filename:    absolute path to file to read motions
    :param parent:      parent to motion object
    :return:
    """
    motions = []
    if os.path.isfile(filename):
        with open(filename, "r") as file_:

            #    motion - predefine parameters
            __motion_name = _motion_name = None
            _body_id = None
            _Rx = _Ry = _theta = _dRx = _dRy = _dtheta = None

            for i, line in enumerate(file_):
                line = line.strip()

                #   remove comment - after # sign
                if "#" in line:
                    line = line[:line.index("#")].strip()

                if line.startswith("MOTIONS FILE-START") or line.startswith("#"):
                    pass

                    # print "_read =", _read
                    # if (_motion_name != None) and (_body_id != None) and _read:
                    #
                    #     #    create motion object here
                    #     motion_ = Motion(_motion_name, _body_id, q0=_q0, q1=_q1, q2=_q2, q3=_q3, q4=_q4, q5=_q5, parent=parent)
                    #     motions.append(motion_)
                    #
                    #     #    motion
                    #     motion_name_ = None
                    #     body_id_ = None
                    #     q0_ = q1_ = q2_ = q3_ = q4_ = q5_ = None
                    #     _read = False


                #   motion name
                elif line.startswith("motion_name") or line.startswith("MOTIONS FILE-END"):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)

                    if _body_id is None:
                        __motion_name = line[match_equal.start() + 1::]

                    else:
                        __motion_name = _motion_name

                    if line.startswith("MOTIONS FILE-END"):
                        pass
                    else:
                        _motion_name = line[match_equal.start() + 1::]

                    if _body_id is not None:
                        motion_ = Motion(__motion_name, _body_id, Rx=_Rx, Ry=_Ry, theta=_theta, dRx=_dRx, dRy=_dRy, dtheta=_dtheta, parent=parent)
                        motions.append(motion_)

                        _body_id = None
                        _Rx = _Ry = _theta = _dRx = _dRy = _dtheta = None

                #   body id
                elif line.startswith("body_id"):
                    match_equal = re.search("=", line)
                    match_newline = re.search("\n", line)
                    # print "testing =", "\n" in line
                    # print "match_newline =", match_newline
                    # print "line:", line
                    # print "start =", line[match_equal.start()+1]
                    # print "end =", line[match_newline.start()]
                    _body_id = np.int(line[match_equal.start() +1::])

                #   Rx
                elif line.startswith("Rx"):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    _Rx = line[match_equal.start() + 1::]

                #   Ry
                elif line.startswith("Ry"):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    _Ry = line[match_equal.start() + 1::]

                #   theta
                elif line.startswith("theta"):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    _theta = line[match_equal.start() + 1::]

                #   dRx
                elif line.startswith("dRx"):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    _dRx = line[match_equal.start() + 1::]

                #   dRy
                elif line.startswith("dRy"):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    _dRy = line[match_equal.start() + 1::]

                #   theta
                elif line.startswith("dtheta"):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    _theta = line[match_equal.start() + 1::]

                else:
                    pass

                if line.startswith("MOTIONS FILE-END"):
                    file_.close()
                    break

    return motions

if __name__ == "__main__":
    motions = create_list("motions.dat")
    print "list of motions:", len(motions)
    for i, object in enumerate(motions):
        print "------------------------"
        print i, object
        pprint(vars(object))

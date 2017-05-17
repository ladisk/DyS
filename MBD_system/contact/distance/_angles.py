"""

created by: lskrinjar
date of creation: 27/03/2017
time of creation: 23:00
"""
import numpy as np


riP = np.array([0, 0], dtype=float)
riR = np.array([1, 0], dtype=float)
rjP = np.array([1.15, -0.1], dtype=float)

riPR = riR - riP
# print "riPR =", riPR, "length =", np.linalg.norm(riPR)

t = riPR / np.linalg.norm(riPR)
print "t =", t

angles = [0., 0.]

for i, sign, riL in zip(range(0, len(angles)), [-1, +1], [riP, riR]):
    print "----------------------------------"
    r_jPiL = rjP - riL
    print "r_jPiL =", r_jPiL
    t_i = sign * t
    print "t_i =", t_i
    # print "r_jPiL =", r_jPiL
    # print "fi_t =", np.rad2deg(np.arctan2(t[1], t[0]))
    # print "fi_r =", np.rad2deg(np.arctan2(r_jPiL[1], r_jPiL[0]))
    # print "test =", np.rad2deg(np.arccos(np.dot(r_jPiL, t)/np.linalg.norm(r_jPiL)))
    # print "test =", np.rad2deg(np.arctan2(r_jPiL[1], r_jPiL[0]) - np.arctan2(t[1], t[0]))
    # angles[i] = np.arctan2(t[1] - r_jPiL[1], t[0] - r_jPiL[0])
    print "angle =", np.rad2deg(np.arctan2(t_i[1] - r_jPiL[1], t_i[0] - r_jPiL[0]))

# print "angles ="
# print np.rad2deg(angles)


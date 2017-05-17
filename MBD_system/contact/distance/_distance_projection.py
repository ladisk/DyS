# coding=utf-8
"""
Created on 23. feb 2017

@author: luka.skrinjar
"""
import numpy as np
from matplotlib import pyplot as plt


# riP = np.array([0, 0], dtype=float)
# riR = np.array([1, 0], dtype=float)
# rjP = np.array([0.5, .11], dtype=float)

riP = np.array([-0.000833658094, -0.000216038301], dtype=float)
riR = np.array([8.482413250022e-05,   2.663517625479e-05], dtype=float)
rjP = np.array([9.999999747379e-06,  -4.199999966659e-05], dtype=float)

# riP = np.array([-8, -2], dtype=float)
# riR = np.array([1, .2], dtype=float)
# rjP = np.array([.1,  -2.1], dtype=float)

R0j = 0.75E-3
h0i = 1.6E-3
R0i = h0i / 2.
riPR = riR - riP
print "riPR =", riPR, "length =", np.linalg.norm(riPR)

t = riPR / np.linalg.norm(riPR)
# print "t =", t

riPjP = rjP - riP
print "riPjP =", riPjP

# rCP = riP + np.dot(riPjP, t) * t
#
# print "rCP =", rCP

#   projection
riPjP_proj = np.dot(riPjP, t) * t
print "riPjP_proj =", riPjP_proj

#   distance vector
d = -riPjP_proj + riPjP
print "distance_vector =", d

#   direction
n = np.array([-t[1], t[0]])
print "n =", n
if np.dot(d, n) > 0.:
    print "CONTACT"
else:
    print "NO CONTACT"

#   normal projection
rijPP_n = np.dot(riPjP, n) * n
print "rijPP_n =", rijPP_n

#   plot
fig = plt.figure(num=1, figsize=(6, 5), dpi=72, facecolor='w', edgecolor='g')
ax = plt.subplot(111, aspect="equal")

#   pin - body j
# plt.plot(rjP[0], rjP[1], color='r')
# circle_j = plt.Circle((rjP[0], rjP[1]), R0j, color='r', fill=False)
# ax.add_artist(circle_j)

#   slot - body i
# circle_iP = plt.Circle((riP[0], riP[1]), R0i, color='b', fill=False)
# ax.add_artist(circle_iP)
# circle_iR = plt.Circle((riR[0], riR[1]), R0i, color='b', fill=False)
# ax.add_artist(circle_iR)

#   plot points
for r in [riP, riR, rjP]:
    plt.plot(r[0], r[1], marker="o", color="black")

#   tangent - center line
plt.plot([riP[0], riR[0]], [riP[1], riR[1]], color="black")
#   rijP
plt.plot([riP[0], riP[0] + riPjP[0]], [riP[1], riP[1] + riPjP[1]], color="red")
#   rijP projection
plt.plot([riP[0], riP[0] + riPjP_proj[0]], [riP[1], riP[1] + riPjP_proj[1]], color="green", linewidth=4)
#   normal distance
plt.plot([riP[0] + riPjP_proj[0], riP[0] + riPjP_proj[0] + rijPP_n[0]], [riP[1] + riPjP_proj[1], riP[1] + riPjP_proj[1] + rijPP_n[1]], color="red", linewidth=4, linestyle="--")
plt.show()


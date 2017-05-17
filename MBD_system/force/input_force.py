# coding=utf-8
import sys
import numpy as np


F0 = -1.
Fx = 0.
if t <= 1.:
	Fy = -1.0 * t
else:
	Fy = -1.0
import numpy as np



a = np.arange(0, 1000, 1)
b = np.zeros(len(a))

for i in range(0, len(a)):

    print "a_i =", a[i]
    # print('i=',i)
    if i < 10:
        b[i] = a[i]

print b


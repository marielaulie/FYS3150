"""
Making a plotting program for project 4.
The program will read the output file and use these values to make a plot.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

n = []
E = []


folder = "../FYS3150/Project4/build-Project4-Debug/output/T1/"
text = sys.argv[1]
F = np.genfromtxt(folder+text, delimiter="")
n = F[:,0]
E = F[:,1]
"""
f = open(folder + text)
lines = f.read().split()
x = lines
n.append(float(x[0]))
E.append(float(x[1]))

f.close()
"""
"""
text = sys.argv[2]
f = open(folder + text)
lines = f.read().split()
y = lines
n.append(float(y[0]))
E.append(float(y[1]))


f.close()

text = sys.argv[3]
f = open(folder + text)
lines = f.read().split()
z = lines
n.append(float(z[0]))
E.append(float(z[1]))


f.close()

text = sys.argv[4]
f = open(folder + text)
lines = f.read().split()
k = lines
n.append(float(k[0]))
E.append(float(k[1]))


f.close()

text = sys.argv[5]
f = open(folder + text)
lines = f.read().split()
l = lines
n.append(float(l[0]))
E.append(float(l[1]))
f.close()

text = sys.argv[6]
f = open(folder + text)
lines = f.read().split()
m = lines
n.append(float(m[0]))
E.append(float(m[1]))
f.close()
"""
"""text = sys.argv[7]
f = open(folder + text)
lines = f.read().split()
q = lines
print(q, text)
n.append(float(q[0]))
E.append(float(q[1]))
f.close()
"""

print(n)
print(E)

plt.title("T= and ordered initial spin matrix")
#plt.axis([0,100000,-0.8,1])
plt.xlabel("MonteCarlo cycles n")
plt.ylabel("Energy E")
plt.plot(n,E)
plt.legend()
plt.show()

#plt.hist(E, bins="auto", normed =True, label="")
#plt.show()
#semilogx
#plt.hist

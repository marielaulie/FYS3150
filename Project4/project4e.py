

import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import interpolate

T_40 = []
meanE_40 = []
meanM_40 = []
Cv_40 = []
chi_40 = []

T_60 = []
meanE_60 = []
meanM_60 = []
Cv_60 = []
chi_60 = []

T_80 = []
meanE_80 = []
meanM_80 = []
Cv_80 = []
chi_80 = []

T_100 = []
meanE_100 = []
meanM_100 = []
Cv_100 = []
chi_100 = []


text40 = sys.argv[1]

f = open(text40)
F40 = np.genfromtxt(text40, delimiter="")
T_40 = F40[:,0]
meanE_40 = F40[:,1]
meanM_40 = F40[:,2]
Cv_40 = F40[:,3]
chi_40 = F40[:,4]
f.close()

text60 = sys.argv[2]

f = open(text60)
F60 = np.genfromtxt(text60, delimiter="")
T_60 = F60[:,0]
meanE_60 = F60[:,1]
meanM_60 = F60[:,2]
Cv_60 = F60[:,3]
chi_60 = F60[:,4]
f.close()

text80 = sys.argv[3]

f = open(text80)

F80 = np.genfromtxt(text80, delimiter="")
T_80 = F80[:,0]
meanE_80 = F80[:,1]
meanM_80 = F80[:,2]
Cv_80 = F80[:,3]
chi_80 = F80[:,4]
f.close()

text100 = sys.argv[4]

f = open(text100)

F100 = np.genfromtxt(text100, delimiter="")
T_100 = F100[:,0]
meanE_100 = F100[:,1]
meanM_100 = F100[:,2]
Cv_100 = F100[:,3]
chi_100 = F100[:,4]
f.close()

n = len(T_40)
k = n
M40 = np.amax(Cv_40)
M60 = np.amax(Cv_60)
M80 = np.amax(Cv_80)
M100 = np.amax(Cv_100)

for i in range (0, k):
    if (Cv_40[i] == M40):
        maxT40 = T_40[i]

for i in range (0, k):
    if (Cv_60[i] == M60):
        maxT60 = T_60[i]

for i in range (0, k):
    if (Cv_80[i] == M80):
        maxT80 = T_80[i]

for i in range (0, k):
    if (Cv_100[i] == M100):
        maxT100 = T_100[i]






print("For L=40 the maximum value of heat capacity is %s" %M40, "and the critical temperature is %.3f" %maxT40)
print("For L=60 the maximum value of heat capacity is %s" %M60, "and the critical temperature is %.3f" %maxT60)
print("For L=80 the maximum value of heat capacity is %s" %M80, "and the critical temperature is %.3f" %maxT80)
print("For L=100 the maximum value of heat capacity is %s" %M100, "and the critical temperature is %.3f" %maxT100)


plt.plot(T_40, Cv_40, "r", label='L=40', linewidth=2)
plt.plot(T_40, Cv_60, "g", label='L=60', linewidth=2)
plt.plot(T_40, Cv_80, "b", label='L=80', linewidth=2)
plt.plot(T_40, Cv_100, "y", label='L=100',linewidth=2)

plt.plot(maxT40, M40, marker='o', markersize=10, label='max L=40')
plt.plot(maxT60, M60, marker='o', markersize=10, label='max L=60')
plt.plot(maxT80, M80, marker='o', markersize=10, label='max L=80')
plt.plot(maxT100, M100, marker='o', markersize=10, label='max L=100')

plt.xlabel("Temperature [kT/J]")
plt.ylabel("Heat capacity")
plt.legend(loc = "upper left")
plt.savefig('Cv_L')
plt.show()







"""

#plt.subplot(3,2,1)
plt.plot(T, meanE, 'g')
plt.title("Accepted configurations as a function of temperature ", fontweight='bold')
plt.suptitle("MCcycles = 1.000.000, L=20, ordered initial matrix")
plt.xlabel("Temperature [kT/J]")
plt.ylabel("Accepted configurations per spin")
"""
"""
plt.subplot(3,2,2)
plt.plot(T, meanM, 'r')
plt.title('Mean absolute magnetization')
plt.xlabel('Temperatur [kT/J]')
plt.ylabel('|M(T)|')

plt.subplot(3,2,5)
plt.plot(T, Cv, 'y')
plt.title('Spesific heat')
plt.xlabel('Temperatur [kT/J]')
plt.ylabel('Cv(T)')

plt.subplot(3,2,6)
plt.plot(T, chi, 'b')
plt.title('Susceptibility ')
plt.xlabel('Temperatur [kT/J]')
plt.ylabel('x(T)')
"""
"""
plt.savefig('acc_T_O')
plt.show()

"""

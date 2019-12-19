import numpy as np, matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# plt.style.use('seaborn-darkgrid')


n = 10
h = (1.0/(n+1))
dt = 0.25*h*h
alpha = dt/h/h










for i in range(0,n+1):
    for j in range(0,n+1):
        utemp = 1.0

x = np.linspace(0, 1, n+1)
y = np.linspace(0, 1, n+1)

u = [[ 0,0,0,0,0,0,0,0,0,0,0],
[0  , 0.3164 ,  0.5280 ,  0.6515  , 0.7124 ,  0.7303 ,  0.7124 ,  0.6515  , 0.5280 ,0.3164,0],
[ 0 ,  0.6793  , 1.1086  , 1.3520  , 1.4699  , 1.5043  , 1.4699 ,  1.3520 ,  1.1086 ,  0.6793     ,   0],
[  0  , 1.0517  , 1.6869  , 2.0368  , 2.2035  , 2.2517  , 2.2035  , 2.0368  , 1.6869  , 1.0517   ,     0],
[        0  , 1.3601 ,  2.1574  , 2.5875  , 2.7894  , 2.8474 ,  2.7894  , 2.5875  , 2.1574  , 1.3601     ,   0],
[    0 ,  1.5309 ,  2.4144 ,  2.8851 ,  3.1043  , 3.1669  , 3.1043  , 2.8851  , 2.4144  , 1.5309     ,   0],
[   0  , 1.5170 ,  2.3891,   2.8521  , 3.0670 ,  3.1283 ,  3.0670  , 2.8521  , 2.3891 ,  1.5170    ,    0],
[      0  , 1.3128  , 2.0717  , 2.4760 ,  2.6640  , 2.7177  , 2.6640 ,  2.4760  , 2.0717  , 1.3128      ,  0],
[        0  , 0.9541 ,  1.5119  , 1.8111  , 1.9510 ,  1.9910  , 1.9510 ,  1.8111 ,  1.5119  , 0.9541     ,   0],
[        0  , 0.5008  , 0.7959  , 0.9551  , 1.0298 ,  1.0512 ,  1.0298  , 0.9551  , 0.7959  , 0.5008   ,     0],
[       0    ,    0    ,    0  ,      0    ,    0   ,     0   ,     0    ,    0    ,    0     ,   0  ,      0]]
print (u)
def f(x, y):
    return np.sin(x*np.pi)*np.sin(np.pi*y)




print u

fig = plt.figure()
ax = fig.add_subplot(111)
x, y = np.meshgrid(x, y)

m = ax.pcolormesh(x, y, u, cmap=plt.cm.viridis)
cax = plt.colorbar(m, ax=ax)
cax.set_label(r'$u(x, y)$')
plt.title("Explicit scheme, h = 0.1")
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
plt.tight_layout()
plt.show()



import numpy, sys, math

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import *




def FE(alpha,u,uPrev,N):


    for x in xrange(1,N+1):
        u[x] = alpha*uPrev[x-1] + (1.0-2*alpha)*uPrev[x] + alpha*uPrev[x+1]

def forward_euler(alpha,u,N,T):
    """
    Implements the forward Euler sheme, results saved to
    array u
    """


    for t in xrange(1,T):
        FE(alpha,u[t],u[t-1],N)

def gaussian(alpha,u,N):

    d = numpy.zeros(N) + (1+2*alpha)
    b = numpy.zeros(N-1) - alpha

    for i in xrange(1,N):

        b[i-1] /= d[i-1];
        u[i] /= d[i-1]
        d[i-1] = 1.0

        u[i+1] += u[i]*alpha
        d[i] += b[i-1]*alpha

    u[N] /= d[N-1]
    d[N-1] = 1.0


    for i in xrange(N,1,-1):
        u[i-1] -= u[i]*b[i-2]


def backward_euler(alpha,u,N,T):

    for t in xrange(1,T):
        u[t] = u[t-1].copy()
        gaussian(alpha,u[t],N)

def CN(alpha,u,N,T):

    for t in xrange(1,T):
        FE(alpha/2,u[t],u[t-1],N)
        gaussian(alpha/2,u[t],N)

def g(x):

    for i in range (1,T):
        G = (-1)**i*2/(i*math.pi)*numpy.sin(i*math.pi*x)
    return x + G


N       =   98

dt      =   0.00001

T       =   100

method  =   1



u = numpy.zeros((T,N+2),numpy.double)
(x,dx) = numpy.linspace (0,1,N+2, retstep=True)

alpha = dt/(dx**2)

#initial and boundary conditions
u[0,:] = g(x)
u[0,0] = u[0,N+1] = 0.0
t = np.zeros(T+2)
for i in xrange(1,T):
    t[i] = i*dt




if   method == 1:
    forward_euler(alpha,u,N,T)
elif method == 2:
    backward_euler(alpha,u,N,T)
elif method == 3:
    crank_nicolson(alpha,u,N,T)
else:
    print "Please select method 1,2, or 3!"
    import sys
    sys.exit(0)
#print(tlist, x)

"""
for i in xrange(1,N+1):
        if x[i] == 1-dx:
            print i
            print(t[i])
            print(u[i][i])

"""

print(u[1])
print(u[T-1])
plt.figure()
plt.title("The forward Euler method")
plt.suptitle("T=1000, dx = 1/100, dt = ")
plt.xlabel('x')
plt.ylabel('u(x,t)')
plt.plot(x, u[1],"r", label='t = 1')
plt.plot(x,u[T-1], "b", label='t = T-1')
plt.savefig('FE_10')
plt.legend(loc="upper left")
plt.show()

#fig.add_axes([0,-3.0,10000,-1.0])

"""
#print(u[1],x)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x, t = np.meshgrid(x, t)
ax.plot_surface(x, t, u)

ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$t$')
ax.set_zlabel(r'$u(x, t)$')
plt.tight_layout()
plt.show()

"""
#fig = plt.figure()
#ax = plt.axes(projection='3d')
#ax.contour3D(X, Y, Z, 50, cmap='binary')
#ax.set_xlabel('x')
#ax.set_ylabel('y')
#ax.set_zlabel('z');

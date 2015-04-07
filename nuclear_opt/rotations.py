from qutip import *
from numpy import *
from numpy.linalg import *
from scipy.optimize import *
from math import sin,cos

sx = sigmax()
sy = sigmay()
sz = sigmaz()


def rotation(phi,theta):
    return (-1j*(cos(phi)*sx + sin(phi)*sy)/2.0).expm()

def rotationx(theta):
    return rotation(0,theta)

def rotationy(theta):
    return rotation(pi,theta)

def decompose2angles(U,rotfunc):
    angels = fmin(rotfunc,ndarray([1,1,1,1,1,1]),args=U)
    return angels.tolist()

def su2X2_rot(theta11,theta12,theta13,theta21,theta22,theta23,targetU):
    theta11,theta12,theta13,theta21,theta22,theta23 = angles
    rotop1 = rotationx(theta11)*rotationy(theta12)*rotationx(theta13)
    rotop2 = rotationx(theta21)*rotationy(theta22)*rotationx(theta23)
    fid = 1-abs((tensor(rotop1,rotop2).dag()*targetU).tr())
    return fid


if __name__=='__main__':
    rotx = (-1j*sx*pi/2).expm()
    targetU = tensor(rotx,rotx)
    
    angles = decompose2angles(targetU,su2X2_rot) 
    print angles
    
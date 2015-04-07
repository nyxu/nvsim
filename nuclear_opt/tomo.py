from qcexp import *
from nv_system import *
from pulse_control import *
from timeutils import *
from nvutils import *
import scipy.io as sio
import os
from numpy import *

from scipy import *
from numpy import *
from msgutils import *
from pulse_optimization_DD import *
from matlab_opt_bridge import *
import random
sx = jmat(1,'x')
sy = jmat(1,'y')
sz = jmat(1,'z')
Id = qeye(3)
sxx = tensor(sx,sx)
syy = tensor(sy,sy)
szz = tensor(sz,sz)
#     print_matrix(szz)
sx1 = Qobj([[0,1,0],[1,0,0],[0,0,0]],dims=sx.dims,shape=sx.shape)/sqrt(2)
sy1 = Qobj([[0,-1j,0],[1j,0,0],[0,0,0]],dims=sx.dims,shape=sx.shape)/sqrt(2)
sx2 = Qobj([[0,0,0],[0,0,1],[0,1,0]],dims=sx.dims,shape=sx.shape)/sqrt(2)
sy2 = Qobj([[0,0,0],[0,0,-1j],[0,1j,0]],dims=sx.dims,shape=sx.shape)/sqrt(2)

r_pi_x1= (-1j*sx1*pi/sqrt(2)).expm()
r_pi_x11=tensor(r_pi_x1,r_pi_x1)
r_hpi_x1= (-1j*sx1*pi/sqrt(2)/2).expm()
r_hpi_x11=tensor(r_hpi_x1,r_hpi_x1)


r_pi_x2= (1j*sx2*pi/sqrt(2)).expm()
r_pi_x22=tensor(r_pi_x2,r_pi_x2)
r_hpi_x2= (1j*sx2*pi/sqrt(2)/2).expm()
r_hpi_x22=tensor(r_hpi_x2,r_hpi_x2)


def DeerOp(U1,U2,H,t):
    Uev = (-1j*H*t).expm();
    return U1*Uev*U2*Uev*U1;

def RefOp(U1,U2,U3,H,t):
    Uev = (-1j*H*t).expm();
    return U3*U1*Uev*U2*Uev*U1;

def tomo_spin1():

    Uev = (-1j*szz*pi/8).expm()
#     rou = ket2dm(tensor(basis(3,2),basis(3,2)))
    rou = ket2dm(tensor(basis(3,0),basis(3,2))+tensor(basis(3,2),basis(3,0)))
    
#     rou = r_pi_x22 * rou * r_pi_x22.dag()
#     rou = r_hpi_x11 * rou * r_hpi_x11.dag()
#     rou = r_pi_x22 * rou * r_pi_x22.dag()
    print_matrix(rou,'rou init')
    
    U_hpi=r_pi_x22 * r_hpi_x11 * r_pi_x22
    U_pi = r_pi_x22*r_pi_x11*r_pi_x22
    U_total = U_hpi*Uev*U_pi*Uev*U_hpi
    rou = U_total*rou*U_total.dag()
    print_matrix(rou,'after evo state')
    print_matrix(U_total,'Utotal')
    
    
def tomo_spin_half():
    sx = jmat(0.5,'x')
    sy = jmat(0.5,'y')
    sz = jmat(0.5,'z')
    Id = qeye(2)
    sxx = tensor(sx,sx)
    syy = tensor(sy,sy)
    szz = tensor(sz,sz)
    
    r_pi_x= (1j*sx*pi).expm()
    r_pi_xx=tensor(r_pi_x,r_pi_x)
    r_hpi_x= (1j*sx*pi/2).expm()
    r_hpi_xx=tensor(r_hpi_x,r_hpi_x)
    

    Uev = (-1j*szz*pi).expm()
    rou00 = ket2dm(tensor(basis(2,0),basis(2,1)))
    print_matrix(rou00,'initial state @ |00>')


    rou = r_hpi_xx*rou00*r_hpi_xx
    rou = Uev*rou*Uev
    rou = r_hpi_xx*rou*r_hpi_xx
    print_matrix(rou,'after evo state')

    
    
if __name__=='__main__':
#     U1 = tensor(r_hpi_x1,Id)
#     U2 = r_pi_x11
#     U3 = tensor(r_pi_x1,Id)
#     U4 = tensor(Id,r_pi_x1)
#     A=0.005
#     t = 12.5
#     H=A*szz
#     op = DeerOp(U1,U2,H,t)
#     rou = ket2dm(tensor(basis(3,0),basis(3,0)))
#     rou =op* rou *op.dag()
#     print_matrix(rou)
#     refop = RefOp(U1,U3,U4,H,t)
#     rou = ket2dm(tensor(basis(3,0),basis(3,0)))
#     rou =refop * rou *refop.dag()
#     print_matrix(rou)
    
#     tomo_spin_half()
    tomo_spin1()
    
    
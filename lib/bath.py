

from qutip import *
from pylab import *
from numpy import *
from math import *
from datetime import *
from nvpara import *
from datetime import *
import random as rand

#####################################
#get bath field in 3D for the Qutip Evolution Function
def bath_x(t,args):
    i = int(round(t / float(args['step'])))
    return args['bath'].get_field_x(i)

def bath_y(t,args):
    i = int(round(t / float(args['step'])))
    return args['bath'].get_field_y(i)

def bath_z(t,args):
    i = int(round(t / float(args['step'])))
    return args['bath'].get_field_z(i)


   

######################################################################

class Bath:
    """class to describe the bath of the spin, currently only simulating a classical random fluctuation"""    
    def __init__(self,BT1,Bb,step):
        
        self.bathT1 = float(BT1)
        self.Bb = float(Bb)
        self.step = float(step)
        self.snum = int(2.0*float(BT1)/step)  
        for kk in range(self.snum):
            self.fx.append(0.0)
            self.fy.append(0.0)
            self.fz.append(0.0)
        self.renew()
    
    fx=[]
    fy=[]
    fz=[]
    
    def renew(self):
        rand.seed(datetime.now())
        x = 0.0
        y = 0.0
        z = 0.0
        for kk in range(self.snum):
            x=x*exp(-self.step/self.bathT1)+self.Bb*sqrt((1-exp(-2*self.step/self.bathT1)))*(rand.random()-0.5)*2*sqrt(3)
            y=y*exp(-self.step/self.bathT1)+self.Bb*sqrt((1-exp(-2*self.step/self.bathT1)))*(rand.random()-0.5)*2*sqrt(3)
            z=z*exp(-self.step/self.bathT1)+self.Bb*sqrt((1-exp(-2*self.step/self.bathT1)))*(rand.random()-0.5)*2*sqrt(3)
            self.fx[kk] = x
            self.fy[kk] = y
            self.fz[kk] = z
       
    def get_field_x(self,i):
        return self.fx[i]*pi *2
    
    def get_field_y(self,i):
        return self.fy[i]*pi *2
    
    def get_field_z(self,i):
        return self.fz[i]*pi *2
    
    def gen_bath_ham_list(self):
        """generator NV hamiltonian in a bath for qutip evolution function""" 
        list=[]
        #bath for electron spin
        list.append([-Ngn*NIx+ge*eSx- Cgn*CIx,bath_x])
        list.append([-Ngn*NIy+ge*eSy- Cgn*CIy,bath_y])
        list.append([-Ngn*NIz+ge*eSz- Cgn*CIz,bath_z])
        
        #bath for N nuclear spi
        
        return list
    
    def analyze_bath(self,field):
        xf=fft.fft(field)
        plot(linspace(0,len(field),len(field)),field,'b')
        xlabel('Freq')
        ylabel('Occupation probability')
        title('Bath Spectrum')
        
        show()

#####################################
#get bath field in 3D for the Qutip Evolution Function
def gauss_x(t,args):
    i = args['sim_index']
    return args['bath'].get_field_x(i) *pi *2

def gauss_y(t,args):
    i = args['sim_index']
    return args['bath'].get_field_y(i)*pi *2

def gauss_z(t,args):
    i = args['sim_index']
    return args['bath'].get_field_z(i)*pi *2

class GaussianBath():
    def __init__(self,miu,sigma,sim_num):
        
        self.miu = float(miu)
        self.sigma = float(sigma)
        
        for kk in range(sim_num):
            self.random_seqx.append(rand.gauss(miu,sigma))
            self.random_seqy.append(rand.gauss(miu,sigma))
            self.random_seqz.append(rand.gauss(miu,sigma))
    random_seqx=[]
    random_seqy=[]
    random_seqz=[]
    def get_field_x(self,i):
        return self.random_seqx[i]
    
    def get_field_y(self,i):
        return self.random_seqy[i]
    
    def get_field_z(self,i):
        return self.random_seqz[i]
    def gen_bath_ham_list(self):
        """generator NV hamiltonian in a bath for qutip evolution function""" 
        list=[]
        #bath for electron spin
        list.append([-Ngn*NIx+ge*eSx- Cgn*CIx,gauss_x])
        list.append([-Ngn*NIy+ge*eSy- Cgn*CIy,gauss_y])
        list.append([-Ngn*NIz+ge*eSz- Cgn*CIz,gauss_z])
        
        return list    
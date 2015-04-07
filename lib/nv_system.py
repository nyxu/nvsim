from qutip import *
from spin import *
from numbers import *
from nvutils import *
from spin_system import *
from numpy import *
#for all, unit is GHz or ns or Tesla
Ngn = 0.003077#GHz/T
Cgn = 0.011#GHz/T
ge = -28.02 #GHz/T
nv_spliting = 2.87
class nv_electron_spin(electron_spin):

    def __init__(self,s=1,miu=ge):
        electron_spin.__init__(self,s,miu)
        self.D = nv_spliting
                
    def get_ham(self,Bx,By,Bz):
        return self.D*self.Szsq\
          - Bz*self.miu*self.jmatz - Bx*self.miu*self.jmatx \
          - By*self.miu*self.jmaty    

class N14_nuclear_spin(nuclear_spin):
    def __init__(self,s=1,miu=Ngn):
        nuclear_spin.__init__(self,s,miu)
        self.P = -0.00496
        
    def get_ham(self,Bx,By,Bz):
        return self.P*self.Szsq\
          -Bz*self.miu*self.jmatz -Bx*self.miu*self.jmatx\
              -By*self.miu*self.jmaty 

class C13_nuclear_spin(nuclear_spin):
    
    def __init__(self,s=0.5,miu=Cgn):
        nuclear_spin.__init__(self,s,miu)
    
    def get_ham(self,Bx,By,Bz):
        return -Bz*self.miu*self.jmatz \
                  -Bx*self.miu*self.jmatx -By*self.miu*self.jmaty


class nv_system(spin_system):
    
    def __init__(self,spins,couplings):
        spin_system.__init__(self,spins,couplings)
        self.electron_spin = -2
        self.nuclear_spins = []
        for i in range(len(self.spins)):
            if isinstance(self.spins[i],nv_electron_spin):
                self.electron_spin = i
            elif isinstance(self.spins[i],nuclear_spin):
                self.nuclear_spins.append(i)
            
    def cal_eff_g(self,j,ms,Bz):
        if j == self.electron_spin or (not j in self.nuclear_spins):
            return eye(3)
        c = self.couplings[self.electron_spin][j]
        if isinstance(c,j_coupling):
            return eye(3)
        i = self.electron_spin
        B = Bz*self.spins[i].miu
        if ms == 1 :          
            k = -1.0/(self.spins[i].D-B)
        elif ms == 0:
            k = 1.0/(self.spins[i].D+B)+1.0/(self.spins[i].D-B)
        elif ms == -1:
            k = -1.0/(self.spins[i].D+B)
            
        ratio = self.spins[i].miu*k/self.spins[j].miu

        m = zeros((3,3))
        m[0][0] = c.get_coupling(0,0)
        m[0][1] = c.get_coupling(0,1)
        m[0][2] = c.get_coupling(0,2)
        m[1][0] = c.get_coupling(1,0)
        m[1][1] = c.get_coupling(1,1)
        m[1][2] = c.get_coupling(1,2)
        return eye(3) - ratio * m

#     def cal_eff_g(self,j,ms,Bz):
#         if j == self.electron_spin or (not j in self.nuclear_spins):
#             return eye(3)
#         c = self.couplings[self.electron_spin][j]
#         if isinstance(c,j_coupling):
#             return eye(3)
#         i = self.electron_spin
#         if ms == 1:
#             ratio = -self.spins[i].miu/(self.spins[j].miu*(self.spins[i].D-self.spins[i].miu*Bz))
#         elif ms==-1:
#             ratio = -self.spins[i].miu/(self.spins[j].miu*(self.spins[i].D+self.spins[i].miu*Bz))
#         else:
#             ratio = self.spins[i].miu*(1.0/(self.spins[i].D-self.spins[i].miu*Bz)+1.0/(self.spins[i].D+self.spins[i].miu*Bz))/self.spins[j].miu
#         m = zeros((3,3))
#         m[0][0] = c.get_coupling(0,0)
#         m[0][1] = c.get_coupling(0,1)
#         m[0][2] = c.get_coupling(0,2)
#         m[1][0] = c.get_coupling(1,0)
#         m[1][1] = c.get_coupling(1,1)
#         m[1][2] = c.get_coupling(1,2)
#         return eye(3) - ratio * m

#     def cal_eff_g(self,j,ms,Bz):
#         if j == self.electron_spin or (not j in self.nuclear_spins):
#             return eye(3)
#         c = self.couplings[self.electron_spin][j]
#         if isinstance(c,j_coupling):
#             return eye(3)
#         electron = self.spins[self.electron_spin]
#         nuclear = self.spins[j]
#         axx = c.get_coupling(0,0)
#         axy = c.get_coupling(0,1)
#         axz = c.get_coupling(0,2)
#         ayx = c.get_coupling(1,0)
#         ayy = c.get_coupling(1,1)
#         ayz = c.get_coupling(1,2)
#  
#         delta_n = 1.0/(electron.D-electron.miu*Bz)
#         delta_p = 1.0/(electron.D+electron.miu*Bz)
#         delta_sum = delta_n + delta_p
#         delta_diff= delta_n - delta_p
#         if ms == 0:
# #             matrix = array([[axx*delta_sum+1j*ayx*delta_diff,axy*delta_sum+1j*ayy*delta_diff,axz*delta_sum+1j*ayz*delta_diff],\
# #                       [ayx*delta_sum-1j*axx*delta_diff,ayy*delta_sum-1j*axy*delta_diff,ayz*delta_sum-1j*axz*delta_diff],\
# #                       [                              0,                               0,                             0]])
#             matrix = array([[axx*delta_sum+1j*ayx*delta_diff,axy*delta_sum+1j*ayy*delta_diff,axz*delta_sum+1j*ayz*delta_diff],\
#                       [ayx*delta_sum-1j*axx*delta_diff,ayy*delta_sum-1j*axy*delta_diff,ayz*delta_sum-1j*axz*delta_diff],\
#                       [                              0,                               0,                             0]])
#         elif ms==1:
#             matrix = array([[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]])
#         else:
#             matrix = array([[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]])
#         return eye(3) - electron.miu * matrix / nuclear.miu


    def cal_electron_spin_projector(self,ms):
        index = get_m_index(self.spins[self.electron_spin].s,ms)
        bas = basis(self.spins[self.electron_spin].s*2+1,index)
        return self.create_matrix(self.electron_spin,Qobj(ket2dm(bas),type='oper'))
    
    def cal_nulclear_eff_control_matrix(self,Bz):
        Ieff = [0,0,0]
        P = {}
        for ms in get_ms_set(self.spins[self.electron_spin].s):
            P[ms] = self.cal_electron_spin_projector(ms)
        for j in self.nuclear_spins:
            for ms in get_ms_set(self.spins[self.electron_spin].s):
                g = self.cal_eff_g(j,ms,Bz)
                Ims = [P[ms]*self.get_jmat(j,'x')*P[ms],P[ms]*self.get_jmat(j,'y')*P[ms],P[ms]*self.get_jmat(j,'z')*P[ms]]
                for i in range(len(g)):
                    for k in range(len(g[i])):
                        Ieff[i] = Ieff[i] + g[i][k]*Ims[k]*self.spins[j].miu
        return Ieff
    
    def cal_electron_eff_control_matrix(self):    
        j = self.electron_spin
        miu = self.spins[j].miu
        return [self.get_jmat(j,'x')*miu,self.get_jmat(j,'y')*miu,self.get_jmat(j,'z')*miu]
#         print'waring:by and bz is not included'
#         return [tensor((basis(3,0)*basis(3,1).dag()+basis(3,1)*basis(3,0).dag())/sqrt(2),qeye(2))*miu,\
#                 tensor(qeye(3),qeye(2))*miu,tensor(qeye(3),qeye(2))*miu]
    
    def cal_global_control_matrix(self,eff_spins=None):
        matx=0
        maty=0
        matz=0
        if eff_spins == None:
            eff_spins = range(len(self.spins))
        for i in eff_spins:#range(len(self.spins)):
            matx = matx + self.spins[i].miu*self.get_jmat(i,'x')
            maty = maty + self.spins[i].miu*self.get_jmat(i,'y')
            matz = matz + self.spins[i].miu*self.get_jmat(i,'z')        
        return [matx,maty,matz] 
     
    def get_electron_spin_local_ham(self,B0):
        if isinstance(B0,list):
            sBx,sBy,sBz = B0
        else:
            sBx=sBy=0
            sBz=B0              
        return self.create_matrix(self.electron_spin, self.spins[self.electron_spin].get_ham(sBx,sBy,sBz))  
            
#default coupling between nuclear and NV electron spin
default_NV_C13_coupling = [[5e-3,-6.3e-3,-2.9e-3],[-6.3e-3,4.2e-3,-2.3e-3],[-2.9e-3,-2.3e-3,8.2e-3]]    
default_NV_N14_coupling =  -0.00214  



    
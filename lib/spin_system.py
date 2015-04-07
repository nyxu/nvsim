from qutip import *
from spin import *
from numbers import *
from nvutils import *


class spin_system:
    
    def __init__(self,spins,coupling_list):
        self.spins = spins
        self.couplings = []
        self.II = Qobj()
        self.jmats = {'x':[],'y':[],'z':[]}
        for spin in self.spins:
            self.jmats['x'].append(None)
            self.jmats['y'].append(None)
            self.jmats['z'].append(None)
 
            self.couplings.append([])
            for spin1 in self.spins:
                self.couplings[-1].append(none_coupling())

        for coup in coupling_list:
            [i,j,c] = coup
            if isinstance(c,list) and len(c)==3:
                if isinstance(c[0],list) and len(c[0]) == 3:
                    self.couplings[i][j]=self.couplings[j][i]=tensor_coupling(c)
                elif isinstance(c[0],Number):
                    self.couplings[i][j] = self.couplings[j][i]=vector_coupling(c)
            elif isinstance(c,Number):
                self.couplings[i][j] = self.couplings[j][i]=j_coupling(c)

              
    def create_matrix(self,i,smat):
        mats = []
        for j in range(len(self.spins)):
            if j == i:
                mats.append(smat)
            else:
                mats.append(self.spins[j].II)
        return tensor(mats)
    
    
    
    def get_jmat(self,i,ss):
        if self.jmats[ss][i]==None:          
            self.jmats[ss][i]=self.create_matrix(i, jmat(self.spins[i].s,ss))
        return self.jmats[ss][i]

    
    def get_static_ham(self,B0):
        return self.get_local_static_ham(B0) + self.get_coupled_ham()
    
    def get_local_static_ham(self,B0,spins=None):
        
            
        if isinstance(B0,list):
            sBx,sBy,sBz = B0
        else:
            sBx=sBy=0
            sBz=B0
        ham = 0  
        if spins == None:
            for i in range(len(self.spins)):
                ham = ham + self.create_matrix(i, self.spins[i].get_ham(sBx,sBy,sBz))   
        else:
            for i in spins:
                ham = ham + self.create_matrix(i, self.spins[i].get_ham(sBx,sBy,sBz))   
               
        return ham

    def get_single_coupled_ham(self,i,j,coups=['xx','yy','zz','xy','yx','xz','zx','yz','zy']): 
        if isinstance(coups,str):
            coups = [coups]
        hams = {}
        if isinstance(self.couplings[i][j],tensor_coupling):
            hams['xx'] = self.get_jmat(i, 'x')*self.get_jmat(j, 'x')*self.couplings[i][j].get_coupling(0,0)
            hams['yy'] = self.get_jmat(i, 'y')*self.get_jmat(j, 'y')*self.couplings[i][j].get_coupling(1,1)
            hams['zz'] = self.get_jmat(i, 'z')*self.get_jmat(j, 'z')*self.couplings[i][j].get_coupling(2,2)
            
            hams['xy'] = self.get_jmat(i, 'x')*self.get_jmat(j, 'y')*self.couplings[i][j].get_coupling(0,1)
            hams['yz'] = self.get_jmat(i, 'y')*self.get_jmat(j, 'z')*self.couplings[i][j].get_coupling(1,2)
            hams['xz'] = self.get_jmat(i, 'x')*self.get_jmat(j, 'z')*self.couplings[i][j].get_coupling(0,2)

            hams['yx'] = self.get_jmat(i, 'y')*self.get_jmat(j, 'x')*self.couplings[i][j].get_coupling(1,0)
            hams['zy'] = self.get_jmat(i, 'z')*self.get_jmat(j, 'y')*self.couplings[i][j].get_coupling(2,1)
            hams['zx'] = self.get_jmat(i, 'z')*self.get_jmat(j, 'x')*self.couplings[i][j].get_coupling(2,0)
        
        elif isinstance(self.couplings[i][j],vector_coupling):
            hams['xx'] = self.get_jmat(i, 'x')*self.get_jmat(j, 'x')*self.couplings[i][j].get_coupling(0)
            hams['yy'] = self.get_jmat(i, 'y')*self.get_jmat(j, 'y')*self.couplings[i][j].get_coupling(1)
            hams['zz'] = self.get_jmat(i, 'z')*self.get_jmat(j, 'z')*self.couplings[i][j].get_coupling(2)
       
        elif isinstance(self.couplings[i][j],j_coupling):
            hams['zz'] =self.get_jmat(i, 'z')*self.get_jmat(j, 'z')*self.couplings[i][j].get_coupling()
        ham = 0 
        for c in coups:
            if c in hams:
                ham = ham + hams[c]
            
        return ham
    def get_coupled_ham(self):
        ham = 0  
        for i in range(len(self.spins)):
            for j in range(i,len(self.spins)):
                ham = ham + self.get_single_coupled_ham(i,j)
        return ham
    
    """get eigenstates of nv system"""
    def get_basis(self,m):
        smats=[]
        for i in range(len(self.spins)):
            smats.append(self.spins[i].get_eigenstate(m[i]))
        return tensor(smats)
    
    """get the lindblad operator from T1 and T2 for spin i"""
    def get_lindblad_op(self,i,T1,T2):
        
        gamma1 =(1/float(T1))/2
        gamma2 = gamma1
        gamma3 = ((1.0/T2)-(gamma1 + gamma2)/2.0)/2.0
        
#         i = self.__get_index(a)
        s = self.spins[i].s
        sp = self.spins[i].jmatx + complex(0,1)*self.spins[i].jmaty
        sn = self.spins[i].jmatx - complex(0,1)*self.spins[i].jmaty
        sz = self.spins[i].jmatz
        
        L1 = self.create_matrix(i,sqrt(gamma1)*sp)
        L2 = self.create_matrix(i,sqrt(gamma2)*sn)
        L3 = self.create_matrix(i,sqrt(gamma3)*sz)
    
        return [L1,L2,L3]
    

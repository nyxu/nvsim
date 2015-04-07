from qutip import *
from nvutils import *
#for all, unit is GHz or ns

class spin:
    def __init__(self,s,miu):
        self.s = s
        self.miu = miu
        self.jmatx = jmat(s,'x')
        self.jmaty = jmat(s,'y')
        self.jmatz = jmat(s,'z')
        self.II = qeye(s*2+1)
        if s > 0.5:
            self.Szsq=self.jmatz**2
            
    def get_eigenstate(self,ms):
        return basis(int(self.s*2+1),get_m_index(self.s,ms))
    

           
class electron_spin(spin):
    def __init__(self,s,gratio):
        spin.__init__(self,s,gratio)
        
class nuclear_spin(spin):
    def __init__(self,s,gratio):
        spin.__init__(self,s,gratio)
        
class spin_coupling:
    def __init__(self,coupling):
        self.coupling=coupling

class tensor_coupling(spin_coupling):
        
    def get_coupling(self,i,j):
        return self.coupling[i][j]
    
class vector_coupling(spin_coupling):
            
    def get_coupling(self,i,j=None):
        if j == None or i == j:
            return self.coupling[i]
        else:
            return 0.0
    
class j_coupling(spin_coupling):
    
    def get_coupling(self,i=None,j=None):
        if j == None or i == None or (j==i and j==2):
            return self.coupling
        else:
            return 0.0

class none_coupling(spin_coupling):   
    
    def __init__(self):
        pass
    
    def get_coupling(self,i=None,j=None):
        return 0.0
   
from qutip import *
from nvutils import *    

class Perturbation():
    def __init__(self,ham,freq,max_field):
        self.ham = ham
        self.frequency = freq
        self.field = max_field
    
    def __decompose(self,P,omega_eval,eff_digits=6):
        Q={}
        diagonal_item = 0
        for a in range(len(P)):
            diagonal_item = diagonal_item + P[a]*self.ham*P[a]
            for b in range(a+1,len(P)):
                omega_slow = round(self.frequency + (omega_eval[a]-omega_eval[b]),eff_digits)
                omega_fast = round(self.frequency - (omega_eval[a]-omega_eval[b]),eff_digits)
                product = P[a]*self.ham*P[b]
                if not omega_slow in Q:
                    Q[omega_slow]=[]
                Q[omega_slow].append(product)                    
                if not omega_fast in Q:
                    Q[omega_fast]=[]
                Q[omega_fast].append(product.dag())
        
        return diagonal_item,Q
    
    def rotating_wave_approximate(self,P,omega_eval,s):
        diagonal_item,Q=self.__decompose(P,omega_eval)
        nondiagonal_items = []
        for omega_freq in  Q:
            prod_sum = 0
            for product in Q[omega_freq]:
                
                if max_val(abs(product.full())) > 0:
                    prod_sum = prod_sum + product
                   
            tmp = [prod_sum,omega_freq]
            if isinstance(prod_sum,Qobj):
                prod_sum = abs(prod_sum.full())
            else:
                prod_sum = abs(prod_sum)
                
            if abs(s * self.field * max_val(prod_sum)) > abs(omega_freq):
                
                nondiagonal_items.append(tmp)
            elif s == 0 and max_val(prod_sum) > 0: # means collect all items
                nondiagonal_items.append(tmp)              
        #print_rotating_frame_form(diagonal_item,nondiagonal_items) 
        return diagonal_item,nondiagonal_items

def print_matrix(mful,label=''):
    if isinstance(mful,Qobj):
        mful = mful.full()
    print label,'[\r\n'
    for i in range(len(mful)):
        text = '\t\t['
        for j in range(len(mful[0])):
            text = text + (str(mful[i][j]) if abs(mful[i][j])!=0 else '0'  ) + ',\t'
        text = text +'],'
        print text
    print ']'
    
def print_rotating_frame_form(diagonal_item,nondiagonal_items):
    diag = diagonal_item.full()
    print 'diagonal items:'
    print_matrix(diag)
    print 'nondiagonal items:'
    n = len(diag)
    text = []
    for i in range(n):
        text.append([])
        for j in range(n):
            text[-1].append([])    
    
    for n in range(len(nondiagonal_items)):
        product,omega = nondiagonal_items[n]
        pX = (product + product.dag()).full()
        pY = -1j*(product - product.dag()).full()
        print 'item ',n,':'
        print '\tomega ',omega
        print '\tmatrix pX \n',
        print_matrix(pX)
        print '\tmatrix pY \n',
        print_matrix(pY)

            
class MWFieldPerturbation(Perturbation):
    def __init__(self,mwchannel):
        Perturbation.__init__(self,mwchannel.control_ham,mwchannel.frequency,mwchannel.field)
                    
class rotating_frame:  
    def __init__(self,sysdef,H0,eff_digits=6):
        self.H0 = H0
        self.sysdef = sysdef
        evals,estates = self.H0.eigenstates()
        self.P=[]
        self.omega = []
        last_val = None
        for  i in range(len(evals)):
            eval = round(evals[i],eff_digits)
            estate = estates[i]
            if last_val != None and eval==last_val:
                self.P[-1]=self.P[-1]+estate*estate.dag()
            else:
                self.P.append(estate*estate.dag())
                self.omega.append(eval)
            last_val = eval        
            
        
if __name__=='__main__':
    from nv_system import *
    nvsys = get_nv_double_spin_system(B0=[0,0,0.5])
    from pulse_control import *
    freq = nvsys.get_transition_freq([1,0],[1,1])
    print 'mwfreq=',freq
    mwchan = MWChannel('C trans',freq,0,1e-2/ge,nvsys.gen_control_hams())
    rframe = rotation_frame(nvsys.get_local_ham(),MWChannelDefinition(mwchan))
    rframe.approximate()
    rframe.print_items()
    
    
    
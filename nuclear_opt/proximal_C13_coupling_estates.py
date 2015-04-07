from qcexp import *
from nv_system import *
from pulse_control import *
from timeutils import *
from nvutils import *
default_NV_C13_coupling = [
                           [0,0,0],
                           [0,0,0],
                           [-2.9e-3,-2.3e-3,8.2e-3]]  
default_NV_C13_coupling = [
                           [  167e-3,   0,  -90.3e-3],
                           [       0,123e-3,    0],
                           [-90.3e-3,   0,  90e-3]]    

# default_NV_C13_coupling = [[5e-3,-6.3e-3,-2.9e-3],[-6.3e-3,4.2e-3,-2.3e-3],[-2.9e-3,-2.3e-3,8.2e-3]]

if __name__=='__main__':
    nvsys = nv_system([nv_electron_spin(),C13_nuclear_spin()],[[0,1,default_NV_C13_coupling]]) 
    B0=[0,0,0.05]
    ham = nvsys.get_static_ham(B0)
    hamz=nvsys.get_local_static_ham(B0[2])
    evals,estates = ham.eigenstates()
    for i in range(len(estates)):
        eval = evals[i]
        estate = estates[i]
        sestate = cal_eigenstate(hamz,estate)
        print 'state ',i,' fidelity:',cal_state_fidelity(estate,sestate)
        
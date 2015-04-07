#
# Rabi oscillations of qubit subject to a classical driving field.
#
#from pathdef import *
from qutip import *
from pulse_control import *
from nvsys import *
from runningutils import *
from comp_basis import *

def C13_Rabi_osc_me(evo_time,Bz,pow,samplingrate,T1=0,T2=0):

    nvsys = get_double_nuclear_spin_system(B0=[0,0,Bz])
    #basis = comp_basis(nvsys,[[0,0.5],[0,-0.5],[1,0.5],[1,-0.5]])
    basis = get_double_nuclear_spin_qubit_mapping(nvsys)
    
    Ht = []
    Ht.append(nvsys.get_ham() * 2*pi)

        
    #calculte the transition frequency
    state1 = nvsys.get_sz_eigenstate([0,0.5])#basis.states[0]
    state2 = nvsys.get_sz_eigenstate([0,-0.5])#basis.states[1]
    ob_ops = [ket2dm(state1)]

    mw_freq=abs(nvsys.get_transition_freq(state1,state2))
    print 'transition frequency is ',mw_freq ,'GHz' 

    #pulse control
    pulse_control = PulseController(MWChannel(mw_freq,0,'RF'),Rabi_pulse_pattern(pow))
    Ht.extend(pulse_control.gen_control_ham_list(nvsys))

    #the initial state
    rou0 = ket2dm(state1)   
    # define the time-dependence of the hamiltonian using the list-string format
    args={'step':step,'PC':pulse_control}    

    
    # callapse operator
    c_op_list=[]
    if T1 and T2:
        c_op_list = nvsys.get_lindblad_op(1,T1,T2)
    
    #options for mesolve
    opts=Odeoptions()
    opts.nsteps=1e6
    
    # evolve and system subject to the time-dependent hamiltonian
    sampling_times = samplingrate * evo_time
    tlist = linspace(0, evo_time, sampling_times)
    output=mesolve(Ht, rou0, tlist, c_op_list, ob_ops, args,opts)
    return tlist,output.expect


if __name__=='__main__':
    pow = 2e-4/Cgn
    evo_time = 1e4
    T1=1e4
    T2=5e3
    Bz = 1 #KGauss   
    tlist,expect=C13_Rabi_osc_me(evo_time,Bz,pow,1e-2,T1,T2)
    result = data2d('Nutation of C13',tlist,expect[0],'time(ns)','occupation probility ')
    show_data_2D(result)
    show_fft_2D(result)
    
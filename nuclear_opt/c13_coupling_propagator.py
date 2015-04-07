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
                           [5e-3,-6.3e-3,-2.9e-3],
                           [-6.3e-3,4.2e-3,-2.3e-3],
                           [-2.9e-3,-2.3e-3,8.2e-3]]    
default_NV_C13_coupling = [
                           [5e-3,-6.3e-3,-2.9e-3],
                           [-6.3e-3,4.2e-3,-2.3e-3],
                           [-2.9e-3,-2.3e-3,8.2e-3]]
if __name__=='__main__':
    evo_time = 50 
    nvsys = nv_system([nv_electron_spin(),C13_nuclear_spin()],[[0,1,default_NV_C13_coupling]]) 
    
    state00 = nvsys.get_basis([0,-0.5])
    state01 = nvsys.get_basis([0,0.5])
    state10 = nvsys.get_basis([1,-0.5])
    state11 = nvsys.get_basis([1,0.5])
    B0=[0,0,0.05]
    
#     static_ham = nvsys.get_electron_spin_local_ham(B0)
    static_ham = nvsys.get_static_ham(B0)   
    nv1freq = abs(get_transition_freq(static_ham,state01,state11))
    c13freq = abs(get_transition_freq(static_ham,state10,state11))
    print 'frequency of RF is ',c13freq
    print 'frequency of MW is ',nv1freq
    #nvfreq = abs(get_transition_freq(static_ham,state10,state11))

    rf1chan = MWChannel(c13freq,0,1e-3/Cgn,nvsys.cal_nulclear_eff_control_matrix())
    mw1chan = MWChannel(nv1freq,0,1e-2/ge,nvsys.cal_electron_eff_control_matrix())
    rabipulse = single_channel_pulse('MW1')
    rframe = rotating_frame(nvsys,nvsys.get_electron_spin_local_ham(B0))
    Hleft = nvsys.get_static_ham(B0)-rframe.H0
    Utar=(-1j*evo_time*tensor(jmat(1,'z'),jmat(0.5,'z'))*pi*2*default_NV_C13_coupling[2][2]).expm()
    #Utar=(-1j*tensor(qeye(3),jmat(0.5,'x'))*pi).expm()
    print 'Utar is'
    print_matrix(Utar.full())
    mwchannels={'MW1':mw1chan}

    pcontrl = pulsecontrol(single_channel_pulse('MW1'))
    s_coup = 300
    s_mw = 30
    max_freq,args,opts,Ht = gen_rotating_frame_setting(nvsys,rframe,Hleft,mwchannels,pcontrl,s_coup,s_mw)
    default_timeline.begin()
    U=propagator(Ht, evo_time, [], args,opts)
    default_timeline.mark('using state evolution')
    evaluate_times = 1
    print 'max frequency:',max_freq
    
    U1=cal_propagator(args,Ht,[0,evo_time],max_freq,thresh_err = 0.05)
    default_timeline.mark('using exponentiation')
    print 'timecost using state evolution:',default_timeline.getduration('using state evolution')
    print 'timecost using exponentiation in '+str(evaluate_times)+' times:',default_timeline.getduration('using state evolution','using exponentiation')
    P = [nvsys.cal_electron_spin_projector(1),nvsys.cal_electron_spin_projector(0)]
    print 'fidelity between U and U1:',cal_fidelity(U,U1,P)
    print 'fidelity between U and Utar:',cal_fidelity(U,Utar,P)
    
    
from qcexp import *
from nv_system import *
from pulse_control import *
from timeutils import *
from nvutils import *
from runningutils import *
from C13_optimize import *
  
default_NV_C13_coupling = [
                           [5e-3,-6.3e-3,-2.9e-3],
                           [-6.3e-3,4.2e-3,-2.3e-3],
                           [-2.9e-3,-2.3e-3,8.2e-3]]
# default_NV_C13_coupling=8.2e-3
if __name__=='__main__':
    debug=True
     
    nvsys = nv_system([nv_electron_spin(),C13_nuclear_spin()],[[0,1,default_NV_C13_coupling]]) 
    
    state00 = nvsys.get_basis([0,-0.5])
    state01 = nvsys.get_basis([0,0.5])
    state10 = nvsys.get_basis([1,-0.5])
    state11 = nvsys.get_basis([1,0.5])
    B0=[0,0,0.05]
    
    P1=nvsys.cal_electron_spin_projector(ms=1)
    H0 = nvsys.get_local_static_ham(B0[2]) + nvsys.get_single_coupled_ham(0,1,'zz')

    c13freq = abs(get_transition_freq(H0,state10,state11))
    print 'frequency of RF is ',c13freq
    RF_pow = 2e-5
    rf1chan = MWChannel(c13freq,0,abs(RF_pow/Cgn),nvsys.cal_nulclear_eff_control_matrix())
#     rf1chan = MWChannel(c13freq,0,abs(RF_pow/Cgn),nvsys.cal_global_control_matrix())
    print_matrix(rf1chan.control_ham.full(),'rf cntl ham:')
    channels={'RF':rf1chan}
#     mw1chan = MWChannel(nv1freq,0,1e-2/ge,nvsys.cal_electron_eff_control_matrix())
    rframe = rotating_frame(nvsys,H0)
    Hleft = nvsys.get_static_ham(B0)-rframe.H0
    
    pcontrl = pulsecontrol(single_channel_pulse('RF'))
    s_coup = 10
    s_mw = 10
    evo_time = float(120)
    max_freqr,argsr,optsr,Htr = gen_rotating_frame_setting(nvsys,rframe,Hleft,channels,s_coup,s_mw)   
    argsr['PULS']=pcontrl
    argsl,optsl,Htl = gen_lab_frame_setting(nvsys.get_static_ham(B0) ,nvsys.cal_global_control_matrix(),channels)
    argsl['PULS']=pcontrl
    phaselist = [0]#linspace(0,2*pi,20)
    fidelities = []
    for phase in phaselist:
        rf1chan.phase = phase
        U_rot = cal_propagator(argsr,Htr, evo_time, int(max(10,evo_time*max_freqr*10)) ,1)
#         U_rot = propagator(Htr, evo_time,[],argsr,optsr)
        rf1chan.phase = 0
#         U_lab = propagator(Htl, evo_time,[],argsl,optsl)
        U_lab = cal_propagator(argsl,Htl, evo_time, evo_time*30 ,min(int(evo_time),6))
        print_matrix(U_rot.full(),'U_rot')
        print_matrix(U_lab.full(),'U_lab before rotating')
        U_h0 = (-1j*rframe.H0*evo_time*pi*2).expm()
        evo_states = [ket2dm(state10),ket2dm(state11)]
        expect_states = [U_h0*U_rot*ket2dm(state10)*U_rot.dag()*U_h0.dag(),U_h0*U_rot*ket2dm(state11)*U_rot.dag()*U_h0.dag()]
        fide = cal_fide_by_state(U_lab,evo_states,expect_states)
        print 'fidelity @ phase ',phase,' is ',fide
        fidelities.append(fide)
    if len(fidelities) > 1 :
        plot(phaselist,fidelities)
        show()
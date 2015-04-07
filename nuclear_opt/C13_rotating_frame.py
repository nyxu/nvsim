from qcexp import *
from nv_system import *
from pulse_control import *
from timeutils import *
from nvutils import *
from runningutils import *
  
default_NV_C13_coupling = [
                           [5e-3,-6.3e-3,-2.9e-3],
                           [-6.3e-3,4.2e-3,-2.3e-3],
                           [-2.9e-3,-2.3e-3,8.2e-3]]
global debug
default_NV_C13_coupling=8.2e-3
if __name__=='__main__':
    debug=True
     
    nvsys = nv_system([nv_electron_spin(),C13_nuclear_spin()],[[0,1,default_NV_C13_coupling]]) 
    
    state00 = nvsys.get_basis([0,-0.5])
    state01 = nvsys.get_basis([0,0.5])
    state10 = nvsys.get_basis([1,-0.5])
    state11 = nvsys.get_basis([1,0.5])
    B0=[0,0,0.05]
    
    P1=nvsys.cal_electron_spin_projector(ms=1)
#     H0 = nvsys.get_electron_spin_local_ham(B0[2]) #+ nvsys.get_single_coupled_ham(0,1,'zz')
    H0 = nvsys.get_local_static_ham(B0[2]) + nvsys.get_single_coupled_ham(0,1,'zz')

    c13freq = abs(get_transition_freq(H0,state10,state11))
    print 'frequency of RF is ',c13freq
    RF_pow = 1e-4
    rf1chan = MWChannel(c13freq,0,abs(RF_pow/Cgn),nvsys.cal_nulclear_eff_control_matrix(B0[2]))
    print_matrix(rf1chan.control_ham.full(),'rf cntl ham:')
    channels={'RF':rf1chan}
#     mw1chan = MWChannel(nv1freq,0,1e-2/ge,nvsys.cal_electron_eff_control_matrix())
    rframe = rotating_frame(nvsys,H0)
    Hleft = nvsys.get_static_ham(B0)-rframe.H0
    
    pcontrl = pulsecontrol(single_channel_pulse('RF'))
    s_coup = 100
    s_mw = 1
#     max_freq,args,opts,Ht = gen_rotating_frame_setting(nvsys,rframe,Hleft,channels,s_coup,s_mw)
    args,opts,Ht = gen_lab_frame_setting(nvsys.get_static_ham(B0) ,nvsys.cal_global_control_matrix(),channels)
    print 'Ht[0] is '
    for i in range(1,len(Ht)):
        print_matrix(Ht[i][0].full(),str(i-1)+'th:')    
    evo_time = 500
    
    sampling_times = 1 * evo_time
    
    tlist = linspace(0, evo_time, sampling_times)
    rou0 = ket2dm(state10)
    args['PULS']=pcontrl
    output=mesolve(Ht, rou0, tlist, [], rou0, args,opts)

    result = data2d('C13 Rabi OSC @ Pow ' + str(RF_pow*1e6) + 'KHz,s=('+str(s_coup)+','+str(s_mw)+')',\
                    tlist,output.expect[0],'time(ns)','overlap @ |1,1/2> ')
    show_data_2D(result)
    show_fft_2D(result)

    
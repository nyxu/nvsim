from qcexp import *
from nv_system import *
from pulse_control import *
from timeutils import *
from nvutils import *
from runningutils import *
  

if __name__=='__main__':

    B0=[0,0,0.5]
    coupling = [2.2e-3,2.2e-3,2.2e-3]
    print coupling
    nvsys = nv_system([nv_electron_spin(),N14_nuclear_spin()],[[0,1,coupling]]) 
    ms = 0
    state0 = nvsys.get_basis([ms,0])
    state1 = nvsys.get_basis([ms,1])
    print nvsys.cal_eff_g(1,0,B0[2])
    print nvsys.cal_eff_g(1,1,B0[2])
    print nvsys.cal_eff_g(1,-1,B0[2])

    
    H0 = nvsys.get_local_static_ham(B0[2],[0]) #+ nvsys.get_single_coupled_ham(0,1,'zz')
    
    N14freq = abs(get_transition_freq(nvsys.get_local_static_ham(B0[2])+ nvsys.get_single_coupled_ham(0,1,'zz'),state0,state1))
    print 'frequency of RF is (MHz)',N14freq*1e3#,N14freq1*1e3
    RF_pow = 2e-5
    RF0 = 'RF0'
    rf0chan = MWChannel(N14freq,0,abs(RF_pow/Ngn),nvsys.cal_nulclear_eff_control_matrix(B0[2]))
    channels={RF0:rf0chan}  
    
    rframe = rotating_frame(nvsys,H0)
    Hleft = nvsys.get_static_ham(B0)-rframe.H0
    s_coup,s_mw = 300,0
    max_freq,args,opts,Ht = gen_rotating_frame_setting(nvsys,rframe,Hleft,channels,s_coup,s_mw,debug=True)

    evo_time = 10000.0
    sampling_times = 10 * evo_time
    tlist = linspace(0, evo_time, sampling_times)

    rabi_pulse = single_channel_pulse(RF0,evo_time)
    args['PULS']=rabi_pulse
    
    rou0 = ket2dm(state0)
    output1=mesolve(Ht, rou0, tlist, [], rou0, args,opts)
    
        
  
    rf0chan = MWChannel(N14freq,0,abs(RF_pow/Ngn),nvsys.cal_nulclear_eff_control_matrix(B0[2]))
    channels={RF0:rf0chan}  
   
    args,opts,Ht = gen_lab_frame_setting(nvsys.get_static_ham(B0) ,nvsys.cal_global_control_matrix(),channels) 
    args['PULS']=rabi_pulse
    output2=mesolve(Ht, rou0, tlist, [], rou0, args,opts)
      
    
    plot(tlist,output1.expect[0],label='rotation frame')
    plot(tlist,output2.expect[0],label='lab frame')
    title('N14 enhancement effect in ms = '+str(ms)+ ' @ Bz='+str(B0[2]*1000)+'Gauss with RF Pow ' + str(RF_pow*1e6)+ 'KHz')
    legend()
    show()

    
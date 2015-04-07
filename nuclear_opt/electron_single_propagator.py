from qcexp import *
from nv_system import *
from pulse_control import *
from timeutils import *
from nvutils import *
default_NV_C13_coupling = [
                           [0,0,0],
                           [0,0,0],
                           [-2.9e-3,-2.3e-3,8.2e-3]]  
none_NV_C13_coupling = 0   
default_NV_C13_coupling = [
                           [5e-3,-6.3e-3,-2.9e-3],
                           [-6.3e-3,4.2e-3,-2.3e-3],
                           [-2.9e-3,-2.3e-3,8.2e-3]]
default_NV_C13_coupling=10e-3

if __name__=='__main__':
    nvsys1 = nv_system([nv_electron_spin(),C13_nuclear_spin()],[[0,1,default_NV_C13_coupling]]) 
    nvsys2 = nv_system([nv_electron_spin(),C13_nuclear_spin()],[[0,1,none_NV_C13_coupling]]) 
    
    state00 = nvsys1.get_basis([0,-0.5])
    state01 = nvsys1.get_basis([0,0.5])
    state10 = nvsys1.get_basis([1,-0.5])
    state11 = nvsys1.get_basis([1,0.5])
    B0=[0,0,0.05]
    
    static_ham1 = nvsys1.get_electron_spin_local_ham(B0)
    static_ham2 = nvsys2.get_electron_spin_local_ham(B0)

    nv1freq = abs(get_transition_freq(static_ham1,state01,state11))

    print 'frequency of MW is ',nv1freq

    mw1chan = MWChannel(nv1freq,0,40e-3/sqrt(2)/ge,nvsys1.cal_electron_eff_control_matrix())
    pcontrol = pulsecontrol(single_channel_pulse('MW1'))
    rframe1 = rotating_frame(nvsys1,nvsys1.get_electron_spin_local_ham(B0))
    rframe2 = rotating_frame(nvsys2,nvsys2.get_electron_spin_local_ham(B0))
    Hleft1 = nvsys1.get_static_ham(B0) - rframe1.H0
    Hleft2 = 0
    channels = {'MW1':mw1chan}

    s_coup = 300
    s_mw = 1
    angles = []
    fidelities = []

    max_freq1,args1,opts,Ht1 = gen_rotating_frame_setting(nvsys1,rframe1,Hleft1,channels,pcontrol,s_coup,s_mw)
    max_freq2,args2,opts,Ht2 = gen_rotating_frame_setting(nvsys2,rframe2,Hleft2,channels,pcontrol,s_coup,s_mw)
    print Ht1[0]/pi/2
    print Ht2[0]/pi/2
    print Ht1[1][0]/pi/2
    print Ht2[1][0]/pi/2  
    print Ht1[2][0]/pi/2
    print Ht2[2][0]/pi/2    
    for evo_time in linspace(0.1,500,500):
          
#         U1=propagator(Ht1, evo_time, [], args1,opts)
#           
#         U2=propagator(Ht2, evo_time, [], args2,opts)
        U1=cal_propagator(args1,Ht1, [0,evo_time],max_freq1 )
            
        U2=cal_propagator(args2,Ht2, [0,evo_time],max_freq2)
      
        P = [nvsys1.cal_electron_spin_projector(1),nvsys1.cal_electron_spin_projector(0)]
  
        fidelity = cal_fidelity(U1,U2)
        print 'fidelity between U and Utar:',fidelity
        angles.append(evo_time )
        fidelities.append(fidelity)
  
    print_matrix(U1.full())
    print_matrix(U2.full())
    plot(angles,fidelities)
    ylabel('fidelity')
    #xlabel('rotation angle (2pi)')
    xlabel('time')
    show()
    
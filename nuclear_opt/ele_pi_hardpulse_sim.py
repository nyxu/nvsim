from qcexp import *
from nv_system import *
from pulse_control import *
from timeutils import *
from nvutils import *


default_NV_C13_coupling = [
                           [5e-3,-6.3e-3,-2.9e-3],
                           [-6.3e-3,4.2e-3,-2.3e-3],
                           [-2.9e-3,-2.3e-3,8.2e-3]]
sx01 = basis(3,0)*basis(3,1).dag()+basis(3,1)*basis(3,0).dag()
Upi=tensor((-1j*sx01*pi/2).expm(),qeye(2))
B0=[0,0,0.05]

def diff_propagator(nvsys,rframe,Hleft,channels,pcontrol,pulse_duration,s_coup=300,s_mw=1):
    #ratating frame approximation
    max_freqr,argsr,optsr,Htr = gen_rotating_frame_setting(nvsys,rframe,Hleft,channels,s_coup,s_mw)    
    argsr['PULS']= pcontrol 
    Ur=propagator(Htr, pulse_duration, [], argsr,optsr)
    fidelity_r = cal_fidelity(Upi,Ur)

    #lab frame settings
    argsl,optsl,Htl = gen_lab_frame_setting(nvsys.get_static_ham(B0) ,nvsys.cal_electron_eff_control_matrix(),channels)
    argsl['PULS']= pcontrol 
    Ulpi=propagator(Htl,pulse_duration,[],argsl,optsl)
    Ulpir = (1j*rframe.H0*pulse_duration*pi*2).expm()*Ulpi
    fidelity_l = cal_fidelity(Upi,Ulpir)
#     fidelity_l = 1
    return fidelity_r,fidelity_l

if __name__=='__main__':
    
    nvsys = nv_system([nv_electron_spin(),C13_nuclear_spin()],[[0,1,default_NV_C13_coupling]]) 
    nv1freq = abs(get_transition_freq(nvsys.get_electron_spin_local_ham(B0),nvsys.get_basis([0,0.5]),nvsys.get_basis([1,0.5])))
    print 'frequency of MW is ',nv1freq
    


    #rotating frame    
    rframe = rotating_frame(nvsys,nvsys.get_electron_spin_local_ham(B0))
    Hleft = nvsys.get_static_ham(B0) - rframe.H0
    s_coup = 300
    s_mw = 1
    evolist = linspace(10e-3,80e-3,50)
    fidelities_r = []
    fidelities_l = []
    for power in evolist:
        #MW settings
        MW_Power = power #GHz       
        mw1chan = MWChannel(nv1freq,0,abs(MW_Power/sqrt(2)/ge),nvsys.cal_electron_eff_control_matrix())
        channels = {'MW1':mw1chan}
        
        #ordinary pulse
        pi_len = 0.5/MW_Power 
        pcontrol = pulsecontrol(single_channel_pulse('MW1'))
    
        fidelity_r,fidelity_l=diff_propagator(nvsys,rframe,Hleft,channels,pcontrol,pi_len)
        print 'fidelity of pi by ordinary pulse in Rotating:',fidelity_r,'; versus Labframe: ',fidelity_l,' @ Pow '+str(MW_Power*1e3)+"MHz"
        fidelities_r.append(fidelity_r)
        fidelities_l.append(fidelity_l)
        
    plot(evolist,fidelities_r,label='rotating frame')
    plot(evolist,fidelities_l,label='lab frame')
    xlabel('MW Power (GHz)')
    ylabel('Fidelity')
    title('Hard Pulse fidelity for pi Operation with Coupling '+str(default_NV_C13_coupling[2][2]*1e3)+"MHz")
    legend()
    show()
    #optimized shaped pulse
#     slotwidth=4.000000e-01
#     x=[-1.595632e+00,-8.383072e-01,9.201577e-01,9.571027e-01,9.789834e-01,9.903642e-01,9.999999e-01,9.999998e-01,9.999997e-01,9.999996e-01,9.999985e-01,9.999988e-01,9.999992e-01,9.999992e-01,9.999992e-01,9.999967e-01,9.999995e-01,9.999981e-01,9.999952e-01,9.999991e-01,9.999988e-01,9.999999e-01,9.999978e-01,9.999964e-01,1.000000e+00,9.999989e-01,9.999969e-01,9.999971e-01,9.999974e-01,9.999966e-01,9.999991e-01,9.999983e-01,9.999975e-01,9.999974e-01,1.000000e+00,9.999987e-01,9.999971e-01,9.999975e-01,9.999993e-01,9.999999e-01,9.838592e-01,9.468232e-01,9.123810e-01,8.803658e-01,8.845804e-01,8.624973e-01,7.984867e-01,7.865275e-01,6.982913e-01,6.405604e-01]
#     y=[-3.778873e-14,-2.008214e-14,5.310598e-15,2.497353e-15,3.754872e-15,4.356397e-15,1.988300e-15,1.388338e-15,1.169134e-15,8.038737e-16,6.404175e-16,6.783290e-16,1.529723e-16,-5.822182e-17,-1.972737e-16,-3.236137e-16,-1.694290e-16,-5.745780e-16,-1.448177e-15,-1.838654e-15,-1.613799e-15,-1.903146e-15,-1.697935e-15,-1.596947e-15,-1.349174e-15,-1.748114e-15,-1.075400e-15,-7.819280e-17,-4.757497e-16,-4.889154e-16,-1.874490e-15,-1.715770e-15,-1.621171e-15,-1.084460e-15,4.696479e-16,9.676666e-16,-1.564331e-15,2.561327e-15,9.208595e-16,2.125328e-15,7.302121e-16,-3.676287e-16,-3.526079e-15,-3.766677e-15,-2.822590e-15,-5.106410e-15,-5.924277e-15,-3.753896e-15,-9.011255e-15,-1.064802e-14]
#     grape_pulse = shaped_pulse(['MW1',x,y],slotwidth)
#     fidelity_r,fidelity_l=diff_propagator(nvsys,rframe,Hleft,channels,grape_pulse,slotwidth*len(x))
#     print 'fidelity of pi by GRAPE pulse in Rotating:',fidelity_r,'; versus Labframe: ',fidelity_l
#     
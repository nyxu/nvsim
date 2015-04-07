from qcexp import *
from nv_system import *
from pulse_control import *
from timeutils import *
from nvutils import *
from ele_pi_hardpulse_sim import *



if __name__=='__main__':
    sx01 = basis(3,0)*basis(3,1).dag()+basis(3,1)*basis(3,0).dag()
    Upi=tensor((-1j*sx01*pi/2).expm(),qeye(2))
    fidelities_r = []
    fidelities_l = []
    evolist=linspace(0,2e-2)
    for coupling in evolist:
        nvsys = nv_system([nv_electron_spin(),C13_nuclear_spin()],[[0,1,coupling]]) 
        nv1freq = abs(get_transition_freq(nvsys.get_electron_spin_local_ham(B0),nvsys.get_basis([0,0.5]),nvsys.get_basis([1,0.5])))
#         print 'frequency of MW is ',nv1freq
        #rotating frame    
        rframe = rotating_frame(nvsys,nvsys.get_electron_spin_local_ham(B0))
        Hleft = nvsys.get_static_ham(B0) - rframe.H0
        s_coup = 300
        s_mw = 1
        #MW settings
        MW_Power = 40e-3 #GHz       
        mw1chan = MWChannel(nv1freq,0,abs(MW_Power/sqrt(2)/ge),nvsys.cal_electron_eff_control_matrix())
        channels = {'MW1':mw1chan}
        
        #ordinary pulse
        pi_len = 0.5/MW_Power 
        pcontrol = pulsecontrol(single_channel_pulse('MW1'))
    
        fidelity_r,fidelity_l=diff_propagator(nvsys,rframe,Hleft,channels,pcontrol,pi_len)
        print 'fidelity of pi by ordinary pulse in Rotating:',fidelity_r,'; versus Labframe: ',fidelity_l,' @ Coupling '+str(coupling*1e3)+"MHz"
        fidelities_r.append(fidelity_r)
        fidelities_l.append(fidelity_l)
        
    plot(evolist,fidelities_r,label='rotating frame')
    plot(evolist,fidelities_l,label='lab frame')
    xlabel('Coupling (GHz)')
    ylabel('Fidelity')
    title('Hard Pulse fidelity for pi Operation @ Power '+str(MW_Power*1e3)+"MHz")
    legend()
    show()

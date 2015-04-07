from qcexp import *
from nv_system import *
from pulse_control import *
from timeutils import *
from nvutils import *
if __name__=='__main__':
    evo_time = 10000
    sampling_rate = 1
    
    nvsys = nv_system([nv_electron_spin(),C13_nuclear_spin()],[[0,1,default_NV_C13_coupling]]) 
    
    state00 = nvsys.get_basis([0,-0.5])
    state01 = nvsys.get_basis([0,0.5])
    state10 = nvsys.get_basis([1,-0.5])
    state11 = nvsys.get_basis([1,0.5])
    B0=[0,0,0.05]
    
    static_ham = nvsys.get_static_ham(B0)   
    c13freq = abs(get_transition_freq(static_ham,state00,state01))
    nvfreq = abs(get_transition_freq(static_ham,state10,state11))

    mw1chan = MWChannel(c13freq,0,1e-4/Cgn)
    pulse = pulse_sequence(pulse_shape())
    rframe = rotating_frame(nvsys,nvsys.get_local_static_ham(B0))
    
    exp = qc_experiment(nvsys,{'MW1':mw1chan},{'MW1':pulse},rframe,B0)
    init_state = ket2dm(state00)
    ob_op = ket2dm(state00)
    s_coup = 300
    s_mw = 300
#     opts,Ht = gen_rotating_frame_setting(exp,s_coup,s_mw)
    opts,Ht = gen_lab_frame_setting(exp)
    
    default_timeline.begin()
    default_timeline.mark('rotation frame ended')
    
    print 'rotation frame simulation cost:',default_timeline.getduration('rotation frame ended')
    title('lab frame Rabi Osc')
    tlist = linspace(0,evo_time,evo_time* sampling_rate)
    output=mesolve(Ht, init_state, tlist, [], [ob_op], exp,opts)
    default_timeline.mark('lab frame ended')
    print 'lab frame simulation cost:',default_timeline.getduration('rotation frame ended','lab frame ended')
    plot(tlist,real(output.expect[0]))
    show()
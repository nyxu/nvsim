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
default_NV_C13_coupling = [[5e-3,-6.3e-3,-2.9e-3],[-6.3e-3,4.2e-3,-2.3e-3],[-2.9e-3,-2.3e-3,8.2e-3]]

if __name__=='__main__':
    evo_time = 40 
    nvsys = nv_system([nv_electron_spin(),C13_nuclear_spin()],[[0,1,default_NV_C13_coupling]]) 
    
    B0=[0,0,0.05]
    rframe = rotating_frame(nvsys,nvsys.get_electron_spin_local_ham(B0))

    exp = exp_setup(nvsys,{},rframe,B0)
    #exp = experiment_settings(nvsys,{},{},rframe,B0)
    evo_range = linspace(0.5,40,40)
    fidelities = []
    for evo_time in evo_range:
        s_coup = 300
        max_freq,args,opts,Ht = gen_rotating_frame_setting(exp,{},s_coup,0)
        U1=cal_propagator(args,Ht,[0,evo_time],max_freq,thresh_err = 0.05)
        s_coup = 0
        max_freq,args,opts,Ht = gen_rotating_frame_setting(exp,{},s_coup,0)
        U2=cal_propagator(args,Ht,[0,evo_time],max_freq,thresh_err = 0.05)
        fid = cal_fidelity(U1,U2)
        print 'fidelity between U1 and U2 @ time ',str(evo_time),':', fid
        fidelities.append(fid)   
    title('evolution in different rotation wave appr')
    xlabel('evo_time')
    ylabel('fidelity')
    plot(evo_range,fidelities)
    
    show()
    print 'fidelity norm is :',cal_fidelity(U1,U1)
    
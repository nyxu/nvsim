from qcexp import *
from nv_system import *
from pulse_control import *
from timeutils import *
from nvutils import *
import scipy.io as sio
import os
from numpy import *
from scipy import fftpack
from numpy import *
from msgutils import *
from C13_optimize import *

coupling = [[5e-3,-6.3e-3,-2.9e-3],
           [-6.3e-3,4.2e-3,-2.3e-3],
           [-2.9e-3,-2.3e-3,8.2e-3]]


    
if __name__=='__main__':
    
    nvsys = nv_system([nv_electron_spin(),C13_nuclear_spin()],[[0,1,coupling]]) 
    state00 = nvsys.get_basis([0,-0.5])
    state01 = nvsys.get_basis([0,0.5])
    state10 = nvsys.get_basis([1,-0.5])
    state11 = nvsys.get_basis([1,0.5])
    B0=[0,0,0.05]
    H0 = nvsys.get_local_static_ham(B0[2]) + nvsys.get_single_coupled_ham(0,1,'zz')
    
    c13freq = abs(get_transition_freq(H0,state10,state11))
    print 'frequency of RF is ',c13freq
    RF_pow = 6e-5
    rf1chan = MWChannel(c13freq,0,abs(RF_pow/Cgn),nvsys.cal_nulclear_eff_control_matrix())
    channels={'RF':rf1chan} 
    
    P1 = nvsys.cal_electron_spin_projector(1)
    Ix = tensor(basis(3,0)*basis(3,0).dag(),sigmax())
    targetU = (-1j*Ix*pi/2).expm()
    
    rframe = rotating_frame(nvsys,H0)
    Hleft = nvsys.get_static_ham(B0)-rframe.H0
    s_coup,s_mw = 30,0.1
    max_freq,args,opts,Ht = gen_rotating_frame_setting(nvsys,rframe,Hleft,channels,s_coup,s_mw)
#     print nvsys.cal_eff_g(1,1)
    print_matrix(Ht[0].full(),'static ham')
    for i in range(1,len(Ht)):
        print_matrix(Ht[i][0].full(),str(i-1)+'th:')
    period = 1.0/coupling[2][2]
    oper_time = int(2000.0/period)*period
    nSlots = 500.0
    x,y=optimize_C13_operator(targetU,Ht,oper_time,nSlots,args,P1)

    print 'in rotating frame'
    grape_pulse = shaped_pulse(['RF',x,y],oper_time/nSlots)
    args['PULS']=grape_pulse
    U_rot = cal_propagator(args,Ht,oper_time,nSlots,1)
    print_matrix(U_rot.full(),'U=')
    fider = cal_fide(U_rot,targetU,P1)

    print 'in lab frame'    
    args,opts,Ht = gen_lab_frame_setting(nvsys.get_static_ham(B0) ,nvsys.cal_global_control_matrix(),channels) 
    args['PULS']=grape_pulse
    U_lab = cal_propagator(args,Ht,oper_time,oper_time*30)
    U_h0 = (-1j*rframe.H0*oper_time*pi*2).expm()
    evo_states = [ket2dm(state10),ket2dm(state11)]
    expect_states = [U_h0*targetU*ket2dm(state10)*targetU.dag()*U_h0.dag(),\
                     U_h0*targetU*ket2dm(state11)*targetU.dag()*U_h0.dag()]
    print_matrix(U_lab.full(),'U_lab=')
    fidel = cal_fide_by_state(U_lab,evo_states,expect_states)
    qsave(U_lab,'U_c13_lab')
    msg = 'fidelity @ power '+str(RF_pow) + ' (r/l):' +str(fider) +','+str(fidel)
    print msg
    send_result(['nyxu@ustc.edu.cn'],msg,'nyxu@mail.ustc.edu.cn', 'trapdoor') 

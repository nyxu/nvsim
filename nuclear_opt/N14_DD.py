from qcexp import *
from nv_system import *
from pulse_control import *
from timeutils import *
from nvutils import *
import scipy.io as sio
import os
from numpy import *

from scipy import *
from numpy import *
from msgutils import *
from pulse_optimization_DD import *
from matlab_opt_bridge import *
import random


def opt_pulse_N14(args):
    targetU,target_fide,time_ratio=args
    coupling = [-2.2e-3,-2.2e-3,-2.2e-3]
    B0=[0,0,0.5]
    RF_pow = 1.5e-6
    MW_pow = 1e-3
    electron_pi_len = 0.5/(MW_pow * sqrt(2))   
    
    total_time = abs(time_ratio*1.0/coupling[2])
    nSlots = 400.0
    
    
    nvsys = nv_system([nv_electron_spin(),N14_nuclear_spin()],[[0,1,coupling]]) 
    e_state0 = 0
    e_state1 = -1
    state00 = nvsys.get_basis([e_state0,0])
    state01 = nvsys.get_basis([e_state0,1])
    state10 = nvsys.get_basis([e_state1,0])
    state11 = nvsys.get_basis([e_state1,1])
    
    P = nvsys.cal_electron_spin_projector(e_state0)+nvsys.cal_electron_spin_projector(e_state1)

    H0 = nvsys.get_local_static_ham(B0[2])# + nvsys.get_single_coupled_ham(0,1,'zz')
#     H0 = nvsys.get_electron_spin_local_ham(B0[2]) 

    N14freq0 = abs(get_transition_freq(H0 + nvsys.get_single_coupled_ham(0,1,'zz'),state00,state01))
    N14freq1 = abs(get_transition_freq(H0 + nvsys.get_single_coupled_ham(0,1,'zz'),state10,state11))
    electron_ns_freq = abs(get_transition_freq(nvsys.get_local_static_ham(B0[2],spins=[0]),state00,state10))
    print 'frequency of RF is (MHz)',N14freq0*1e3,N14freq1*1e3
    print 'frequency of MW is (MHz)',electron_ns_freq*1e3

    rf0chan = MWChannel(N14freq0,0,abs(RF_pow/Ngn),nvsys.cal_nulclear_eff_control_matrix(B0[2]))
    rf1chan = MWChannel(N14freq1,0,abs(RF_pow/Ngn),nvsys.cal_nulclear_eff_control_matrix(B0[2]))
    mwchan = MWChannel(electron_ns_freq,0,abs(MW_pow/ge),nvsys.cal_electron_eff_control_matrix())
    
    optimization_channels={RF0:rf0chan,RF1:rf1chan} 
    sim_channels={RF0:rf0chan,RF1:rf1chan,MW:mwchan} 

 
    Sx01 = MW_pow*pi*2*\
        tensor(basis(3,1)*basis(3,2).dag()+basis(3,2)*basis(3,1).dag(),basis(3,1)*basis(3,1).dag()+basis(3,0)*basis(3,0).dag())\
        /sqrt(2)
    Sy01 = MW_pow*pi*2j*\
        tensor(basis(3,1)*basis(3,2).dag()-basis(3,2)*basis(3,1).dag(),basis(3,1)*basis(3,1).dag()+basis(3,0)*basis(3,0).dag())\
        /sqrt(2)
    print_matrix(Sy01,'Sy01=')
 
       
    rframe = rotating_frame(nvsys,H0)
    Hleft = nvsys.get_static_ham(B0)-rframe.H0
    s_coup,s_mw = 300,300

    tau,tlist,tslices,init_controls,control_mask,electron_pulse,H_eslices=cal_xy4_settings(total_time,electron_pi_len,nSlots,optimization_channels,Sx01,Sy01)
     
    max_freq_opt,args_opt,opts_opt,Ht_opt = gen_rotating_frame_setting(nvsys,rframe,Hleft,optimization_channels,s_coup,s_mw,debug=True)
       
    H_static, ham_static_nondiagnal, chamsx, chamsy  = gen_ham_series_for_matlab(nvsys,rframe,Hleft,optimization_channels,s_coup,s_mw,tlist,P)
 
    Ht_opt[0] = P*Ht_opt[0]*P
    for i in range(1,len(Ht_opt)):
        Ht_opt[i][0] = P*Ht_opt[i][0]*P
    
     
       
    x,y,dict=optimize_operator_DD(targetU,target_fide,H_static, ham_static_nondiagnal,H_eslices, chamsx, chamsy,P,init_controls,control_mask,tslices)
    grape_pulse = sliced_pulse([[RF0,x[RF0],y[RF0]],[RF1,x[RF1],y[RF1]]],tslices,tlist)
    total_pulse = multichannel_pulse([grape_pulse,electron_pulse])
    print 'in rotating frame'
    s_coup,s_mw = 300,300
    max_freq_sim,args_sim,opts_sim,Ht_sim = gen_rotating_frame_setting(nvsys,rframe,Hleft,sim_channels,s_coup,s_mw,debug=True)
    args_sim['PULS']=total_pulse      
            
    U_rot = cal_propagator(args_sim,Ht_sim,tlist,1)
    print_matrix(U_rot.full(),'U=')
    fider = real((U_rot.dag()*targetU).tr()/9.0)
    print 'fidelity by tr is ',fider

    return fider  
       
if __name__=='__main__':
    thelist = [pi/2,pi/3]
    timelist = [40.0,40]
    Ix = tensor(basis(3,1)*basis(3,1).dag()+basis(3,2)*basis(3,2).dag(),basis(3,0)*basis(3,1).dag()+basis(3,1)*basis(3,0).dag())
    arglist=[]
    for i in range(len(thelist)):
        targetU = (-1j*Ix*thelist[i]).expm()
        arglist.append([targetU,0.994,timelist[i]])
    fider=opt_pulse_N14(arglist[0])
#     reslist=parfor(opt_pulse_N14,arglist)
#     result = ''
#     for i in range(len(thelist)):
#         fider,fide_diff = reslist[i]
#         result = result+str(thelist[i])+':'+str(fider)+','+str(fide_diff)
# 	send_result(['nanyang.xu@gmail.com'],result,'nyxu@mail.ustc.edu.cn', 'trapdoor')
#
#     send_result(['nanyang.xu@gmail.com'],result,'nyxu@mail.ustc.edu.cn', 'trapdoor')
    send_result(['nanyang.xu@gmail.com'],str(fider),'nyxu@mail.ustc.edu.cn', 'trapdoor')
    
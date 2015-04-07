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
from pulse_optimization import *
from matlab_opt_bridge import *
import random

RF0 = 'RF0'
RF1 = 'RF1'
MW  = 'MW'
import StringIO          
def optimize_operator_DD(targetU,target_fide,H_static,H_nondiagnal,H_electron,ctlx,ctly,P,init_controls,control_mask,time_slices,max_power=1.0):
    
    matlabcmd = '/Applications/MATLAB_R2012a.app/bin/matlab -nodesktop -nosplash -nojvm -wait -r pulse_optimization_DD  '
    matdir = '/Users/nyxu/Dropbox/workspace/MATLAB/custemed_DYNAMO_v_1_2/'
    opfile = matdir + 'pulse_optimization.mat'
    pulsefile =  matdir + 'optimized_pulse.mat'
    
    
    Hnd = []
    for H in H_nondiagnal:
        Hnd.append(H.trans().full() if isinstance(H,Qobj) else zeros(targetU.shape))
    Hnd = transpose(array(Hnd))
    
    Hfree = []
    for H in H_electron:
        Hfree.append(H.trans().full() if isinstance(H,Qobj) else zeros(targetU.shape))
    Hfree = transpose(array(Hfree))
    
    ham_ctl = []
    for i in ctlx:
        hx = []
        hy = []
        for j in range(len(ctlx[i])):
            hx.append(ctlx[i][j].trans().full())
            hy.append(ctly[i][j].trans().full())
        ham_ctl.append(hx)
        ham_ctl.append(hy)
    ham_ctl_array=transpose(array(ham_ctl))
    dhams = {}
        
    dhams['ham_nondiagnal'] = Hnd
    dhams['ham_free'] = Hfree
    dhams['initial_controls']=init_controls
    dhams['ext_control_mask']=control_mask
    dhams['ndims']=targetU.shape[0]
    dhams['ham_control'] = ham_ctl_array
    dhams['Utar']=targetU.full()
    dhams['target_fide']=target_fide
    dhams['tslices'] = time_slices
    dhams['H_static']=H_static.full()
    dhams['total_time']=sum(time_slices)
    dhams['nSlots']=len(time_slices)
    dhams['max_power'] = max_power
    sio.savemat(opfile,dhams)
    ret = os.system(matlabcmd)
    x = {}
    y = {}
    if ret == 0 and os.path.isfile(pulsefile):
        dict = sio.loadmat(pulsefile)
#         print dict.keys()
        controls = dict['controls']  
        k = 0      
        for name in ctlx:
#             name=channel_ctl_hams[k][0]
            x[name] = []
            y[name] = []
            for i in range(len(controls)):
                x[name].append(controls[i][k*2])
                y[name].append(controls[i][k*2+1])
            k = k + 1
#     for name in ctlx:
#         clear_mfile(matdir,name+'x')
#         clear_mfile(matdir,name+'y')
    return x,y,dict

def cal_xy4_settings(total_time,electron_pi_len,nSlots,optimization_channels,Sx,Sy,power_limit=1.0):
    num_pi_pulse = 4.0
    total_slots = nSlots+num_pi_pulse
    tau = (total_time - electron_pi_len*num_pi_pulse)/num_pi_pulse/2
    slices = int(nSlots/num_pi_pulse/2.0)
    H_eslices = []
    hempty = Qobj(zeros(Sx.shape),dims = Sx.dims,shape=Sx.shape)
    pi_positions = [slices,slices*3+1,slices*5+2,slices*7+3]
    tlist=linspace(0,tau,slices+1).tolist()
    H_eslices.extend([hempty] * (slices))
    H_eslices.append(Sx)
    tlist.extend(linspace(tlist[-1]+electron_pi_len,tlist[-1]+electron_pi_len+tau*2,slices*2+1).tolist())
    H_eslices.extend([hempty] * (slices*2))
    H_eslices.append(Sy)  
    tlist.extend(linspace(tlist[-1]+electron_pi_len,tlist[-1]+electron_pi_len+tau*2,slices*2+1).tolist())
    H_eslices.extend([hempty] * (slices*2))
    H_eslices.append(Sy)  
    tlist.extend(linspace(tlist[-1]+electron_pi_len,tlist[-1]+electron_pi_len+tau*2,slices*2+1).tolist())
    H_eslices.extend([hempty] * (slices*2))
    H_eslices.append(Sx)  
    tlist.extend(linspace(tlist[-1]+electron_pi_len,tlist[-1]+electron_pi_len+tau,slices+1).tolist())
    H_eslices.extend([hempty] * (slices+1))
    print 'tlist len:',len(tlist)
    tslices = []
    for n in range(1,len(tlist)):
        tslices.append(tlist[n] - tlist[n-1])
    
    control_mask = []
    init_controls = []
    random.seed()
    random.uniform(-abs(power_limit),abs(power_limit))
    
    for i in range(len(tslices)):
        tmpm = []
        tmpc = []
        for j in range(len(optimization_channels)*2):
            if i in pi_positions:
                tmpm.append(0.0)
                tmpc.append(0.0)
            else:   
                tmpm.append(1.0)
#                 tmpc.append(random.random())
                tmpc.append(0.0)
            
        control_mask.append(tmpm)
        init_controls.append(tmpc)
    
#     for i in range(optimization_channels*2):  
#         control_mask[slices][i] = 0
#         control_mask[slices][i] = 0
#         control_mask[slices][i] = 0
#         control_mask[slices][i] = 0
#         
        
    
    wait_pulse = single_channel_pulse(MW,tau,0.0)
    wait_pulse2 = single_channel_pulse(MW,tau*2,0.0)
    xpi = single_channel_pulse(MW,electron_pi_len,0.0,0)
    ypi = single_channel_pulse(MW,electron_pi_len,0.0,pi/2)
    electron_pulse = pulse_sequence([wait_pulse,xpi,wait_pulse2,ypi,wait_pulse2,xpi,wait_pulse2,ypi,wait_pulse])
    
    return tau,tlist,tslices,init_controls,control_mask,electron_pulse,H_eslices

def test_opt_pulse_N14(target_fide):
    coupling = [-2.2e-3,-2.2e-3,-2.2e-3]
    B0=[0,0,0.5]
    RF_pow = 1.5e-6
    MW_pow = 1e-2 
    electron_pi_len = 0.5/(MW_pow * sqrt(2))   
    
    total_time = abs(40*1.0/coupling[2])
    nSlots = 400.0
    
#     slice_width = total_time/nSlots
#     nSlots = ceil(total_time/slice_width)
#     total_time = slice_width * nSlots
    
    
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

    Ix = tensor(basis(3,1)*basis(3,1).dag()+basis(3,2)*basis(3,2).dag(),basis(3,0)*basis(3,1).dag()+basis(3,1)*basis(3,0).dag())
    targetU = (-1j*Ix*pi/2).expm()
    print_matrix(targetU,'target U')
 
    Sx01 = MW_pow*pi*2*\
        tensor(basis(3,1)*basis(3,2).dag()+basis(3,2)*basis(3,1).dag(),basis(3,1)*basis(3,1).dag()+basis(3,0)*basis(3,0).dag())\
        /sqrt(2)
    Sy01 = MW_pow*pi*2j*\
        tensor(basis(3,1)*basis(3,2).dag()-basis(3,2)*basis(3,1).dag(),basis(3,1)*basis(3,1).dag()+basis(3,0)*basis(3,0).dag())\
        /sqrt(2)
    print_matrix(Sy01,'Sy01=')
 
#     pix = (-1j*Ix*pi/2).expm()
#     print_matrix(targetU.full(),'target U')
       
    rframe = rotating_frame(nvsys,H0)
    Hleft = nvsys.get_static_ham(B0)-rframe.H0
    s_coup,s_mw = 300,300

    
    max_freq_sim,args_sim,opts_sim,Ht_sim = gen_rotating_frame_setting(nvsys,rframe,Hleft,sim_channels,s_coup,s_mw,debug=True)
    
    tau,tlist,tslices,init_controls,control_mask,electron_pulse,H_eslices=cal_xy4_settings(total_time,electron_pi_len,nSlots,optimization_channels,Sx01,Sy01)
    empty_pulse0 = single_channel_pulse(RF0,total_time,0.0,0.0)
    empty_pulse1 = single_channel_pulse(RF1,total_time,0.0,0.0)
    test_phase_pulse = multichannel_pulse([empty_pulse0,empty_pulse1,electron_pulse])
#     total_pulse = multichannel_pulse([empty_pulse0,empty_pulse1,])
    args_sim['PULS']=test_phase_pulse      
    print len(tslices)
    U_phase = cal_propagator(args_sim,Ht_sim,tlist,1)
    print_matrix(U_phase,'U test phase:')
 
    max_freq_opt,args_opt,opts_opt,Ht_opt = gen_rotating_frame_setting(nvsys,rframe,Hleft,optimization_channels,s_coup,s_mw,debug=True)
       
    H_static, ham_static_nondiagnal, chamsx, chamsy  = gen_ham_series_for_matlab(nvsys,rframe,Hleft,optimization_channels,s_coup,s_mw,tlist,P)
 
    Ht_opt[0] = P*Ht_opt[0]*P
    for i in range(1,len(Ht_opt)):
        Ht_opt[i][0] = P*Ht_opt[i][0]*P
     
       
    x,y,dict=optimize_operator_DD(targetU,target_fide,H_static, ham_static_nondiagnal,H_eslices, chamsx, chamsy,P,init_controls,control_mask,tslices)
    grape_pulse = sliced_pulse([[RF0,x[RF0],y[RF0]],[RF1,x[RF1],y[RF1]]],tslices,tlist)
   
    print 'in rotating frame'
    total_pulse = multichannel_pulse([grape_pulse,electron_pulse])
    args_sim['PULS']=total_pulse      
            
    U_rot = cal_propagator(args_sim,Ht_sim,tlist,1)
#     U_rot = propagator(Ht,total_time,[],args,opts)
    print_matrix(U_rot.full(),'U=')
    fider = real((U_rot.dag()*targetU).tr()/9.0)
    print 'fidelity by tr is ',fider
   
   
    print 'in lab frame'    
    args,opts,Ht = gen_lab_frame_setting(nvsys.get_static_ham(B0) ,nvsys.cal_global_control_matrix(),sim_channels) 
    args['PULS']=total_pulse
#     U_lab = cal_propagator_deprecated(args,Ht,total_time,total_time*30,num_cpus=6)
    U_lab = cal_propagator_deprecated(args,Ht,total_time,nSlots,num_cpus=6)
    U_h0 = (-1j*rframe.H0*total_time*pi*2).expm()
    evo_states = [ket2dm(state10),ket2dm(state11),ket2dm(state00),ket2dm(state01)]
    expect_states = [U_h0*targetU*ket2dm(state10)*targetU.dag()*U_h0.dag(),\
                     U_h0*targetU*ket2dm(state11)*targetU.dag()*U_h0.dag(),\
                     U_h0*targetU*ket2dm(state00)*targetU.dag()*U_h0.dag(),\
                     U_h0*targetU*ket2dm(state01)*targetU.dag()*U_h0.dag()]
    print_matrix(U_lab.full(),'U_lab=')
    fidel = cal_fide_by_state(U_lab,evo_states,expect_states)
    expect_states = [targetU*ket2dm(state10)*targetU.dag(),\
                     targetU*ket2dm(state11)*targetU.dag(),\
                     targetU*ket2dm(state00)*targetU.dag(),\
                     targetU*ket2dm(state01)*targetU.dag()]
    fider = cal_fide_by_state(U_rot,evo_states,expect_states)
    qsave(U_lab,'U_c13_lab')
    msg = 'fidelity @ power '+str(RF_pow) + ' (r/l):' +str(fider) +','+str(fidel)
    print msg
     
if __name__=='__main__':
    test_opt_pulse_N14(0.99)
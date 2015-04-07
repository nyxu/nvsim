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
# coupling = 8.2e-3
def write_mfile(dir,funcname,funclines):
    fhandle = open(dir+funcname+'.m','w')
    fhandle.write(funclines)
    fhandle.close()
    
def clear_mfile(dir,funcname):
    os.remove(dir+funcname+'.m')

def diff_hams(Ht,args,total_time,nslots,Hmats,Umats):
    dt = float(total_time)/nslots
#     print len(Hmats),len(Hmats[0])
    print 'starting diff ...'
    unequal_num = 0
    if len(Hmats)!=nslots:
        Hmats = Hmats[0]
        Umats = Umats[0]
    for t in range(int(nslots)):
        Hth = Ht[0]
        for i in range(1,len(Ht)):
            Hth = Hth + Ht[i][0]*Ht[i][1](t*dt,args)
        
        Uth = (-1j*Hth*dt).expm()
        Hmat = Qobj(Hmats[t],dims=Ht[0].dims,shape=Ht[0].shape)
        Umat = Qobj(Umats[t],dims=Ht[0].dims,shape=Ht[0].shape)
 
        Hdifference = (Hth - Hmat).tidyup()
        Udifference = (Uth - Umat).tidyup()
        if(Hdifference.full().any()):
            print t,': hamiltonian unequal'
            unequal_num = unequal_num + 1
        if(Udifference.full().any()):
            print t,': Operator unequal'
            unequal_num = unequal_num + 1
    print 'complete diff'
    return (unequal_num==0)
        
def optimize_operator(targetU,target_fide,hamsettings,total_time,nSlots,args,P1,max_power=1.0):
    
    H_static,static_nondiagnal_hams,channel_ctl_hams = hamsettings
    matlabcmd = '/Applications/MATLAB_R2012a.app/bin/matlab -nodesktop -nosplash -nojvm -wait -r pulse_optimization  '
    matdir = '/Users/nyxu/Dropbox/workspace/MATLAB/custemed_DYNAMO_v_1_2/'
    opfile = matdir + 'pulse_optimization.mat'
    pulsefile =  matdir + 'optimized_pulse.mat'
    chamfuncs = []
    nondiaghamfuncs=[]
    for name in static_nondiagnal_hams:
        write_mfile(matdir,name,static_nondiagnal_hams[name])
        nondiaghamfuncs.append(name)
    for item in channel_ctl_hams:
        name,cx,cy = item
        write_mfile(matdir,name+'x',cx)
        write_mfile(matdir,name+'y',cy)
        chamfuncs.append(name+'x')
        chamfuncs.append(name+'y')
        
    dhams = {}
    dhams['ctl_funcs'] = chamfuncs
    dhams['nd_funcs'] = nondiaghamfuncs
    dhams['ext_mode']=1
    dhams['ndims']=targetU.shape[0]
#     for m in matrix:
#         dhams[m]=matrix[m].full()
    dhams['Utar']=targetU.full()
    dhams['target_fide']=target_fide
    dhams['H_static']=H_static.full()
    dhams['total_time']=total_time
    dhams['nSlots']=nSlots
    dhams['max_power'] = max_power
#     print_matrix(Ht[4][0].full())
    sio.savemat(opfile,dhams)
    ret = os.system(matlabcmd)
    x = {}
    y = {}
    if ret == 0 and os.path.isfile(pulsefile):
        dict = sio.loadmat(pulsefile)
#         print dict.keys()
        controls = dict['controls']        
        for k in range(len(channel_ctl_hams)):
            name=channel_ctl_hams[k][0]
            x[name] = []
            y[name] = []
            for i in range(len(controls)):
                x[name].append(controls[i][k*2])
                y[name].append(controls[i][k*2+1])
            
    for i in range(len(channel_ctl_hams)):
        clear_mfile(matdir,channel_ctl_hams[i][0]+'x')
        clear_mfile(matdir,channel_ctl_hams[i][0]+'y')
    return x,y,dict
    


def test_opt_pulse_c13(target_fide):
    coupling = [[5e-3,-6.3e-3,-2.9e-3],
               [-6.3e-3,4.2e-3,-2.3e-3],
               [-2.9e-3,-2.3e-3,8.2e-3]]

    nvsys = nv_system([nv_electron_spin(),C13_nuclear_spin()],[[0,1,coupling]]) 
    state00 = nvsys.get_basis([0,-0.5])
    state01 = nvsys.get_basis([0,0.5])
    state10 = nvsys.get_basis([1,-0.5])
    state11 = nvsys.get_basis([1,0.5])
    B0=[0,0,0.05]
    H0 = nvsys.get_local_static_ham(B0[2]) + nvsys.get_single_coupled_ham(0,1,'zz')
    
    c13freq = abs(get_transition_freq(H0,state10,state11))
    print 'frequency of RF is ',c13freq
    RF_pow = 1e-5
    rf1chan = MWChannel(c13freq,0,abs(RF_pow/Cgn),nvsys.cal_nulclear_eff_control_matrix(B0[2]))
    channels={'RF':rf1chan} 
    

    P1 = nvsys.cal_electron_spin_projector(1)
    Ix = tensor(basis(3,0)*basis(3,0).dag(),sigmax())
    targetU = (-1j*Ix*pi/2).expm()
    
    rframe = rotating_frame(nvsys,H0)
    Hleft = nvsys.get_static_ham(B0)-rframe.H0
    s_coup,s_mw = 30,0.1
    max_freq,args,opts,Ht = gen_rotating_frame_setting(nvsys,rframe,Hleft,channels,s_coup,s_mw,debug=False)
    print nvsys.cal_eff_g(1,1,B0[2])

    total_time = 10000.0
    nSlots = 1000.0

    H_static,static_nondiagnal_hams,channel_ctl_hams  = gen_ham_for_matlab(nvsys,rframe,Hleft,channels,s_coup,s_mw,P1)
#     print static_nondiagnal_hams.values()[0]
#     print '======================'
#     for chname in channel_ctl_hams:
#         print channel_ctl_hams[chname][0]
#         print '-----------------------'
#         print channel_ctl_hams[chname][1]
#         print '======================'
    Ht[0] = P1*Ht[0]*P1
    for i in range(1,len(Ht)):
        Ht[i][0] = P1*Ht[i][0]*P1
    hamsetings = [H_static,static_nondiagnal_hams,channel_ctl_hams]
    x,y,dict=optimize_operator(targetU,target_fide,hamsetings,total_time,nSlots,args,P1)
    print 'in rotating frame'
    grape_pulse = shaped_pulse(['RF',x['RF'],y['RF']],total_time/nSlots)
    args['PULS']=grape_pulse
    
#     if you want to diff the hamiltonians here and matlab,please uncomment the following
#     diff_hams(Ht,args,total_time,nSlots,dict['H'],dict['U'])
       
        
    U_rot = cal_propagator(args,Ht,total_time,nSlots,1)
    print_matrix(U_rot.full(),'U=')
    fider = real((U_rot.dag()*targetU).tr()/6.0)
    print 'fidelity by tr is ',fider


    print 'in lab frame'    
    args,opts,Ht = gen_lab_frame_setting(nvsys.get_static_ham(B0) ,nvsys.cal_global_control_matrix(),channels) 
    args['PULS']=grape_pulse
    U_lab = cal_propagator(args,Ht,total_time,total_time*30)
    U_h0 = (-1j*rframe.H0*total_time*pi*2).expm()
    evo_states = [ket2dm(state10),ket2dm(state11)]
    expect_states = [U_h0*targetU*ket2dm(state10)*targetU.dag()*U_h0.dag(),\
                     U_h0*targetU*ket2dm(state11)*targetU.dag()*U_h0.dag()]
    print_matrix(U_lab.full(),'U_lab=')
    fidel = cal_fide_by_state(U_lab,evo_states,expect_states)
    expect_states = [targetU*ket2dm(state10)*targetU.dag(),\
                     targetU*ket2dm(state11)*targetU.dag()]
    fider = cal_fide_by_state(U_rot,evo_states,expect_states)
    qsave(U_lab,'U_c13_lab')
    msg = 'fidelity @ power '+str(RF_pow) + ' (r/l):' +str(fider) +','+str(fidel)
    print msg
#     send_result(['nyxu@ustc.edu.cn'],msg,'nyxu@mail.ustc.edu.cn', 'trapdoor') 
    
def test_opt_pulse_N14(target_fide):
    coupling = [-2.2e-3,-2.2e-3,-2.2e-3]

    nvsys = nv_system([nv_electron_spin(),N14_nuclear_spin()],[[0,1,coupling]]) 
    e_state0 = 0
    e_state1 = -1
    state00 = nvsys.get_basis([e_state0,0])
    state01 = nvsys.get_basis([e_state0,1])
    state10 = nvsys.get_basis([e_state1,0])
    state11 = nvsys.get_basis([e_state1,1])
    P = nvsys.cal_electron_spin_projector(e_state0)+nvsys.cal_electron_spin_projector(e_state1)
    print_matrix( P.full(),'P=')
    B0=[0,0,0.5]
    print nvsys.cal_eff_g(1,0,B0[2])
    print nvsys.cal_eff_g(1,-1,B0[2])
    H0 = nvsys.get_local_static_ham(B0[2]) + nvsys.get_single_coupled_ham(0,1,'zz')
#     print_matrix(H0.full())
    N14freq0 = abs(get_transition_freq(H0,state00,state01))
    N14freq1 = abs(get_transition_freq(H0,state10,state11))
    print 'frequency of RF is (MHz)',N14freq0*1e3,N14freq1*1e3
    RF_pow = 2e-6
    RF0 = 'RF0'
    RF1 = 'RF1'
    rf0chan = MWChannel(N14freq0,0,abs(RF_pow/Ngn),nvsys.cal_nulclear_eff_control_matrix(B0[2]))
    rf1chan = MWChannel(N14freq1,0,abs(RF_pow/Ngn),nvsys.cal_nulclear_eff_control_matrix(B0[2]))
    channels={RF0:rf0chan,RF1:rf1chan} 
#     channels={RF0:rf0chan} 
    

    Ix = tensor(basis(3,1)*basis(3,1).dag()+basis(3,2)*basis(3,2).dag(),basis(3,0)*basis(3,1).dag()+basis(3,1)*basis(3,0).dag())
    targetU = (-1j*Ix*pi/2).expm()
    print_matrix(targetU.full(),'target U')
    rframe = rotating_frame(nvsys,H0)
    Hleft = nvsys.get_static_ham(B0)-rframe.H0
    s_coup,s_mw = 300,100
    max_freq,args,opts,Ht = gen_rotating_frame_setting(nvsys,rframe,Hleft,channels,s_coup,s_mw,debug=True)
#     print nvsys.cal_eff_g(1,1)

    total_time = 15000.0
    nSlots = 300.0

    H_static,static_nondiagnal_hams,channel_ctl_hams  = gen_ham_for_matlab(nvsys,rframe,Hleft,channels,s_coup,s_mw,P)
    print static_nondiagnal_hams.values()[0]
    print '======================'
    for item in channel_ctl_hams:
        print item[0],':'
        print item[1]
        print '-----------------------'
        print item[2]
        print '======================'
    Ht[0] = P*Ht[0]*P
    for i in range(1,len(Ht)):
        Ht[i][0] = P*Ht[i][0]*P
    hamsetings = [H_static,static_nondiagnal_hams,channel_ctl_hams]
    x,y,dict=optimize_operator(targetU,target_fide,hamsetings,total_time,nSlots,args,P)
    print 'in rotating frame'
#     grape_pulse0 = shaped_pulse([RF0,x[RF0],y[RF0]],total_time/nSlots)
    grape_pulse = shaped_pulse([[RF0,x[RF0],y[RF0]],[RF1,x[RF1],y[RF1]]],total_time/nSlots)
    args['PULS']=grape_pulse
    
#     if you want to diff the hamiltonians here and matlab,please uncomment the following
#     diff_hams(Ht,args,total_time,nSlots,dict['H'],dict['U'])
       
        
    U_rot = cal_propagator(args,Ht,total_time,nSlots,1)
    print_matrix(U_rot.full(),'U=')
    fider = real((U_rot.dag()*targetU).tr()/9.0)
    print 'fidelity by tr is ',fider


    print 'in lab frame'    
    args,opts,Ht = gen_lab_frame_setting(nvsys.get_static_ham(B0) ,nvsys.cal_global_control_matrix(),channels) 
    args['PULS']=grape_pulse
    U_lab = cal_propagator(args,Ht,total_time,total_time*30)
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
#     send_result(['nyxu@ustc.edu.cn'],msg,'nyxu@mail.ustc.edu.cn', 'trapdoor') 
     
if __name__=='__main__':
    test_opt_pulse_N14(0.99)
from qcexp import *
from nv_system import *
from pulse_control import *
from timeutils import *
from nvutils import *
import scipy.io as sio
import os
from numpy import *
coupling = [[5e-3,-6.3e-3,-2.9e-3],
           [-6.3e-3,4.2e-3,-2.3e-3],
           [-2.9e-3,-2.3e-3,8.2e-3]]
# coupling = dot(coupling,0.125).tolist()
from scipy import fftpack
from numpy import *
from msgutils import *

# coupling = 8.2e-3
def optimize_C13_operator(targetU,Ht,total_time,nSlots,args,P1):
    matlabcmd = '/Applications/MATLAB_R2012a.app/bin/matlab -nodesktop -nosplash -nojvm -wait -r optimize_c13  '
    matdir = '/Users/nyxu/Dropbox/workspace/MATLAB/custemed_DYNAMO_v_1_2/'
    opfile = matdir + 'c13_opt_op.mat'
    pulsefile =  matdir + 'optimized_c13_pulse.mat'
    dhams = {}
    dhams['Utar']=targetU.full()
    
    Ix = tensor(basis(3,0)*basis(3,0).dag(),sigmax())
    Iy = tensor(basis(3,0)*basis(3,0).dag(),sigmay())
    c = Ht[-1][0].matrix_element(tensor(basis(3,0),basis(2,0)).dag(),tensor(basis(3,0),basis(2,1)))
    a = float(real(c))
    b = float(imag(c))
    print a,b
    max_power = sqrt((a**2+b**2)/2.0)
    dhams['control_ham1']=(Ix).full()
    dhams['control_ham2']=(Iy).full()
    
    Hdrift = []
    dt = float(total_time)/nSlots
    
    for i in range(int(nSlots)):
        if i % 100 == 0 :print i
        h = Ht[0]
        for j in range(2,len(Ht)-2):     
            h = h +  Ht[j][0]*Ht[j][1](dt*i,args)
#         Hdrift.append(zeros((6,6)))#(P1*h*P1).full())
#         print 
        Hdrift.append((P1*h*P1).full())
    dhams['Hdrift']=Hdrift
    dhams['total_time']=total_time
    dhams['nSlots']=nSlots
    dhams['max_power'] = max_power
#     print_matrix(Ht[4][0].full())
    sio.savemat(opfile,dhams)
    ret = os.system(matlabcmd)
    if ret == 0 and os.path.isfile(pulsefile):
        dict = sio.loadmat(pulsefile)
        x1,y1 = dict['x'],dict['y']
        x = []
        y = []
        normab = a**2+b**2
        for i in range(len(x1)):
            x.append((a*x1[i][0]-b*y1[i][0])/normab)
            y.append((b*x1[i][0]+a*y1[i][0])/normab)
        return x,y
def cal_fide(U,Utar,P1):
    U=(P1*U*P1)
    U = U - U.tr()/2.0
    U=(P1*U*P1)
    fide = abs((U.dag()*(P1*Utar*P1)).tr())/2.0
#     print 'fidelity by Propagator is ',fide
#     print_matrix(U.full(),'U=')
#     print_matrix (targetU.full(),'targetU=')
#     fide = abs((P1*(U.dag()*Utar)*P1).tr())/2.0
    return fide

from scipy.linalg import *

def cal_fide_by_state(U,evo_states,ideal_states):
#     fidelities = []
#     for i in range(len(evo_states)):
#         output = U * evo_states[i] * U.dag()
#         fidelities.append((ideal_states[i] * output * ideal_states[i].dag()).tr())
#         print 'fidelity ',fidelities[-1]
#     return average(fidelities)
#     fidelities = []
#     for i in range(len(evo_states)):
#         rou_exp = U * evo_states[i] * U.dag()
#         rou_ideal_sqrt = Qobj(sqrtm(ideal_states[i].full()),dims=rou_exp.dims,shape=rou_exp.shape)
#         tmp = rou_ideal_sqrt*rou_exp*rou_ideal_sqrt
#         fidelities.append(trace(sqrtm(tmp.full()))**2)
#         print 'fidelity ',fidelities[-1]
    fidelities = []
    for i in range(len(evo_states)):
        output = U * evo_states[i] * U.dag()
        fidelities.append((ideal_states[i] * output * ideal_states[i].dag()).tr())
        print 'fidelity ',fidelities[-1]
    return average(fidelities)


    
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
    rf1chan = MWChannel(c13freq,0,abs(RF_pow/Cgn),nvsys.cal_nulclear_eff_control_matrix(B0[2]))
    channels={'RF':rf1chan} 
    
    P1 = nvsys.cal_electron_spin_projector(1)
    Ix = tensor(basis(3,0)*basis(3,0).dag(),sigmax())
    targetU = (-1j*Ix*pi/2).expm()
    
    rframe = rotating_frame(nvsys,H0)
    Hleft = nvsys.get_static_ham(B0)-rframe.H0
    s_coup,s_mw = 30,0.1
    max_freq,args,opts,Ht = gen_rotating_frame_setting(nvsys,rframe,Hleft,channels,s_coup,s_mw)
    print nvsys.cal_eff_g(1,1,B0[2])
#     Ht[0] = P1*Ht[0]*P1
#     Ht[1][0] = Qobj(zeros((6,6)),dims = Ht[0].dims,shape=Ht[0].shape)
#     for i in range(2,len(Ht)):
#         Ht[i][0] = P1*Ht[i][0]*P1
    total_time = 4000.0
    nSlots = 1000.0
    x,y=optimize_C13_operator(targetU,Ht,total_time,nSlots,args,P1)

    print 'in rotating frame'
    grape_pulse = shaped_pulse(['RF',x,y],total_time/nSlots)
    args['PULS']=grape_pulse
    U_rot = cal_propagator(args,Ht,total_time,nSlots,1)
    print_matrix(U_rot.full(),'U=')
    fider = cal_fide(U_rot,targetU,P1)

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
    qsave(U_lab,'U_c13_lab')
    msg = 'fidelity @ power '+str(RF_pow) + ' (r/l):' +str(fider) +','+str(fidel)
    print msg
    send_result(['nyxu@ustc.edu.cn'],msg,'nyxu@mail.ustc.edu.cn', 'trapdoor') 

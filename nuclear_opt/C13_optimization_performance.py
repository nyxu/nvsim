from qcexp import *
from nv_system import *
from pulse_control import *
from timeutils import *
from nvutils import *
from C13_optimize import *
default_NV_C13_coupling = [
                           [5e-3,-6.3e-3,-2.9e-3],
                           [-6.3e-3,4.2e-3,-2.3e-3],
                           [-2.9e-3,-2.3e-3,8.2e-3]]
global debug
if __name__=='__main__':
    
    nvsys = nv_system([nv_electron_spin(),C13_nuclear_spin()],[[0,1,default_NV_C13_coupling]]) 
    
    state00 = nvsys.get_basis([0,-0.5])
    state01 = nvsys.get_basis([0,0.5])
    state10 = nvsys.get_basis([1,-0.5])
    state11 = nvsys.get_basis([1,0.5])
    B0=[0,0,0.05]
    
    P1=nvsys.cal_electron_spin_projector(ms=1)
    H0 = nvsys.get_local_static_ham(B0[2]) + nvsys.get_single_coupled_ham(0,1,'zz')

    c13freq = abs(get_transition_freq(H0,state10,state11))
    print 'frequency of RF is ',c13freq
    RF_pow = 1e-4
    rf1chan = MWChannel(c13freq,0,abs(RF_pow/Cgn),nvsys.cal_nulclear_eff_control_matrix())
    channels={'RF':rf1chan}
    rframe = rotating_frame(nvsys,H0)
    Hleft = nvsys.get_static_ham(B0)-rframe.H0
    evo_time = 2500 
    s_coup = 30
    s_mw = 0.1
    max_freq,args,opts,Ht = gen_rotating_frame_setting(nvsys,rframe,Hleft,channels,s_coup,s_mw)
    pcontrl = pulsecontrol(single_channel_pulse('RF'))

    
    rou0 = ket2dm(state10)
    Ix = tensor(basis(3,0)*basis(3,0).dag(),sigmax())
    tlist = linspace(100,evo_time,40)
    result = []
    Us = []
    gpulses = []
    for t in tlist:
        targetU = (-1j*Ix*4e-4*t*pi).expm()
        total_time = max(2000.0,t*2)
        nSlots = max(10.0,floor(total_time/20))
        x,y=optimize_C13_operator(targetU,Ht,total_time,nSlots,args,P1)
        print x
        grape_pulse = shaped_pulse(['RF',x,y],total_time/nSlots)
        args['PULS']=grape_pulse
         
        U = cal_propagator(args,Ht,total_time,nSlots)
        overlap = real((U*rou0*U.dag()*rou0).tr())
        print 'overlap @time ',t,' is ',overlap
        result.append(overlap)
        
        Us.append(U)
        gpulses.append([x,y,total_time/nSlots])
        qsave(Us,'Us')
        qsave(gpulses,'gpulses')
        
    title('C13 Rabi OSC using GRAPE Method')
    xlabel('time(ns)')
    ylabel('overlap @ |1,1/2>')
    plot(tlist,real(result),'r.',label='GRAPE')
    tlist = linspace(0, evo_time, 100)
    targets=[]
    for t in tlist:
        targetU = (-1j*Ix*4e-4*t*pi).expm()
        targets.append(real((targetU*rou0*targetU.dag()*rou0).tr()))
    plot(tlist,real(targets),label='Theoritical')
    legend()
    show()

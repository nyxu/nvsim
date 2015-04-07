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

def cal_electron_propagator(static_ham,ctl_hams,channels,evo_time):
    s_coup,s_mw = 100,1
    pcontrl = pulsecontrol(single_channel_pulse('MW'))
    args,opts,Ht = gen_lab_frame_setting(static_ham ,ctl_hams,channels)
    args['PULS']=single_channel_pulse('MW')
    prop=cal_propagator(args,Ht, [0,evo_time], [], args,opts)
    return prop
    
def cal_eff_pow(static_ham,ctl_hams,init_state,channels,evo_time = 5000):
    s_coup,s_mw = 100,1
    pcontrl = pulsecontrol(single_channel_pulse('MW'))
    args,opts,Ht = gen_lab_frame_setting(static_ham ,ctl_hams,channels)
    sampling_times = 1 * evo_time
    tlist = linspace(0, evo_time, sampling_times)
    rou0 = (ket2dm(init_state) if init_state.type=='ket' else init_state)
    args['PULS']=single_channel_pulse('MW')
    output=mesolve(Ht, rou0, tlist, [], rou0, args,opts)
    plot(tlist,abs(output.expect[0]))
    show()
    deltat = tlist[1]-tlist[0]
    freqs = fftpack.fftfreq(len(tlist),deltat)
    values = output.expect[0]     
    values = values - average(values)        
  
    fft = abs(fftpack.fft(values))

#     print 'fft:',fft
    maxfft = 0
    maxfftpos = -2
    minfft = 1000
    minfftpos = -2
    for i  in range(len(fft)):
        f = fft[i]
        if f > maxfft:
            maxfft = f
            maxfftpos = i
        elif f < minfft:
            minfft = f
            minfftpos = i
    
    if abs(maxfft) > abs(minfft):
        eff_pow = abs(freqs[maxfftpos])
    else:
        eff_pow = abs(freqs[minfftpos])
    print 'effective pow:',eff_pow
    plot(freqs,fft)
    show()    
    return eff_pow

if __name__=='__main__':
    coupling = [-2.2e-3,-2.2e-3,-2.2e-3]

    nvsys = nv_system([nv_electron_spin(),N14_nuclear_spin()],[[0,1,coupling]]) 
    e_state0 = 0
    e_state1 = -1
    state00 = nvsys.get_basis([e_state0,0])
    state01 = nvsys.get_basis([e_state0,1])
    state10 = nvsys.get_basis([e_state1,0])
    state11 = nvsys.get_basis([e_state1,1])

    B0=[0,0,0.5]

    H0 = nvsys.get_local_static_ham(B0[2],spins=[0])
    electron_ns_freq = abs(get_transition_freq(H0,state00,state10))

    print 'frequency of electron non-selective pulse is (MHz)',electron_ns_freq
    MW_pow = 2e-2
    MW = 'MW'

    mw0chan = MWChannel(electron_ns_freq,0,abs(MW_pow/ge),nvsys.cal_electron_eff_control_matrix())
    channels={MW:mw0chan} 

    init_state = tensor(ket2dm(basis(3,1)),qeye(3))

#     pow = cal_eff_pow(nvsys.get_static_ham(B0),nvsys.cal_global_control_matrix(),init_state,channels, 1000)
#     print 'pow is ',pow
    U = cal_electron_propagator(nvsys.get_static_ham(B0),nvsys.cal_global_control_matrix(),channels, 0.5/(MW_pow*sqrt(2)))
    print_matrix(U,'U=')
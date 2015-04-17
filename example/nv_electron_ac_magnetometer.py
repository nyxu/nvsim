#
# Sensing of AC magnetic field with NV centers, scheme from Nature 455,644 (2008)
#
from qutip import *
from nv_system import *
from pulse_control import *
from rotating_frame import *
from scipy import fftpack    
from qcexp import *  
from threadreporter import *

import os

def electron_magnetometer(nvsys,evo_ac_fields,points,acfreq,B0,mwfreq,B1,T1,T2,ncpus=1,reportfile = 'magnetometer_ac_process.txt'):
    if os.path.exists(reportfile):
        os.remove(reportfile)
    H0 = nvsys.get_static_ham(B0)
    control_matrix = nvsys.cal_global_control_matrix()
    c_op_list = []
    if T1 and T2:
        c_op_list.extend(nvsys.get_lindblad_op(0, T1, T2))

    rou0 = nvsys.create_matrix(0,ket2dm(basis(3,1)))     # initial state 
    ob_op = rou0*3.0/rou0.shape[0]
    
    
    
    mw_channel = MWChannel(mwfreq,0.0,B1)
    
    evo_list = linspace(0,evo_ac_fields,points).tolist()[1:]
    
    half_pi = simple_pulse(1.0/B1/ge/sqrt(2)/4.0,['MW','AC'])
    pi_pulse = simple_pulse(1.0/B1/ge/sqrt(2)/2.0,['MW','AC'])

    #print 1.0/B1/ge/sqrt(2)/4.0
    freetime = (0.5/acfreq)-pi_pulse.duration
    print freetime
    pulse_slices = [half_pi,simple_pulse(freetime,['AC']),pi_pulse,simple_pulse(freetime,['AC']),half_pi]
    ret_list = []
    if ncpus > 1:
        paras = []
        for acfield in evo_list:
            ac_rf_channel = MWChannel(acfreq,pi/2,acfield,pi/2)
#             ac_rf_channel = ArbWave('cos('+str(acfreq*pi*2)+'*(t+0.0)+'+str(pi/2)+')',acfield*2,pi/2)
            mwchannels = {'MW':mw_channel,'AC':ac_rf_channel}
            paras.append(['Test AC field @ '+str(acfield)+'T',H0,control_matrix,B0,rou0,mwchannels,pulse_slices,c_op_list,acfield,ob_op,reportfile])
        ret_list = parfor(simulate_lab_frame_exp_thread,paras,num_cpus=ncpus)
    else:
        for acfield in evo_list:
#             ac_rf_channel = MWChannel(acfreq,pi/2,acfield,pi/2)
            ac_rf_channel = ArbWave('cos('+str(acfreq*pi*2)+'*(t+offset)+'+str(pi/2)+')',acfield*2,pi/2)
            mwchannels = {'MW':mw_channel,'AC':ac_rf_channel}
#             pulse_slices = [half_pi,simple_pulse(freetime,['AC']),pi_pulse,simple_pulse(freetime,['AC']),half_pi]
            expres = simulate_lab_frame_exp_thread(['Test AC field @ '+str(acfield)+'T',H0,control_matrix,B0,rou0,mwchannels,pulse_slices,c_op_list,acfield,ob_op,reportfile])    
            ret_list.append(expres)
#     result_list = []    
#     for state in ret_list:
#         result_list.append(abs((state*ob_op).tr()))
    return evo_list,ret_list
 
if __name__=='__main__':
#     nvsys = nv_system([nv_electron_spin(),N14_nuclear_spin()],[[0,1,default_NV_N14_coupling]]) 

    nvsys = nv_system([nv_electron_spin(),N14_nuclear_spin()],[[0,1,default_NV_N14_coupling]]) 
    
        
    evo_ac_fields=1.0e-3/300.0
    acfreq = 10e-5 #GHz
    B0=5e-3
    mwfreq = nv_spliting - B0*ge 
    points=12
    B1=1e-3
    T1=1e6
    T2=2e5
    ncpus=3
    evo_list,result_list=electron_magnetometer(nvsys,evo_ac_fields,points,acfreq,B0,mwfreq,B1,T1,T2,ncpus)

    # Plot the result
    plot(evo_list, real(result_list))
    xlabel('Free Evolution Time')
    ylabel('Occupation probability')
    title('AC Sensing of Singnal @ '+str(acfreq*1e9)+'Hz')  

    
         
    show()

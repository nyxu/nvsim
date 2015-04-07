#
# Rabi oscillations of qubit subject to a classical driving field.
#
from qutip import *
from nv_system import *
from pulse_control import *
from rotating_frame import *
from scipy import fftpack    
from qcexp import *  


def electron_fid(nvsys,evo_time,points,B0,mwfreq,B1,T1=0,T2=0,ncpus=1,reportfile = 'fid_process.txt'):
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
    mwchannels = {'MW':mw_channel}
    
    evo_list = linspace(0,evo_time,points).tolist()[1:]
    
    half_pi = simple_pulse(1.0/B1/ge/sqrt(2)/4.0,['MW'])
    #print 1.0/B1/ge/sqrt(2)/4.0
    
    ret_list = []
    if ncpus > 1:
        paras = []
        for time in evo_list:
            pulse_slices = [half_pi,simple_pulse(time,[]),half_pi]
            paras.append(['FID @ '+str(time)+'ns',H0,control_matrix,B0,rou0,mwchannels,pulse_slices,c_op_list,time,ob_op,reportfile])
        ret_list = parfor(simulate_lab_frame_exp_thread,paras,num_cpus=ncpus)
    else:
        for time in evo_list:
            pulse_slices = [half_pi,simple_pulse(time,[]),half_pi]
            state = simulate_lab_frame_exp_thread(['FID @ '+str(time)+'ns',H0,control_matrix,B0,rou0,mwchannels,pulse_slices,c_op_list,time,ob_op,reportfile])    
            ret_list.append(state)
    result_list = []    
    for state in ret_list:
        result_list.append(abs((state*ob_op).tr()))
    return evo_list,result_list
 
if __name__=='__main__':
    nvsys = nv_system([nv_electron_spin(),N14_nuclear_spin()],[[0,1,default_NV_N14_coupling]]) 
    
        
    evo_time=1000
    B0=5e-3
    mwfreq = nv_spliting + B0*ge 
    points=50
    B1=1e-3
    T1=20000
    T2=2000
    ncpus=2
    evo_list,result_list=electron_fid(nvsys,evo_time,points,B0,mwfreq,B1,T1,T2,ncpus)

    # Plot the result
    subplot(211)
    plot(evo_list, real(result_list))
    xlabel('Rotation Time')
    ylabel('Occupation probability')
    title('FID (B0 @ '+str(B0)+'T)')  
    
    deltat=evo_list[1]-evo_list[0]
    # remove the DC component
    data = result_list
    data = data - average(data) 
    fft= fftpack.fft(data)
    freqs = fftpack.fftfreq(len(data),deltat)
    
    subplot(212)
    plot(freqs,real(fft),'r')
    xlabel('Freq')
    ylabel('Arbitrary Unit')
    
         
    show()

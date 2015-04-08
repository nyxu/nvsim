from qutip import *
from nv_system import *
from numpy import *
from rotating_frame import *
import StringIO

from threading import *
from threadreporter import *            
        
            
def simulate_lab_frame_exp_thread(paras):
    thread_name,H0,control_matrix,Bz,init_state_rou,mwchannels,pslices,c_op_list,param,observ_op,reportfile = paras

    print 'Thread ',thread_name, ' start' 
    import time
    starttime = time.clock()

    retval= simulate_lab_frame_exp_strfunc(H0,control_matrix,Bz,init_state_rou,mwchannels,pslices,c_op_list)
    print 'Thread ',thread_name, "cost time :", str(time.clock()-starttime)
    report(reportfile,[param,real((retval*observ_op).tr())])
    return retval

def simulate_lab_frame_exp_strfunc(H0,control_matrix,Bz,init_state_rou,channels,pulse_slices,c_op_list):
    #the time-independent hamiltonian
    matx,maty,matz = control_matrix
    
    #options for mesolve
    opts=Odeoptions()
    opts.nsteps=1e9
    

#     current_phase = 0.0
    current_time = 0.0
    current_state = init_state_rou
    
    for pslice in pulse_slices:
        Ht = [H0*pi*2]
        
        #x direction fields
        strfuncx,strfuncy,strfuncz = ['','','']
        empty_slice_x = empty_slice_y=empty_slice_z=True
        for ch in channels:
            fieldx = 2*channels[ch].field*pslice.get_amplitude(ch)*cos(channels[ch].theta)*cos(channels[ch].phi)
            fieldy = 2*channels[ch].field*pslice.get_amplitude(ch)*cos(channels[ch].theta)*sin(channels[ch].phi)
            fieldz = 2*channels[ch].field*pslice.get_amplitude(ch)*sin(channels[ch].theta)
            if (abs(fieldx) > 0):
                strfuncx += (str(fieldx) + '*cos('+str(2*pi*channels[ch].frequency)+'*(t+'+str(current_time) + ')+'+str(pslice.get_phase(ch)+channels[ch].phase)+')+')
                empty_slice_x = False
            if (abs(fieldy) > 0):
                strfuncy += (str(fieldy) + '*cos('+str(2*pi*channels[ch].frequency)+'*(t+'+str(current_time) + ')+'+str(pslice.get_phase(ch)+channels[ch].phase)+')+')
                empty_slice_y = False
            if (abs(fieldz) > 0):
                strfuncz += (str(fieldz) + '*cos('+str(2*pi*channels[ch].frequency)+'*(t+'+str(current_time) + ')+'+str(pslice.get_phase(ch)+channels[ch].phase)+')+')
                empty_slice_z = False
        strfuncx += '0'
        strfuncy += '0'
        strfuncz += '0'
        if not empty_slice_x:
            Ht.append([-matx*pi*2,strfuncx])
        if not empty_slice_y:
            Ht.append([-maty*pi*2,strfuncy])
        if not empty_slice_z:
            Ht.append([-matz*pi*2,strfuncz])



        tlist = linspace(0.0, pslice.get_duration(), 2)
        
        output = mesolve(Ht, current_state, tlist, c_op_list, [], {},opts)      
        current_state = output.states[-1]
        current_time = current_time + pslice.get_duration()
        
    return current_state

def simulate_lab_frame_exp_callback(nvsys,Bz,init_state_rou,mwchannels,pulse_slices,T1,T2):
    #the time-independent hamiltonian
    H0 = nvsys.get_static_ham(Bz)
    Ht = [H0*pi*2]
    control_matrix = nvsys.cal_global_control_matrix()
    matx,maty,matz = control_matrix
    Ht.append([-matx*pi*2,Bdx])
#     Ht.append([-maty*pi*2,Bdy])
#     Ht.append([-matz*pi*2,Bdz])
    # define the time-dependence of the hamiltonian using the list-string format
    args = {'MWCH':mwchannels}
    
    # unitary evolution
    c_op_list = []
    if T1 and T2:
        c_op_list.extend(nvsys.get_lindblad_op(0, T1, T2))
    
    #options for mesolve
    opts=Odeoptions()
    opts.nsteps=1e6
    

#     current_phase = 0.0
    current_time = 0.0
    current_state = init_state_rou
    
    import time
    starttime = time.clock()
    for pslice in pulse_slices:
        
        args['ST_TIME'] = current_time
        args['PSLICE'] = pslice
        tlist = linspace(0.0, pslice.get_duration(), 2)
        
        output = mesolve(Ht, current_state, tlist, c_op_list, [], args,opts)      
        current_state = output.states[-1]
        current_time = current_time + pslice.get_duration()
        
    print "cost time :",str(time.clock()-starttime)
    return current_state



def gen_rotating_frame_setting(sysdef,rframe,Hleft,channels,s_coup,s_mw,threshold = 1e-6,debug=True):    
    Ht = []
    rh0 = 0
    max_freq = 0
    diagonal_ham_items = []
    nondiagonal_ham_items = []
    import StringIO
    
    sperturb = Perturbation(Hleft,0,1.0)
    diagonal_item,nondiagonal_items = \
        sperturb.rotating_wave_approximate(rframe.P,rframe.omega,s_coup)               
    rh0 = rh0 + diagonal_item
    
    for nondiagonal_item in nondiagonal_items:
        if(abs(nondiagonal_item[1])<threshold):
            rh0 = rh0+nondiagonal_item[0]
        else:
            nondiagonal_item.append('static_item')
            nondiagonal_ham_items.append(nondiagonal_item)
        
    
    for chname in channels:
        channel = channels[chname]
        mwperturb=MWFieldPerturbation(channel)            
        diagonal_item,nondiagonal_items = \
            mwperturb.rotating_wave_approximate(rframe.P,rframe.omega,s_mw)

        if diagonal_item.full().any():
            diagonal_ham_items.append([diagonal_item,chname])
        
        for nondiagonal_item in nondiagonal_items:
            nondiagonal_item.append(chname)
            nondiagonal_ham_items.append(nondiagonal_item)
    
    if rh0 != 0:
        Ht.append(rh0*pi*2)   
    hamNO=0    
    for item in diagonal_ham_items:
        product,chname = item
        hamfunc='rotation_ham_' +str(hamNO)
        hamNO=hamNO+1
        output = StringIO.StringIO()
        output.write('def '+ hamfunc + '(t,args):\n')
        output.write('\tpulses,channels = args[\'PULS\'],args[\'MWCH\']\n')
        output.write('\tamp,phase = pulses.get_shape(\''+chname+'\',t)\n')
        
        output.write('\treturn cos(('+str(channels[chname].frequency)+\
                     '*t*pi*2+channels[\''+chname+'\'].phase + phase))  * amp\n')
        output.write('Ht.append([product*pi*2*'+str(channels[chname].field*2)+', '+hamfunc+'])\n')
        exec output.getvalue()
        print output.getvalue()
        print_matrix(product,'product=')
        output.close()
        max_freq = max(max_freq,abs(channels[chname].frequency))
        
    for item in nondiagonal_ham_items:
        product,omega_freq,chname = item
        max_freq = max(max_freq,abs(omega_freq))
        #print a,b,freq,sign
        hamfunc='rotation_ham_' +str(hamNO)
        hamNO=hamNO+1
        output = StringIO.StringIO()
        output.write('def '+ hamfunc + '(t,args):\n')
        if chname == 'static_item':
#             output.write('\tprint t\n')
            output.write('\treturn exp(1j*('+str(omega_freq)+'*t*pi*2))\n')
            output.write('Ht.append([product*pi, '+hamfunc+'])\n')
        else:      
            output.write('\tpulses,channels = args[\'PULS\'],args[\'MWCH\']\n')
            output.write('\tamp,phase = pulses.get_shape(\''+chname+'\',t)\n')   
#             output.write('\tprint amp\n')      
#             output.write('\tprint t,amp\n')
            output.write('\treturn exp(1j*('+str(omega_freq)+'*t*pi*2' + \
                     '+channels[\''+chname+'\'].phase + phase))  * amp\n')
            output.write('Ht.append([product*pi*'+str(channels[chname].field*2)+', '+hamfunc+'])\n')

        exec output.getvalue()
        if debug :
            print output.getvalue()
            print_matrix(product,'product=')

        hamfunc='rotation_ham_' +str(hamNO)
        hamNO=hamNO+1
        output = StringIO.StringIO()
        output.write('def '+ hamfunc + '(t,args):\n')

        if chname=='static_item' :
            output.write('\treturn exp(-1j*('+str(omega_freq)+'*t*pi*2))\n')
        else:      
            output.write('\tpulses,channels = args[\'PULS\'],args[\'MWCH\']\n')
            output.write('\tamp,phase = pulses.get_shape(\''+chname+'\',t)\n')         
            output.write('\treturn exp(-1j*('+str(omega_freq)+'*t*pi*2' + \
                     '+channels[\''+chname+'\'].phase + phase))  * amp\n')

        output.write('Ht.append([product.dag()*pi'+('*'+str(channels[chname].field*2) if chname!='static_item' else '')+', '+hamfunc+'])\n')
        if debug :
            print output.getvalue()
            print_matrix(product,'product=')
            
        exec output.getvalue()
        output.close()

    #options for mesolve
    opts=Odeoptions()
    opts.nsteps=1e6
    args={'MWCH':channels}    
    return max_freq,args,opts,Ht

    
def cal_static_ham_rframe_items(sysdef,rframe,ham,s,threshold):
    H_static = 0
    H_nondiagnal=[]
    diagonal_ham_items = []
    nondiagonal_ham_items = []
    
    sperturb = Perturbation(ham,0,1.0)
    diagonal_item,nondiagonal_items = \
        sperturb.rotating_wave_approximate(rframe.P,rframe.omega,s)               
    H_static = H_static + diagonal_item
    
    for nondiagonal_item in nondiagonal_items:
        if(abs(nondiagonal_item[1])<threshold):
            H_static = H_static+nondiagonal_item[0]
        else:
            H_nondiagnal.append(nondiagonal_item)
            
    return H_static,H_nondiagnal 
 
def cal_mw_ham_rframe_items(sysdef,rframe,channel,s):
    H_diagnal = 0
    H_nondiagnal=[]
    mwperturb=MWFieldPerturbation(channel)            
    diagonal_item,nondiagonal_items = \
        mwperturb.rotating_wave_approximate(rframe.P,rframe.omega,s)

    if diagonal_item.full().any():
        H_diagnal=diagonal_item
    
    for nondiagonal_item in nondiagonal_items:
#         nondiagonal_item.append(chname)
        H_nondiagnal.append(nondiagonal_item)
    return H_diagnal,H_nondiagnal

from numbers import Number
def conver_matrix_to_matlab(m):
    if isinstance(m,Qobj):
        m = m.full()
    s = '['
    for i in range(len(m)):
        if isinstance(m[i],Number):
            s = s + str(m[i]) + ','
        else:
            for j in range(len(m[i])):
                s = s + str(m[i][j]) + ','
            s = s + ';'
    return s+']'



def propagator_thread_deprecated(pargs):
    args,Ht,trange,steps,show_process = pargs
    dt = (trange[1]-trange[0])/float(steps)
    U = 1
    report_loops = int(steps/10.0)
    for n in range(int(steps)):
        if show_process and n % report_loops == 0 : print 'calcuating ',n*100.0/float(steps),'% complete'
        t = trange[0] + n * dt
         
        ham = 0
        for i in range(0, len(Ht)):
            if isinstance(Ht[i],Qobj):
                ham = ham + Ht[i]
            else:
                ham = ham + Ht[i][0] * Ht[i][1](t, args) 
        V = (-1j*dt*ham).expm()
#         if(ham.full().any()):
#             print_matrix(V,'V=@'+str(n))
        U = V*U
    return U      
#     return tensor(qeye(3),qeye(2))


def cal_propagator_deprecated(args,Ht,trange,steps,num_cpus = 4):        
    if not isinstance(trange,list):
        trange = [0,trange]
#     dt=(trange[1]-trange[0])/float(steps)#thresh_err/freq if abs(freq)> 0 else (trange[1]-trange[0])
#     nslices = int(floor((trange[1]-trange[0])/dt))
#     print 'dt = ',dt 
    if num_cpus>1:
        if int(steps) % int(num_cpus) !=0:
            steps = ((int(steps) / int(num_cpus))+1)*int(num_cpus)
        
        arglists=[]
        substeps = steps / num_cpus
        t = 0
        subrange = (trange[1]-trange[0])/float(num_cpus)
        for i in range(num_cpus):
            if i == 0:
                arglists.append([args,Ht,[i*subrange+trange[0],(i+1)*subrange+trange[0]],substeps,True])
#                 arglists.append([])
            else:
#                 arglists.append([])
                arglists.append([args,Ht,[i*subrange+trange[0],(i+1)*subrange+trange[0]],substeps,False])
        Us = parfor(propagator_thread_deprecated,arglists,num_cpus)
        U=1
        for u in Us:
            U = u * U
        return U
    else:
        return propagator_thread_deprecated([args,Ht,trange,steps,True])

def propagator_thread(pargs):
    args,Ht,tlist,show_process = pargs
    U = 1
    steps = len(tlist)-1
    report_loops = int(steps/10.0)
    for n in range(int(steps)):
        if show_process and n % report_loops == 0 : print 'calcuating ',round(n*100.0/float(steps)),'% complete'
        t = tlist[n]
        dt = tlist[n+1]-t
        ham = 0
        for i in range(0, len(Ht)):
            if isinstance(Ht[i],Qobj):
                ham = ham + Ht[i]
            else:
                ham = ham + Ht[i][0] * Ht[i][1](t, args) 
#         V = (-1j*dt*ham).expm()
        U = (-1j*dt*ham).expm()*U
    return U      

def cal_propagator(args,Ht,tlist,num_cpus = 4):        
    steps = len(tlist)-1
    
    if num_cpus > 1:
        arglists=[]
        substeps = ceil(steps / num_cpus)
        range_begin = 0
        i = 0
        while i < num_cpus:
            range_end = range_begin + substeps
            if range_end >= len(tlist):
                range_end = len(tlist)-1
            if i == 0:
                arglists.append([args,Ht,tlist[range_begin:range_end+1],True])
            else:
                arglists.append([args,Ht,tlist[range_begin:range_end+1],False])
            i = i + 1
            range_begin = range_end
            
        Us = parfor(propagator_thread,arglists,num_cpus)
        U=1
        for u in Us:
            U = u * U
        return U
        
    else:
        return propagator_thread([args,Ht,tlist,True])
        

def Bdx(t,args):
    channels = args['MWCH']
    pslice = args['PSLICE']
    starting_time = args['ST_TIME']
#     print t
    b=0
    for ch in channels:
#         ampl,phase = pulses.get_shape(ch,t)
        b+=2*channels[ch].field*pslice.get_amplitude(ch)*cos(channels[ch].theta)*cos(channels[ch].phi)*cos(2*pi*channels[ch].frequency*(t+starting_time) + pslice.get_phase(ch)+channels[ch].phase)  
#     print pslice.get_amplitude('MW')
    return b

def Bdy(t,args):
    channels = args['MWCH']
    pslice = args['PSLICE']
    starting_time = args['ST_TIME']
    b=0
    for ch in channels:
#         ampl,phase = pulses.get_shape(ch,t)
        b+= 2*channels[ch].field*pslice.get_amplitude(ch)*cos(channels[ch].theta)*sin(channels[ch].phi)*cos(2*pi*channels[ch].frequency*(t+starting_time) +pslice.get_phase(ch)+ channels[ch].phase)
    return b

    
def Bdz(t,args):
    channels = args['MWCH']
    pslice = args['PSLICE']
    starting_time = args['ST_TIME']
    b=0
    for ch in channels:
#         ampl,phase = pulses.get_shape(ch,t)
        b+= 2*channels[ch].field*pslice.get_amplitude(ch)*sin(channels[ch].theta)*cos(2*pi*channels[ch].frequency*(t+starting_time) + pslice.get_phase(ch)+channels[ch].phase)
    return b

  


    
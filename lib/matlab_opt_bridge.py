from qutip import *
from nv_system import *
from numpy import *
from rotating_frame import *
from qcexp import *

import StringIO


def gen_ham_func_for_matlab(sysdef,rframe,Hleft,channels,s_coup,s_mw,P=None,threshold = 1e-6,debug=True):     
    static_nondiagnal_hams={}
    channel_ctl_hams=[]    
    H_static,H_nondiagnal = cal_static_ham_rframe_items(sysdef,rframe,Hleft,s_coup,threshold)   
    static_nondiagnal_func_name = 'ham_static_nondiagnal'
    sndout = StringIO.StringIO()
    sndout.write('function hamsn='+static_nondiagnal_func_name+'(t)')    
    sndout.write('\r\n')
    sndout.write('hamsn=0'+';')
    sndout.write('\r\n')
    
    for i in range(len(H_nondiagnal)):
        product,omega_freq= H_nondiagnal[i]
        if P :
            product = P*product*P
#         matrix_name = 'static_nondiag_ham'+str(i)  
        if product.full().any():
            sndout.write('product = '+ conver_matrix_to_matlab((product)*pi)+';')
            sndout.write('\r\n')
            sndout.write('hamsn=hamsn + exp(1j*('+str(omega_freq)+'*t*pi*2))*product+exp(-1j*('+str(omega_freq)+'*t*pi*2))*ctranspose(product)'+';')
            sndout.write('\r\n')
    static_nondiagnal_hams[static_nondiagnal_func_name]=sndout.getvalue()
    sndout.close()
    for chname in channels:
        channel = channels[chname]
        matfuncoutx = StringIO.StringIO()
        matfuncouty = StringIO.StringIO()
        xfuncname = str(chname)+'x'
        yfuncname = str(chname)+'y'
        matfuncoutx.write('function x='+xfuncname+'(t)')    
        matfuncouty.write('function y='+yfuncname+'(t)')    
        matfuncoutx.write('\r\n')    
        matfuncouty.write('\r\n')    
        matfuncoutx.write('x=0;')
        matfuncouty.write('y=0;')
        matfuncoutx.write('\r\n')    
        matfuncouty.write('\r\n')
        H_diagnal,H_nondiagnal = cal_mw_ham_rframe_items(sysdef,rframe,channel,s_mw)        
        if P :
            H_diagnal = P * H_diagnal * P
#         matrix_name = 'ctlop_diag_'+chname  
        if H_diagnal.full().any():
            matfuncoutx.write('product = '+ conver_matrix_to_matlab(H_diagnal*pi*2*channel.field*2)+';')
            matfuncoutx.write('\r\n')
            matfuncouty.write('product = '+ conver_matrix_to_matlab(H_diagnal*pi*2*channel.field*2)+';')
            matfuncouty.write('\r\n')
            matfuncoutx.write('x=x + cos('+str(channels[chname].frequency)+'*t*pi*2+'+str(channel.phase)+')*product'+';')
            matfuncouty.write('y=y - sin('+str(channels[chname].frequency)+'*t*pi*2+'+str(channel.phase)+')*product'+';')
            matfuncoutx.write('\r\n')    
            matfuncouty.write('\r\n')    

#         max_freq = max(max_freq,abs(channel.frequency))
    
        for i in range(len(H_nondiagnal)):
            product,omega_freq= H_nondiagnal[i]
            if P :
                product = P*product*P

#             max_freq = max(max_freq,abs(omega_freq))
#             matrix_name = 'ctlop_'+chname+'_'+str(i)  
  
            if product.full().any():
                matfuncoutx.write('product = '+ conver_matrix_to_matlab(product *pi*2*channel.field)+';')
                matfuncoutx.write('\r\n')
                matfuncouty.write('product = '+ conver_matrix_to_matlab(product *pi*2*channel.field)+';')
                matfuncouty.write('\r\n')
    
                matfuncoutx.write('x=x+exp(1j*('+str(omega_freq)+'*t*pi*2+'+str(channel.phase)+'))*product'+';')
                matfuncouty.write('y=y+1j*exp(1j*('+str(omega_freq)+'*t*pi*2+'+str(channel.phase)+'))*product'+';')      
                matfuncoutx.write('\r\n')    
                matfuncouty.write('\r\n')    
    
                matfuncoutx.write('x=x+exp(-1j*('+str(omega_freq)+'*t*pi*2+'+str(channel.phase)+'))*ctranspose(product)'+';')
                matfuncouty.write('y=y-1j*exp(-1j*('+str(omega_freq)+'*t*pi*2+'+str(channel.phase)+'))*ctranspose(product)'+';')
                matfuncoutx.write('\r\n')    
                matfuncouty.write('\r\n')    

        channel_ctl_hams.append([chname,matfuncoutx.getvalue(),matfuncouty.getvalue()]) 
        matfuncoutx.close()
        matfuncouty.close()

 
    return H_static*pi*2,static_nondiagnal_hams,channel_ctl_hams

def gen_ham_series_for_matlab(sysdef,rframe,Hleft,channels,s_coup,s_mw,tlist,P=None,threshold = 1e-6,debug=True):     
    static_nondiagnal_hams={}
    channel_ctl_hams=[]    
    H_static,H_nondiagnal = cal_static_ham_rframe_items(sysdef,rframe,Hleft,s_coup,threshold)   
    
    ham_static_nondiagnal = []
    for hi in range(len(tlist)):
        t = tlist[hi]
        ham = 0
        for i in range(len(H_nondiagnal)):
            product,omega_freq= H_nondiagnal[i]
            if P :
                product = P*product*P
            ham = ham + exp(1j*omega_freq*t*pi*2)*product*pi + exp(-1j*omega_freq*t*pi*2)*product.dag()*pi
        ham_static_nondiagnal.append(ham)
#         sndout.write('H_static_nondiagnal(:,:,'+str(hi+1)+')='+conver_matrix_to_matlab(ham)+';')
#         sndout.write('\r\n')
    chamsx = {}
    chamsy = {}
    for chname in channels:
        channel = channels[chname]
        chamsx[chname] = []
        chamsy[chname] = []
        H_diagnal,H_nondiagnal = cal_mw_ham_rframe_items(sysdef,rframe,channel,s_mw)        
             
        for t in tlist:
#             print t
            hamx =  cos(channel.frequency*t*pi*2 + channel.phase )*H_diagnal*pi*channel.field*4
            hamy =  -sin(channel.frequency*t*pi*2 + channel.phase )*H_diagnal*pi*channel.field*4
    
        
            for i in range(len(H_nondiagnal)):
                product,omega_freq= H_nondiagnal[i]
                product = product *pi*2*channel.field
    
                hamx = hamx + exp(1j*(omega_freq*t*pi*2+channel.phase))*product + exp(-1j*(omega_freq*t*pi*2 + channel.phase))*product.dag()
                hamy = hamy + 1j*exp(1j*(omega_freq*t*pi*2+channel.phase))*product -1j*exp(-1j*(omega_freq*t*pi*2+channel.phase))*product.dag()
            if P :
                hamx = P * hamx * P
                hamy = P * hamy * P
                
            chamsx[chname].append(hamx)
            chamsy[chname].append(hamy)
 
    return H_static*pi*2,ham_static_nondiagnal,chamsx,chamsy

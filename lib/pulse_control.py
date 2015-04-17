import scipy.io as sio
from qutip import *
from pylab import *
from numpy import *
from math import *
from datetime import *
from nv_system import *



axis_phase = {'x':0,'y':pi/2,'-x':pi,'-y':-pi/2}
        
class MWChannel:
    def __init__(self,freq,init_phase,field,theta=0,phi=0):
        self.frequency = freq
        self.phase = init_phase
        self.field = field
        self.phi = phi
        self.theta = theta        

        
class ArbWave:
    def __init__(self,wavestr,ratio=1.0,theta=0,phi=0):
        self.ratio = ratio
        self.waveform_str = wavestr
        self.phi = phi
        self.theta = theta 
                
    def genFuncStr(self,toffset):
        tempstr = str(self.waveform_str).replace('offset',str(toffset))
        return str(self.ratio)+'*('+tempstr+')'
    

    
class pulse_slice:

    def get_amplitude(self,channel):
        pass
            
    def get_phase(self,channel):
        pass
    
    def get_duration(self):
        pass
    
class simple_pulse(pulse_slice):
    
    def __init__(self,duration,pulses,axes={}):
        self.pulses = pulses
        self.duration = abs(duration)
        self.axes = axes


#     def get_shape(self,channel,t):
# #         print t
#         for pulse in self.pulses:
#             if channel in pulse.channel_names:
#                 return pulse.get_shape(channel,t)
#         return 0.0,0.0
    
    def get_amplitude(self,channel):
        if channel in self.pulses:
            return 1.0
        else:
            return 0.0
            
    def get_phase(self,channel):
        axis = self.axes[channel] if channel in self.axes else 'x'
        return axis_phase[axis]
       
    def get_duration(self):
        return self.duration 
         
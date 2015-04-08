import scipy.io as sio
from qutip import *
from pylab import *
from numpy import *
from math import *
from datetime import *
from nv_system import *


        
class MWChannel:
    def __init__(self,freq,init_phase,field,theta=0,phi=0):
        self.frequency = freq
        self.phase = init_phase
        self.field = field
        self.phi = phi
        self.theta = theta

#         Ix,Iy,Iz = eff_momentums
#         self.control_ham = -(Ix*cos(self.theta)*cos(self.phi) + \
#             Iy*cos(self.theta)*sin(self.phi)+\
#             Iz*sin(self.theta))          
    
class pulse_slice:

    def get_amplitude(self,channel):
        pass
            
    def get_phase(self,channel):
        pass
    
    def get_duration(self):
        pass
    
class simple_pulse(pulse_slice):
    
    def __init__(self,duration,pulses):
        self.pulses = pulses
        self.duration = abs(duration)


    def get_shape(self,channel,t):
#         print t
        for pulse in self.pulses:
            if channel in pulse.channel_names:
                return pulse.get_shape(channel,t)
        return 0.0,0.0
    
    def get_amplitude(self,channel):
        if channel in self.pulses:
            return 1.0
        else:
            return 0.0
            
    def get_phase(self,channel):
        return 0.0
       
    def get_duration(self):
        return self.duration 
    
# class single_channel_pulse(pulse_shape):   
#     def __init__(self,channame,duration,ampratio = 1.0,phase = 0.0):
#         self.duration = duration
#         self.ampl = ampratio
#         self.phase = phase
#         self.channel_names = [channame]
#         
#     def get_shape(self,channel,t):
#         if channel not in self.channel_names:
#             return 0.0,0.0
#         if t < self.duration and t >= 0:
#             return self.ampl,self.phase
#         else:
#             return 0.0,0.0
#             
#     def get_duration(self):
#         return self.duration

# class pulsecontrol:
#     def __init__(self,pulse_slices):
#         self.pulse_slices=pulse_slices
#         if not isinstance(pulses,list):
#             self.pulses = [pulses]
#         else:
#             self.pulses = pulses
#         self.start_time_index = []
#         start_time = 0
#         if len(self.pulses)>1:
#             for i in range(len(self.pulses)):
#                 self.start_time_index.append(start_time)
#                 start_time = start_time + self.pulses[i].get_duration()
#         self.duration = start_time
#         self.expected_operators = expected_operators
#         self.pointer = 0
#     
#     
#     def get_shape(self,channame,t):
#         
#         if self.duration != 0 :
#             if t > self.duration:
#                 return 0.0,0.0
#         else:
#             return self.pulses[0].get_shape(channame,t)
#         
#         while t < self.start_time_index[self.pointer]:
#             self.pointer = self.pointer + 1
#         while t >= self.start_time_index[self.pointer] +self.pulses[self.pointer].get_duration():
#             self.pointer = self.pointer - 1
#         
#         return self.pulses[self.pointer].get_shape(channame,t) 
      
    

# class shaped_pulse(pulse_shape):
#     def __init__(self,waveforms,slot_width):
#         #waveforms in format [[channelname1,x1,y1]...]
#         self.ampls={}
#         self.phas={}
#         self.durs={}
#         self.slot_width = slot_width
#         if not isinstance(waveforms[0],list):
#             waveforms = [waveforms]
#             
#         for waveform in waveforms:
#             chname,x,y=waveform
#             self.ampls[chname]=[]
#             self.phas[chname]=[]
#             self.durs[chname]=float(slot_width*len(x))
#             print 'duration of channel ',chname,' is ',self.durs[chname],' ns'
#             if not isinstance(x,list):
#                 x = [x]
#                 y = [y]
#             for i in range(len(x)):
#                 c = x[i]+1j*y[i]
#                 self.ampls[chname].append(norm(c))
#                 self.phas[chname].append(angle(c))
#     
#         self.duration = max(self.durs.values())
#         self.channel_names = self.ampls.keys()
        
#     def get_shape(self,channel,t):
# #         print t
#         if channel in self.ampls:
#             if (t < self.durs[channel] and t >= 0):
#                 n = int(floor(t/self.slot_width))
#                 return self.ampls[channel][n],self.phas[channel][n]
#             else:
#                 return 0.0,0.0
#             
#     def get_duration(self):
#         return self.duration

# class sliced_pulse(shaped_pulse):
#     def __init__(self,waveforms,slices,index=None):
#         #waveforms in format [[channelname1,x1,y1]...]
#         self.ampls={}
#         self.phas={}
#         self.durs={}
#         self.slices = slices
#         if index:
#             self.index = index
#         else:
#             self.index = [0.0]
#             for t in slices:
#                 self.index.append[t+self.index[-1]]
#             
#         if not isinstance(waveforms[0],list):
#             waveforms = [waveforms]
#             
#         for waveform in waveforms:
#             chname,x,y=waveform
#             self.ampls[chname]=[]
#             self.phas[chname]=[]
#             self.durs[chname]=sum(slices[:len(x)])
#             print 'duration of channel ',chname,' is ',self.durs[chname],' ns'
#             if not isinstance(x,list):
#                 x = [x]
#                 y = [y]
#             for i in range(len(x)):
#                 c = x[i]+1j*y[i]
#                 self.ampls[chname].append(norm(c))
#                 self.phas[chname].append(angle(c))
#     
#         self.duration = max(self.durs.values())
#         self.channel_names = self.ampls.keys()
#         self.pointer = 0
#         
#     def get_shape(self,channel,t):
# #         print t
#         if channel in self.ampls:
#             if (t < self.durs[channel] and t >= 0):
#                 while(t >= self.slices[self.pointer]+self.index[self.pointer]):
#                     self.pointer = self.pointer + 1
#                 while(t < self.index[self.pointer]):
#                     self.pointer = self.pointer - 1               
# #                 n = int(floor(t/self.slot_width))
#                 return self.ampls[channel][self.pointer],self.phas[channel][self.pointer]
#             else:
#                 return 0.0,0.0


# class pulse_sequence(pulse_shape):
#     def __init__(self,sequences):
#         self.sequences = []
#         self.channel_names = []
#         self.duration = 0
#         self.index = []
#         for pulse in sequences:
#             if isinstance(pulse,single_channel_pulse):
#                 self.sequences.append(pulse)
#                 self.channel_names.extend(pulse.channel_names)
#                 self.index.append(self.duration)
#                 self.duration = self.duration + pulse.get_duration()
#         self.pointer = 0
#             
#     def get_shape(self,channame,t):
#         if channame in self.channel_names :
#             if t >= 0 and t < self.duration:
#                 while(t >= self.sequences[self.pointer].get_duration()+self.index[self.pointer]):
#                     self.pointer = self.pointer + 1
#                 while(t < self.index[self.pointer]):
#                     self.pointer = self.pointer - 1               
#                 ret =  self.sequences[self.pointer].get_shape(channame,t-self.index[self.pointer])
#                 return ret
#             else:
#                 return 0.0,0.0
#     def get_duration(self):
#         return self.duration
        
                
# class gaussian_pulse:
#     def __init__(self,a,b,c):
#         self.a = a
#         self.b = b
#         self.c = c
#         
#     def cal_gaussian_function(self):
#         return self.a*exp(-(x-b)**2/c**2/2)
#     
#     def get_shape(self,channel,t):
#         if not channel in self.ampls:
#             return 0,0
#         if t > self.durs[channel]:
#             return 0,0
#         n = int(floor(t / self.slot_width))
#         return self.ampls[channel][n],self.phas[channel][n]
#                 

       
# class pulse_slice:
#     def __init__(self,*args):
#         self.shapes = {}
#         #args = chname1,pulse1,chname2,pulse2...target_operator
#         for i in range(len(args)/2):
#             self.shapes[args[i*2]]=args[i*2+1]
#         self.duration = 0
#         for pshape in self.shapes.values():
#             self.duration = max(self.duration,pshape.get_duration())
#             
#         if len(args) % 2 == 1:
#             self.target_operator = args[len(args)-1]
#         else:
#             self.target_operator = None
#     
#     def get_shape(self,ch,t):
#         return self.shapes[ch].get_shape(t)
#         
#     def get_duration(self):
#         return self.duration
        
# class pulse_controller:
#     def __init__(self,pulse=None):
#         self.times = [0]
#         self.pointer = 0
#         self.sequence = []
#         if pulse:
#             self.append(pulse)
#         
#     def append(self,pthing):
#         if isinstance(pthing,list) :
#             for pulse in pthing:
#                 self.sequence.append(pulse)               
#         elif isinstance(pthing,pulse_shape):
#             self.sequence.append(pthing)       
#             self.times.append(self.times[-1] + pthing.get_duration())
#             
#     def get_total_duration(self):
#         tduration = 0
#         for pulse in self.sequence:
#             tduration = tduration + pulse.duration
#         return tduration
#     
#     def get_shape(self,t):
#         if len(self.sequence) == 1 and self.sequence[0].get_duration()==0:
#             return self.sequence[0].get_shape(0)
#         
#         while t <= self.times[self.pointer]:
#             self.pointer = self.pointer - 1
#         while t > self.times[self.pointer] + self.sequence[self.pointer].get_duration() and self.pointer < len(self.times):
#             self.pointer = self.pointer + 1
#             
#         if self.pointer == len(self.times):
#             return 0.0,0.0
#         else:
#             return self.sequence[self.pointer].get_shape(t-self.times[self.pointer])
#             
            
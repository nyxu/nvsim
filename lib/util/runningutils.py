
from qutip import *
import matplotlib.pyplot as plt
from msgutils import *
file_suffix = '.nv'

class running_structure:
    """class for storing running information ,including argurments and results"""

    def __init__(self,intro,nametitle):
        self.argv = {}
        self.arg_intro ={}
        self.result = {}  
        self.result_intro = {}
        self.disp_data = [] 
        self.intro = intro
        self.title = nametitle
       

    def add_arg(self,key,value,intro=None):
        self.argv[key]=value
        if(intro != None):
            self.arg_intro[key] = intro
        else:
            self.arg_intro[key] = key
    
    def get_arg(self,key,default_val = None):
        val = self.argv[key]
        if val==None:
            return default_val
        else:
            return val
        
    def add_result(self,key,value,intro=None):
        self.result[key]=value
        if(intro != None):
            self.result_intro[key] = intro
        else:
            self.result_intro[key] = key
            
    def add_visual_data(self,data):
        self.disp_data.append(data)
    def get_visual_data(self):
        return self.disp_data
    
    splittor = '_'
    separator = ':'
    
    def get_file_name(self):
        name = self.title
        for key,val in self.argv.items():
            if isinstance(val,int) or isinstance(val,float):
                name =  name + self.splittor+ key + self.separator + str(self.argv[key])
                
        return str(name)
    
    def check_argv(self,arglist):
        for arg in arglist:
            if not (arg in self.argv):
                return False
        return True
    
    def __str__(self):
        info = self.intro + ' with parameters:\n'
        for key,val in self.argv.items():
            if isinstance(val,int) or isinstance(val,float):
                info = info +  key + '\t' + str(self.argv[key]) + '\n'
        return info


def nsave(running_struct,path):
    filename = os.path.join( path , running_struct.get_file_name()+file_suffix)
    if os.path.exists(filename):
        i = 1
        fname1 = os.path.join( path , running_struct.get_file_name() + '_' + str(i) +file_suffix) 
        while(os.path.exists(fname1)):
            i = i + 1
            fname1 = os.path.join( path , running_struct.get_file_name() + '_' + str(i) +file_suffix) 
        filename = fname1
    print 'save data to ',filename
    dirname = os.path.dirname(filename)
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    f = file(filename,'w')
    pickle.dump(running_struct,f)
    f.close()

def nload(path,iterative=False):
    if os.path.isfile(path):
        return __load_file(path)
    
    running_structs = []
    for file in os.listdir(path):
        if iterative and os.path.isdir(file):
            running_structs.extend(nload(file,iterative))
        elif os.path.isfile(file):
            running_structs.append(__load_file(file))
    return running_structs

def __load_file(filename):
    fname,fext = os.path.splitext(filename)
    if fext == file_suffix:
        f = file(filename,'r')
        obj = pickle.load(f)
        if isinstance(obj,running_structure):
           return obj   

class data2d:
    title = ''
    dlabel = 'x'
    vlabel = 'y'
    def __init__(self,title,domain,values,xlabel=None,ylabel=None,legends=None):
        self.domain = domain
        if len(domain) == len(values):
            self.values = [values]
        else:
            self.values = values
        self.title = title
        self.dlabel = xlabel
        self.vlabel = ylabel
        self.legends = legends
        
def show_data_2D(datum,savefile=None):
    if not isinstance(datum,list):
        datum = [datum]
    nfigs = len(datum) 
    if nfigs > 0:  
        nrows = int(sqrt(nfigs))
        ncols = nfigs/nrows + (1 if nfigs % nrows > 0 else 0)
        for i in range(nfigs):
            plt.subplot(nrows,ncols,i)
            for  j in range(len(datum[i].values)):
                val = datum[i].values[j]
                plt.plot(datum[i].domain, real(val),label = (datum[i].legends[j] if datum[i].legends != None else ('data ' + str(j))))
                
            plt.xlabel(datum[i].dlabel)
            plt.ylabel(datum[i].vlabel)
            plt.title(datum[i].title)
            plt.legend()
        if savefile == None:     
            plt.show()
        else:
            plt.savefig(savefile)
        
from scipy import fftpack
from numpy import *
def show_fft_2D(datum,remove_dc = True,saveImg=None):
    if not isinstance(datum,list):
        datum = [datum]
    nfigs = len(datum) 
    if nfigs > 0:  
        nrows = int(sqrt(nfigs))
        ncols = nfigs/nrows + (1 if nfigs % nrows > 0 else 0)
        for i in range(nfigs):
            data = datum[i]
            plt.subplot(nrows,ncols,i)
            deltat = data.domain[1]-data.domain[0]
            freqs = fftpack.fftfreq(len(data.domain),deltat)             
            for value in data.values:
                if remove_dc :
                    value = value - average(value)
                fft= fftpack.fft(value)     
                plt.plot(freqs, abs(fft))
                
            plt.xlabel('frequency of ' + data.dlabel)
            plt.ylabel('fft of ' + data.vlabel)
            plt.title('Spectrum of ' + data.title)
        if saveImg == None:
            plt.show()
        else:
            plt.savefig(saveImg)
        
def inform_result_2d(rstruct,email_addr):
    datum = rstruct.get_visual_data()
    show_fft_2D(datum,remove_dc = True,saveImg='temp.png')
    info = str(rstruct)
    print 'sending information to addr :',email_addr     
    send_result(email_addr,info,'temp.png')    

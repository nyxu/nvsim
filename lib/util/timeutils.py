from datetime import *

beginning = 'beginning'

class timeline:
    
    def __init__(self):
        self.begin()
    def begin(self):
        #self.begining = datetime.now()
        self.marks[beginning]=datetime.now()
        
    marks={}
    
    def mark(self,markname,time=None):
        if time != None:
            self.marks[markname]=time
        else:
            self.marks[markname]=datetime.now()
    
    def getmark(self,markname):
        return self.marks[markname]
    
    def getduration(self,*argv):
             
        if len(argv)==1:
            duration = self.marks[argv[0]]-self.marks[beginning]
        else :
            duration = self.marks[argv[1]]-self.marks[argv[0]]
        
        return '('+(str(duration.days)+' days ' if duration.days !=0 else '') + str(duration.seconds) + ' seconds)'
    
default_timeline = timeline()    
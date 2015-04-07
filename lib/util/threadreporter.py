
import Queue

import time
import os



        
        
def report(destfile,something):
    trytimes = 0
    if os.path.exists(destfile) :
        while (not os.access(destfile, os.W_OK)) and trytimes < 3:
            trytimes +=1 
            time.sleep(1)
    if trytimes >= 10:
        print 'cannot open file ',destfile,' for write'
        return
    
    f = open(destfile, 'a')
    

    if not isinstance(something,list):
        something = [something]
    for s in something:
        f.write(str(s)+'\t')
    f.write('\n')
    f.close()
        
            
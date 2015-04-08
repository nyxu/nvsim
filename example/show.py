import matplotlib.pyplot as plt
# from nv_system import *
# from pulse_control import *
# from rotating_frame import *
# from scipy import fftpack    
# from qcexp import *  
# from threadreporter import * 
# import pyplot as plot



if __name__=='__main__':
#     nvsys = nv_system([nv_electron_spin(),N14_nuclear_spin()],[[0,1,default_NV_N14_coupling]]) 
    filename = 'magnetometer_ac_process.txt'
    f = open(filename)
    content = f.read()
    f.close()
    
    lines = content.split('\n')
    res = {}
    for line in lines:
        words = line.split('\t')
#         print words[0]
        if len(words) and len(words[0]):
#             evo_list.append(float(words[0].strip()))
#             result_list.append(float(words[1].strip()))
            res[float(words[0].strip())]=double(words[1].strip())
            
    evo_list = res.keys()
    evo_list.sort()
    result_list = [] 
    for key in evo_list:
        result_list.append(res[key])

    # Plot the result
#     subplot(211)
    plt.plot(evo_list, result_list)
    plt.xlabel('Free Evolution Time')
    plt.ylabel('Occupation probability')
    plt.title('Hahn Echo in process')  
         
    plt.show()
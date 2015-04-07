from qutip import *
from pylab import *
from math import *
from nvsys import *
from scipy import *

def gap_no_resize(gap):
    return gap

def find_all_transitions(energies):
    result = [[] for k in range(len(energies))]
    for i in range(len(energies)):
        for j in range(i,len(energies)):
            deltam1 = []
            deltam2 = []
            for k in range(3):
                if abs(energies[i][k]-energies[j][k])== 1.0:
                    deltam1.append(k)
                else:
                    if abs(energies[i][k]-energies[j][k])== 2.0:
                        deltam2.append(k)
            if len(deltam1)==1 and len(deltam2)==0:
                result[i].append([j,deltam1[0],abs(energies[i][3]-energies[j][3])])
    return result

def convert_spectrum_to_tex(
                            energies,names=None,
                            maxlen = 40,
                            gap_resize_func=gap_no_resize,
                            level_st_x=2,
                            level_sp_x=9,
                            min_level_gap=2,
                            transgap = 3,
                            transoffset=0.55):
    gaps = []
    for i in range(1,len(energies)):
        gap = abs(energies[i][3]-energies[i-1][3])
        gaps.append(gap_resize_func(gap))       
    gaps = numpy.sort(gaps)
    #print gaps
    gapnorm = numpy.sum(gaps)
    gapratio = maxlen / gapnorm
    #regaps = [gaps[0]]
    
    level = 0 
    level_y_dims=[level]
    item = energies[0] 
    tex = '%draw the levels\n'
    
    tex = tex +  "\draw[level] ("+str(level_st_x)+"cm,"+ str(level)+\
        "em) -- ("+str(level_sp_x)+"cm,"+ str(level)+"em) node[right] {\ket{"+("0" if names==None else names[0])+"}=\ket{"+\
        str(item[0])+','+str(item[1])+','+('' if item[2]== 0.5 else '-')+\
        "\\frac{1}{2}}};\n" 
    for i in range(1,len(energies)):
        item = energies[i]
        gap = gap_resize_func(abs(energies[i][3]-energies[i-1][3]))
        if gap*gapratio > min_level_gap:
            level = level + gap*gapratio
        else:
            level = level + min_level_gap
            
        level_y_dims.append(level)
        line = "\draw[level] ("+str(level_st_x)+"cm,"+ str(level)+\
        "em) -- ("+str(level_sp_x)+"cm,"+ str(level)+"em) node[right] {\ket{"+(str(i) if names==None else names[i])+"}=\ket{"+\
        str(item[0])+','+str(item[1])+','+('' if item[2]== 0.5 else '-')+\
        "\\frac{1}{2}}};\n"
        tex = tex + line
        #cpos = cpos+4
        
    trans = find_all_transitions(energies)
    
    #print 'transitions is ',trans
    tex = tex + '%draw the transitions\n'
    colors = ['black', 'red', 'green', 'blue', 'cyan', 'magenta', 'yellow']
    for i in range(len(trans)):
        subtrans = trans[i]
        trans_x = level_st_x + i*transoffset
        format = '%.6g'
        for j in range(len(subtrans)):
            if subtrans[j][2] < 1 and subtrans[j][2] > 1e-4:
                transfreq = str(format % (subtrans[j][2]*1e3))+'MHz'
            else:
                if subtrans[j][2] < 1e-2 and subtrans[j][2] > 1e-7:
                    transfreq = str(format % (subtrans[j][2]*1e6))+'KHz'
                else:
                    transfreq = str(format % (subtrans[j][2]*1e9))+'Hz'
            line = '\draw[trans,color = '+colors[i%len(colors)]+'] ('+str(trans_x)+'cm,'+str(level_y_dims[i])+\
            'em) -- ('+str(trans_x)+'cm,'+str(level_y_dims[subtrans[j][0]])+'em) node[midway,'+\
            ('left' if (i % 2 == 0) else 'right')+'] {'\
            +transfreq+'};%('+str(energies[i][0])+','+str(energies[i][2])+')->('+str(energies[subtrans[j][0]][0])+','+\
            str(energies[subtrans[j][0]][2])+')\n'
            tex = tex + line
            trans_x = trans_x + transgap
    return tex



from texutils import *
from nvsys import *


def calEigenVals(nvsys,ms=0,estates=None,state_discrim_threshold=0.7):
    ham = nvsys.get_ham()
    if estates == None:
   	    evals,estates = ham.eigenstates()
    else:
        evals=[]
        for state in estates:
			evals.append((state.dag()*ham*state).tr())
    
    nv=ms
    C13=[0.5,-0.5]
    N14 = [1,0,-1]
    eres = []
    #print len(estates)
    for ns in N14:
        for cs in C13:
            eval,bas = nvsys.get_exact_eigenstate([ns,nv,cs],state_discrim_threshold)    
            eres.append([ns,nv,cs,eval])        
    return eres

def getHeaderAndTail(Bz,selms):
    
    header = """% Paste the following text into a tex file and use it in latex
% Nanyang Xu 2013

\\documentclass[10pt]{article}

% I only need the arrows for this one.
\\usepackage{tikz}
\\usetikzlibrary{arrows}

% Nice captions.
\\usepackage[hang,small,bf]{caption}
\\setlength{\\captionmargin}{25pt}

% New commands to keep things tidy.
\\newcommand{\\ket}[1]{$\\left|#1\\right\\rangle$}
\\newcommand{\\Om}[1]{\\small $\\omega_{#1}$}
\\newcommand{\\De}[1]{$\\Delta_{#1}$}
\\newcommand{\\Ga}[1]{$\\Gamma_{#1}$}
\\begin{document}

% Place the TikZ picture in a figure environment.
\\begin{figure}
\\centerline{
  % Resize it to 5cm wide.
  \\resizebox{10cm}{!}{
    \\begin{tikzpicture}[
      scale=1,
      level/.style={thick},
      virtual/.style={thick,densely dashed},
      trans/.style={thick,<->,shorten >=2pt,shorten <=2pt,>=stealth},
      used/.style={thick,double,<->,shorten >=4pt,shorten <=4pt,>=stealth}
    ]
"""
    tail = """\\end{tikzpicture}
  }
}
\\caption{Energy spectrum and transitions in ms="""+str(selms)+""" subspace under magnetic field """+str(Bz*1000)+"""G}
\\end{figure}
\\end{document} 
"""
    return header,tail

if __name__=='__main__':
    
    #gen nv lab frame options
    Bz=2 #kGauss
    nvsys = default_nv_triple_spin_system(B0=[0,0,Bz])

    #the selected ms
    selms=1
    
    energies = calEigenVals(nvsys,selms)
     
    tex = convert_spectrum_to_tex(energies)
    header,tail=getHeaderAndTail(Bz,selms)
    print 'spectrum for ms = ',str(selms),':'
    print header + tex + tail

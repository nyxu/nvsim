from qcexp import *
from nv_system import *
from pulse_control import *
from timeutils import *
from nvutils import *
from electron_single_operator import *

def test_spin_half(evo_time=0.5,evo_points=50):
    sx = sigmax()
    sy = sigmay()
    sz = sigmaz()
    II = qeye(2)
    
    omega = 40 #MHz
    Azz = 10 #MHz
    H0 = pi*2*omega*tensor(sx,II)/2
    print_matrix((H0/pi/2).full())

    H1 = pi*2*(omega*tensor(sx,II)/2 + Azz*tensor(sz,sz)/4)
    print_matrix((H1/pi/2).full())
    evolist = linspace(0,evo_time,evo_points)
    fides = []
    for t in evolist:
        U0 = (-1j*H0*t).expm()
        U1 = (-1j*H1*t).expm()
        fidelity = (U0.dag()*U1).tr()/4.0
        fides.append(fidelity)
        
    plot(evolist,fides)
    title('test hardpulse in spin 1/2')
    xlabel('evo_time(us)')
    ylabel('fidelity')
    show()

def test_spin1(evo_time=0.1,evo_points=50):
    sx = (basis(3,0)*basis(3,1).dag()+basis(3,1)*basis(3,0).dag())/sqrt(2)
#     sx = jmat(1,'x')
    sy = jmat(1,'y')
    sz = jmat(1,'z')
    Is = qeye(3)
    
    Ix = jmat(0.5,'x')
    Iy = jmat(0.5,'y')
    Iz = jmat(0.5,'z')
    II = qeye(2)
    omega = 40 #MHz
    Azz,Azx,Azy = 10,2.2,2.9 #MHz

    H0 = pi*2*omega*tensor(sx,II)
    print_matrix((H0/pi/2).full())

    H1 = pi*2*(omega*tensor(sx,II) + Azz*tensor(sz,Iz) + Azx*tensor(sz,Ix)+Azy*tensor(sz,Ix))
    print_matrix((H1/pi/2).full())
    evolist = linspace(0,evo_time,evo_points)
    fides = []
    for t in evolist:
        U0 = (-1j*H0*t).expm()
        U1 = (-1j*H1*t).expm()
        fidelity = (U0.dag()*U1).tr()/6.0
        fides.append(sqrt(real(fidelity)**2 + imag(fidelity)**2))
    
    plot(evolist,fides)
    title('test hardpulse in spin 1')
    xlabel('evo_time(us)')
    ylabel('fidelity')
    show()
    

if __name__=='__main__':   
#     test_spin_half()

    test_spin1()
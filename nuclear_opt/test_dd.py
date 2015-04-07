from qutip import *

from qcexp import *

if __name__=='__main__':
    es0 = basis(3,1)
    es1 = basis(3,2)
    ns0 = basis(3,1)
    ns1 = basis(3,0)
    coupling = 1e3
    tau = 5e-4
    ham = tensor(jmat(1,'z'),jmat(1,'z'))
    states = [tensor(es0,ns0),tensor(es0,ns1),tensor(es1,ns0),tensor(es1,ns1)]
    
    ex01 = tensor((es0*es1.dag() + es1*es0.dag())/2.0,qeye(3))
    ey01 = tensor(1j*(es0*es1.dag() - es1*es0.dag())/2.0,qeye(3))
    
    xpi = (-1j*ex01*pi).expm()
    ypi = (-1j*ey01*pi).expm()
    
    print_matrix(xpi,'xpi=')
    print_matrix(ypi,'ypi=')
    Ue = (-1j*coupling * ham*tau).expm()
    Ue2 = Ue*Ue
    
    print 'test hahn echo:'
    U = xpi*Ue*xpi*Ue
    print_matrix(U,'U = ')
        
    print 'test xy4:'
    U = Ue*xpi*Ue2*ypi*Ue2*ypi*Ue2*xpi*Ue
    print_matrix(U,'U = ')

    print 'test xy8:'
    U = Ue*xpi*Ue2*ypi*Ue2*xpi*Ue2*ypi*Ue2*ypi*Ue2*xpi*Ue2*ypi*Ue2*xpi*Ue
    print_matrix(U,'U = ')        
    
    
    
    
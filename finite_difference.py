#Finite difference
import numpy as np 
from scipy.linalg import solve_banded 

def bs_put_implicits(S0,K,T,r,sigma,M,N,Smax):
    """
    European put option pricing using implicit finite difference.

    M: number of spatial intervals
    N: number of time intervals

    Returns 
    -------------
    price: float 
        Numerical European put option price.
    """

    if Smax is None:
        Smax=4*max(S0,K)

    dS=Smax/M
    dt=T/N

    S=np.linspace(0,Smax,M+1)

    # Interior spatial indices
    j=np.arange(1,M,dtype=float)

    # Hull's implicit coefficients
    a=0.5*dt*(r*j-sigma**2*j**2)
    b=1+dt*(sigma**2*j**2+r)
    c=-0.5*dt*(r*j+sigma**2*j**2)

    # Construct the tridiagonal bounded matrix

    ab=np.zeros((3,M-1))
    ab[0,1:]=c[:-1] #upper diagonal
    ab[1,:]=b       #main diagonal
    ab[2,:-1]=a[1:] #lower diagonal

    # Terminal payoff
    V=np.maximum(K-S,0)

    # Backward time iteration
    for i in range(N-1,-1,-1):
        t=i*dt 
        #European put boundary conditions
        left_boundary=K*np.exp(-r*(T-t))
        right_boundary=0.0

        rhs=V[1:-1].copy()

        #Boundary contributions
        rhs[0]-=a[0]*left_boundary
        rhs[-1]-=c[-1]*right_boundary

        #Solve the tridiagonal linear system 
        V[1:-1]=solve_banded((1,1),ab,rhs)

        V[0]=left_boundary
        V[-1]=right_boundary

    return np.interp(S0,S,V)


def bs_put_CrankNicolson(S0,K,T,r,sigma,M,N,Smax):
    '''
    M: number of spatial intervals
    N: number of time intervals

    ---------------
    return : float
        European put price from Crank Nicolsom method
    '''

    if Smax is None:
        Smax=4*max(S0,K)

    dS=Smax/M
    dt=T/N

    S=np.linspace(0,Smax,M+1)
    
    # Interior price segments
    j=np.arange(1,M,dtype=float)

    # Coefficients
    a=dt*(r*j-sigma**2*j**2)/4
    b=1+dt*(sigma**2*j**2+r)/2
    c=-dt*(sigma**2*j**2+r*j)/4
    b_tilde=1-dt*(sigma**2*j**2+r)/2

    # Construct tridiagonal bounded matrix
    ab=np.zeros((3,M-1))
    ab[0,1:]=c[:-1]
    ab[1,:]=b 
    ab[2,:-1]=a[1:]

    V=np.maximum(K-S,0)

    for i in range(N-1,-1,-1):
        t=i*dt 

        # Boundary conditoins
        left_boundary=K*np.exp(-r*(T-t))
        right_boundary=0

        rhs=-a*V[:-2]+b_tilde*V[1:-1]-c*V[2:]

        # Boundary contributors
        rhs[0]-=a[0]*left_boundary
        rhs[-1]-=c[-1]*right_boundary

        # Solve the tridiagonal matrix equations
        V[1:-1]=solve_banded((1,1),ab,rhs)

        V[0]=left_boundary
        V[-1]=right_boundary

    return np.interp(S0,S,V)


def bs_put_explicit(S0,K,T,r,sigma,M,N,Smax):
    '''
    M: number of spatical intervals
    N: number of time intervals

    ------------
    return: float
        European put price by explicit finite difference PDE solver
    '''
    if Smax is None:
        Smax=4*max(S0,K)

    dS=Smax/M
    dt=T/N 

    if sigma**2*(M-1)**2*dt>1:
        print('Warning: time segment may be too large.')

    t=np.arange(N+1)*dt
    S=np.arange(M+1)*dS
    j=np.arange(M+1)

    # construct space-time price grid 
    c=np.zeros((N+1,M+1))
    c[N,:]=np.maximum(K-S,0)
    c[:,0]=K*np.exp(-r*(T-t))
    c[:,-1]=0

    # Coefficients
    a=(-r*j*dt+sigma**2*j**2*dt)/(2*(1+r*dt))
    b=(1-sigma**2*j**2*dt)/(1+r*dt)
    cc=(r*j*dt+sigma**2*j**2*dt)/(2*(1+r*dt))

    for i in range(N-1,-1,-1):
        c[i,1:-1]=a[1:-1]*c[i+1,:-2]+b[1:-1]*c[i+1,1:-1]+cc[1:-1]*c[i+1,2:]

    return np.interp(S0,S,c[0,:])










    

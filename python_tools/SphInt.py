import numpy as np
def fib(n,m ):
    """
    return the m largest values in a sequence of n fibonacci numbers in a sequence
    """
    seq_out=np.zeros(n)
    seq_out[0]=1
    seq_out[1]=1
    for idx in range(1,(n-1)):
        seq_out[idx+1]=seq_out[idx-1]+seq_out[idx]
    return seq_out[-m:]
def calc_z_Phi(n,phi_bound_2pi=True):
    """
    calulate z and Phi base on the perscription of
    https://doi.org/10.1088/0305-4470/37/48/005
    n is order of the fibronacci sequence
    returns
       z = cos(theta)
       phi either on range of 2pi or unbounded
    """
    Fs=fib(n,2)
    dz=2/Fs[1]
    z=np.arange(-1,1,dz)
    phi=(np.pi*Fs[0]/Fs[1])*np.arange(len(z))
    n2p=floor(phi/(np.pi))
    if phi_bound_2pi:
         phi_out=phi-np.pi*n2p
    else:
         phi_out=phi
    return (z,phi_out,dz)

def rtp2xyz(r,t,p):
    """
    input r theta and phi and output x y and z
    """
    z=r*np.cos(t)
    x=r*np.sin(t)*np.cos(p)
    y=r*np.sin(t)*np.sin(p)
    return x,y,z

def weight_sample(z):
    """ perform weighting of z to expedite error reduction
    """
    return(z + np.sin(np.pi*z)/np.pi)

def sph_integrate(func,n=17,nE=1,*fargs):
    """
    perform a spherical integration on the given
    func
    n is the sequence number of the fibronacci sequence
    Reference equation 8 in https://doi.org/10.1088/0305-4470/37/48/005
    """
    z,phi,dz = calc_z_Phi(n)
    zp = weight_sample(z)
    integrand_array=np.zeros((len(zp),nE))
    phi_pp=phi+np.pi # add pi to phi before calculating each function
    # use for loop to allow for non vectorizable functions
    for idx,zpj in enumerate(zp):
        integrand_array[idx,:]=func(phi[idx],zpj,*fargs)+func(phi_pp[idx],zpj,*fargs)
    zt=np.tile(z,(nE,1)).T
    integrand_array*=(1+np.cos(np.pi*zt))
    integral_out=integrand_array.sum(axis=0)*np.pi*dz
    return integral_out
    
def I_qtp_E(p,z,q,E,cell):
    """
    return an intensity based on q, z=cos(theta), and phi (p)
    """
    t=np.arccos(z)
    qv = np.zeros(3)
    qv[0],qv[1],qv[2] = rtp2xyz(q,t,p)
    bv = cell.getBasisVectors()
    hkl = qv.T*bv/2/np.pi

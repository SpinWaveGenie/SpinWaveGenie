import numpy as np

def ngauss(x,A,c,w):
    """
    produces n 1d gaussians
    x is the arange
    A,c, and w are parameters for the gaussians
    each one should be n in length
    A is the amplitude
    c is the center
    w is FWHM
    """
    sig=w/(2*np.sqrt(2*np.log(2)))
    Amat=np.tile(A,(len(x),1))
    cmat=np.tile(c,(len(x),1))
    sigmat=np.tile(sig,(len(x),1))
    xmat=np.tile(x,(len(A),1)).T
    #for idx in range(len(A)):
    tmp=Amat*np.exp(-((xmat-cmat)**2)/2/sigmat**2)
    #    out=out+tmp
    out=tmp.sum(axis=1)
    return out

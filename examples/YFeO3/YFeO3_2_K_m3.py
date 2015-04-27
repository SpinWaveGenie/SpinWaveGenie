import matplotlib
matplotlib.use('Agg')
import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt

rc('font',**{'family':'serif','serif':['Times'],'size':16.0})
rc('text', usetex=True)
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(4.0*2,3.0*2)
ax = plt.subplot(1,1,1)

#dispersion = np.transpose(np.loadtxt('FMChain.txt'))
image_array = np.loadtxt('YFeO3_2_K_m3.mat')
E_array = np.loadtxt('YFeO3_2_K_m3.y')
x_axis = np.loadtxt('YFeO3_2_K_m3.x')

print len(image_array)
print len(image_array[0])
print np.max(image_array)
im = ax.pcolorfast(x_axis[:,1],E_array,image_array,vmin=0.0,vmax=0.1)
#ax.scatter(dispersion[0],dispersion[3],c='white')
ax.set_ylabel(r'Energy (meV)')
ax.set_xlabel(r'distance $(\xi,0,0)$')
#ax.set_xlim(left=0.0,right=3.0)
#ax.set_ylim(bottom=0.0,top=5.0)
fig.colorbar(im);
plt.savefig('YFeO3_2_K_m3.png',dpi=400,bbox_inches='tight')
plt.close()

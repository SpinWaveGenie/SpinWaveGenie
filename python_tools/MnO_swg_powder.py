import ase.spacegroup as sg
import ase.neighborlist as nl
from ase import Atoms
import SpinWaveGenie as swg
import numpy as np
import matplotlib.pyplot as plt
import time
# The following three lines need to be updated once
# python_tools is fully integrated in the package.
import sys
sys.path.append('..')
import python_tools.SphInt as SInt

def calc_points(points):
    """
    """
    in_res=[]
    for idxh, point in enumerate(points):
        genie.createMatrix(*point)
        genie.calculate()
        in_res.append(genie.getPoints())
    return in_res

def get_neighbor_idx(m_at,latt,n_order,max_radii=3.):
    """
    get the indexes neighbors of order n_order
    0 is first order
    """
    n_at=latt.get_number_of_atoms()
    num_n=nl.NeighborList(np.ones(n_at)*max_radii,self_interaction=False,bothways=True)
    num_n.update(latt)
    idxs,vec=num_n.get_neighbors(m_at)
    dist_lst=latt.get_distances(m_at,idxs).round(4)
    distances=np.unique(dist_lst)
    nnidxs=idxs[np.where(dist_lst==distances[n_order])[0]]
    return nnidxs

# lattice parameters
a=4.46
# exchange parameters from ...
# spin wave genie defines the it as J/2
J1m=2*0.324
J1p=2*0.428
J2=2*0.420
spin=2.5
# use ase to define the lattice.  For spin wave genie one needs to define all
# atoms in the unit cell.
mn=sg.crystal(["Mn"],[(0,0,0)], spacegroup=225,
        cellpar=[a,a,a,90.,90.,90.],size=(2,2,2))
n_at=mn.get_number_of_atoms()
# define the magnetic sublattices to be FCC-II antiferromagnetism
# define 1st three atoms by hand and then loop through all NNN to be sure all
# spins are defined
# you can use ase.visualize.view to check if it is correct
mm=np.zeros(n_at)
mm[0:3]=spin
mm[3]=-spin
at_num=0
while np.prod(abs(mm))<1:
    nnnidxs=get_neighbor_idx(at_num,mn,1,max_radii=3.)
    mm[nnnidxs]=-mm[at_num]
    at_num+=1
mn.set_initial_magnetic_moments(mm)
inpos=mn.get_scaled_positions()
#Define spin direction as in the plane
# component perp to plane
Sx_vec=np.cross(inpos[1,:],inpos[2,:])/np.linalg.norm(np.cross(inpos[1,:],inpos[2,:]))
# component in plane and perpendicular to moment direction.
Sy_vec=inpos[1,:]/np.linalg.norm(inpos[1,:])
# Moment direction
Sz_vec=np.cross(Sx_vec,Sy_vec)
# start setting up spin wave genie
cell=swg.Cell()
# use unit cell from ase
incell=mn.get_cell_lengths_and_angles()
cell.setBasisVectors(*incell)
# generate sublattice with + spin direction
Mn1=swg.Sublattice()
Mn1.setName("Mn1")
Mn1.setType("MN2")
Mn1.setMoment(2.5,1.15026,5.8195)
#Mn1.setMoment(2.5,0,0)
#generate sublattice with - spin direction
Mn2=swg.Sublattice()
Mn2.setName("Mn2")
Mn2.setType("MN2")
Mn2.setMoment(2.5,1.99133,2.6779)
# add the subblatices to the cells
cell.addSublattice(Mn1)
cell.addSublattice(Mn2)
atom_pos=np.zeros(mn.positions.shape)
spin_dict={2.5:"Mn1",-2.5:"Mn2"}
for idx in range(len(mn.positions)):
    pos=inpos[idx]
    tlst=list(pos)
    tlst.insert(0,spin_dict[mn.arrays['magmoms'][idx]])
    print(tlst)
    cell.addAtom(*tlst)
    atom_pos[idx,:]=pos

factory=swg.InteractionFactory()
exchange1=factory.getExchange("J1m",J1m,"Mn1","Mn1",3.,3.2)
exchange1a=factory.getExchange("J1m",J1m,"Mn2","Mn2",3.,3.2)
exchange2=factory.getExchange("J1p",J1p,"Mn1","Mn2",3.,3.2)
exchange3=factory.getExchange("J2",J2,"Mn1","Mn2",4.3,4.6)
isotropy1=factory.getAnisotropy("D11",0.059,Sx_vec,"Mn1")
isotropy2=factory.getAnisotropy("D12",0.059,Sx_vec,"Mn2")
isotropy3=factory.getAnisotropy("D21",0.004,Sy_vec,"Mn1")
isotropy4=factory.getAnisotropy("D12",0.004,Sy_vec,"Mn2")
builder = swg.SpinWaveBuilder(cell)
builder.addInteraction(exchange1)
builder.addInteraction(exchange1a)
builder.addInteraction(exchange2)
builder.addInteraction(exchange3)
builder.addInteraction(isotropy1)
builder.addInteraction(isotropy2)
#builder# .addInteraction(1isotropy4)
genie = builder.createElement()
qpts=np.linspace(0,4,301)
#qpts=np.linspace(0,4,31)
energies=np.linspace(0.1,30.1,601)
#energies=swg.Energies(0.1,30.1,301)
I=np.zeros((len(qpts),len(energies)))
for idx,qpt in enumerate(qpts):
    I[idx,:]=SInt.sph_integrate(SInt.I_qtp_E,17,len(energies),qpt,energies,0.1,cell,genie)
    print ("%d %g \n"%(idx, time.clock()))
fig,ax=plt.subplots()
qptsm,em=np.meshgrid(qpts,energies)
ax.pcolormesh(qptsm,em,I.T,vmin=0,vmax=1)
fig.show()

import SpinWaveGenie as swg
import numpy as np

cell=swg.Cell()
cell.setBasisVectors(1.0,10.0,10.0,90.0,90.0,90.0);
S=1.0
J=-1.0
D=0.1
xhat=np.array([1.,0.,0.])
spin0=swg.Sublattice()
name0 = "Spin0"
spin0.setName(name0)
spin0.setType("NONE")
spin0.setMoment(S,0.0,0.0)
cell.addSublattice(spin0)
cell.addAtom(name0,0.0,0.0,0.0)

spin1=swg.Sublattice()
name1 = "Spin1"
spin1.setName(name1)
spin1.setType("NONE")
spin1.setMoment(S,np.pi,0.0)
cell.addSublattice(spin1)
cell.addAtom(name1,0.5,0.0,0.0)

factory = swg.InteractionFactory()
exchange = factory.getExchange("J",J,name0,name1,0.9,1.1)
anisotropy0 = factory.getAnisotropy("D",D,xhat,name0)
anisotropy1 = factory.getAnisotropy("D",D,xhat,name1)

builder=swg.SpinWaveBuilder(cell)
builder.addInteraction(exchange)
builder.addInteraction(anisotropy0)
builder.addInteraction(anisotropy1)

SW = builder.createElement()

Line=swg.PointsAlongLine()
Line.setFirstPoint(0.0,0.0,0.0)
Line.setFinalPoint(0.0,0.0,3.0*2.0*np.pi)
Line.setNumberPoints(201)
kPoints = Line.getPoints()
energies = Energies(0.2, 3.0, 201)

ODfactory = swg.OneDimensionalFactory()
gauss=ODfactory.getGaussian(0.15,1.0e-1)

SWP=swg.SpinWavePlot()
res=SWP.EnergyResolutionFunction(gauss,SW,energies)
cut=SWP.IntegrateThetaPhi(res,1e-1)

twodimcut=TwoDimensionalCut()
twodimcut.setFilename("AFMPowderAverage")
twodimcut.setPlotObject(cut)
twodimcut.setPoints(kPoints)
twodimcut.save()

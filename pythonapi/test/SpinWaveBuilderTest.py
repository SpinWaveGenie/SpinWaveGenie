import unittest
import SpinWaveGenie as swg
import numpy as np

class TestSpinWaveBuilderClass(unittest.TestCase):

    def test_(self):
        cell = swg.Cell()
        cell.setBasisVectors(1.0,10.0,10.0,90.0,90.0,90.0)
    
        spin = swg.Sublattice()
        name = "Spin0";
        spin.setName(name);
        spin.setType("NONE");
        spin.setMoment(1.0,0.0,0.0);
        cell.addSublattice(spin);
        cell.addAtom(name,0.0,0.0,0.0);

        factory = swg.InteractionFactory()
        exchange = factory.getExchange("J",1.0,name,name,0.9,1.1)
     
        builder = swg.SpinWaveBuilder(cell)
        builder.addInteraction(exchange)
        genie = builder.createElement()  

        genie.createMatrix(0.0, 0.0, 0.0);
        genie.calculate();
        pts = genie.getPoints();
        print len(pts)

        for pt in pts:
           print pt.frequency, pt.intensity 

if __name__ == '__main__':
    unittest.main()


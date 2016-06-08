import unittest
import SpinWaveGenie as swg
import numpy as np

class TestSpinWaveBuilderClass(unittest.TestCase):

    def test_(self):
        cell = swg.Cell()
        cell.setBasisVectors(1.0,10.0,10.0,90.0,90.0,90.0)

        J = 1.0
        S = 1.0

        spin = swg.Sublattice()
        name = "Spin0";
        spin.setName(name);
        spin.setType("NONE");
        spin.setMoment(S,0.0,0.0);
        cell.addSublattice(spin);
        cell.addAtom(name,0.0,0.0,0.0);

        factory = swg.InteractionFactory()
        exchange = factory.getExchange("J",J,name,name,0.9,1.1)

        builder = swg.SpinWaveBuilder(cell)
        builder.addInteraction(exchange)
        genie = builder.createElement()
        for k in np.linspace(0.0,3.0,num=61):
            genie.createMatrix(k, 0.0, 0.0)
            genie.calculate()
            pts = genie.getPoints()

            #analytical solution for frequency and intensity
            frequency = 2.0 * J * S * (1.0 - np.cos(2.0 * np.pi * k))
            intensity = S / 4.0
            for pt in pts:
                self.assertAlmostEqual(pt.frequency, frequency)
                if np.abs(pt.frequency) > 1.0e-5:
                    self.assertAlmostEqual(pt.intensity, intensity)


if __name__ == '__main__':
    unittest.main()


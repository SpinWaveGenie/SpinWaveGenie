import unittest
import PySpinWaveGenie as SWG
import numpy as np
try:
    import itertools.izip as zip
except ImportError:
    pass

class TestSublatticeClass(unittest.TestCase):

    def test_cubic_basis_vectors(self):
        ActualBasisVectors = 2.0*np.identity(3)
        ActualReciprocalVectors = np.pi*np.identity(3)

        test = SWG.Cell()
        test.setBasisVectors(2.0,2.0,2.0,90.0,90.0,90.0)
        CalculatedBasisVectors = test.getBasisVectors()

        for actual,calc in zip(ActualBasisVectors.flatten(),test.getBasisVectors().flatten()):
            self.assertAlmostEquals(calc,actual)
        for actual,calc in zip(ActualReciprocalVectors.flatten(),test.getReciprocalVectors().flatten()):
            self.assertAlmostEquals(calc,actual)

    def test_hexagonal_basis_vectors(self):
        ActualBasisVectors = np.array([[1.0,0.0,0.0],[-0.5,0.5*np.sqrt(3.0),0.0],[0.0,0.0,1.0]])
        ActualReciprocalVectors = np.array([[2.0*np.pi,2.0*np.pi/np.sqrt(3.0),0.0],[0.0,4.0*np.pi/np.sqrt(3.0),0.0],[0.0,0.0,2.0*np.pi]])

        test = SWG.Cell()
        test.setBasisVectors(1.0,1.0,1.0,90.0,90.0,120.0)

        for actual,calc in zip(ActualBasisVectors.flatten(),test.getBasisVectors().flatten()):
            self.assertAlmostEquals(calc,actual)
        for actual,calc in zip(ActualReciprocalVectors.flatten(),test.getReciprocalVectors().flatten()):
            self.assertAlmostEquals(calc,actual)

    def test_add_sublattice(self):
        test = SWG.Sublattice()
        test.setMoment(2.0,np.pi/2.0,np.pi)
        test.setName("SL1")
        test.setType("Fe3")

        SLtest = SWG.Cell()
        SLtest.addSublattice(test)
        
        self.assertEqual(len(SLtest),1)
        #self.assertRaises(ValueError,SLtest.addSublattice(test))

    def test_sublattice_position(self):
        test = SWG.Sublattice()
        test.setMoment(2.0,np.pi/2.0,np.pi)
        test.setName("SL1")
        test.setType("Fe3")

        test2 = SWG.Sublattice()
        test2.setMoment(2.0,np.pi/2.0,np.pi)
        test2.setName("SL2")
        test2.setType("Fe3")

        SLtest = SWG.Cell()
        SLtest.addSublattice(test)
        SLtest.addSublattice(test2)
        self.assertEqual(SLtest.getPosition("SL1"),0)
        self.assertEqual(SLtest.getPosition("SL2"),1)

    def test_add_atom(self):
        cell = SWG.Cell()
        cell.setBasisVectors(2.0,2.0,2.0,90.0,90.0,90.0)
        test = SWG.Sublattice()
        test.setMoment(2.0,np.pi/2.0,np.pi)
        test.setName("SL1")
        test.setType("Fe3")
        cell.addSublattice(test)
        cell.addAtom("SL1",0.0,0.0,0.0)
        cell.addAtom("SL1",0.5,0.5,0.5)

        self.assertEqual(len(cell.getSublattice("SL1")), 2)

        atom0 = np.array([0.0,0.0,0.0])
        FoundAtom0 = False
        atom1 = np.array([1.0,1.0,1.0])
        FoundAtom1 = False
        for point in cell.getSublattice("SL1"):
            test0 = point - atom0
            test1 = point - atom1
            if np.linalg.norm(test0) < 0.01:
                FoundAtom0 = True
            if np.linalg.norm(test1) < 0.01:
                FoundAtom1 = True
        self.assertTrue(FoundAtom0)
        self.assertTrue(FoundAtom1)

    def test_iterator(self):
        cell = SWG.Cell()
        test = SWG.Sublattice()
        test.setMoment(2.0,np.pi/2.0,np.pi)
        test.setName("SL1")
        test.setType("Fe3")

        test2 = SWG.Sublattice()
        test2.setMoment(2.0,np.pi/2.0,np.pi)
        test2.setName("SL2")
        test2.setType("Fe3")

        cell.addSublattice(test)
        cell.addSublattice(test2)

        foundSL1 = False
        foundSL2 = False
        foundBlank = False
        for elem in cell:
            if elem.getName() == "SL1":
                foundSL1 = True
            elif elem.getName() == "SL2":
                foundSL2 = True
            elif elem.getName() == "":
                foundBlank = True

        self.assertTrue(foundSL1)
        self.assertTrue(foundSL2)
        self.assertFalse(foundBlank)


if __name__ == '__main__':
    unittest.main()

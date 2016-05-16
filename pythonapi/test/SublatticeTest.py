import unittest
import PySpinWaveGenie as SWG
import numpy as np

class TestSublatticeClass(unittest.TestCase):

    def test_constructor(self):
        test = SWG.Sublattice()
        self.assertEqual(test.getName(),"")
        self.assertEqual(test.getType(),"None")
        self.assertAlmostEqual(test.getMoment(),0.0)
        self.assertAlmostEqual(test.getTheta(),0.0)
        self.assertAlmostEqual(test.getPhi(),0.0)

    def test_moment(self):
        test = SWG.Sublattice()
        test.setMoment(2.0,np.pi/2.0,np.pi)
        self.assertAlmostEqual(test.getMoment(),2.0)
        self.assertAlmostEqual(test.getTheta(),np.pi/2.0)
        self.assertAlmostEqual(test.getPhi(),np.pi)

    def test_rotation_matrices(self):
        test = SWG.Sublattice()
        RotationMatrix = test.getRotationMatrix()
        InverseMatrix = test.getInverseMatrix()
        self.assertTrue((RotationMatrix == np.identity(3)).all())
        self.assertTrue((InverseMatrix == np.identity(3)).all())

    def test_atoms(self):
        test = SWG.Sublattice()
        test.addAtom(0.0,0.0,0.0)
        test.addAtom(0.5,0.5,0.5)
        self.assertEqual(len(test), 2)

        atom0 = np.array([0.0,0.0,0.0])
        FoundAtom0 = False
        atom1 = np.array([0.25,0.25,0.25])
        FoundAtom1 = False
        for point in test:
            test0 = point - atom0
            test1 = point - atom1
            if np.linalg.norm(test0) < 0.01:
                FoundAtom0 = True
            if np.linalg.norm(test1) < 0.01:
                FoundAtom1 = True
        self.assertTrue(FoundAtom0)
        self.assertFalse(FoundAtom1)

if __name__ == '__main__':
    unittest.main()

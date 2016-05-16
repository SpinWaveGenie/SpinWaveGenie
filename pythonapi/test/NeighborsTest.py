import unittest
import PySpinWaveGenie as SWG
import numpy as np

def createMnBiCell():
    cell = SWG.Cell()
    cell.setBasisVectors(4.2827,4.2827,6.1103,90.0,90.0,120.0)

    Spin0 = SWG.Sublattice()
    Spin0.setName("Spin0")
    Spin0.setType("MN2")
    Spin0.setMoment(2.0,np.pi/2.0,0.0)
    cell.addSublattice(Spin0)
    cell.addAtom("Spin0",0.0,0.0,0.0)
    cell.addAtom("Spin0",0.0,0.0,0.5)

    return cell

class TestNeighborsClass(unittest.TestCase):
    def test_first_neighbors(self):
        cell = createMnBiCell()
        neighborList = SWG.Neighbors()
        neighborList.findNeighbors(cell,"Spin0", "Spin0", 3.0, 3.1)
        self.assertEqual(len(neighborList),2)
        for pt in neighborList:
            self.assertAlmostEqual(np.linalg.norm(pt),3.05515,places=6)

    def test_second_neighbors(self):
        cell = createMnBiCell()
        neighborList = SWG.Neighbors()
        neighborList.findNeighbors(cell,"Spin0", "Spin0", 4.2, 4.3)
        self.assertEqual(len(neighborList),6)
        for pt in neighborList:
            self.assertAlmostEqual(np.linalg.norm(pt),4.2827,places=6)

    def test_third_neighbors(self):
        cell = createMnBiCell()
        neighborList = SWG.Neighbors()
        neighborList.findNeighbors(cell,"Spin0", "Spin0", 5.2, 5.3)
        self.assertEqual(len(neighborList),12)
        for pt in neighborList:
            self.assertAlmostEqual(np.linalg.norm(pt),5.260747,places=6)

    def test_fourth_neighbors(self):
        cell = createMnBiCell()
        neighborList = SWG.Neighbors()
        neighborList.findNeighbors(cell,"Spin0", "Spin0", 6.0, 6.2)
        self.assertEqual(len(neighborList),2)
        for pt in neighborList:
            self.assertAlmostEqual(np.linalg.norm(pt),6.1103,places=6)

    def test_fifth_neighbors(self):
        cell = createMnBiCell()
        neighborList = SWG.Neighbors()
        neighborList.findNeighbors(cell,"Spin0", "Spin0", 7.40, 7.42)
        self.assertEqual(len(neighborList),6)
        for pt in neighborList:
            self.assertAlmostEqual(np.linalg.norm(pt),7.41785399,places=6)

    def test_sixth_neighbors(self):
        cell = createMnBiCell()
        neighborList = SWG.Neighbors()
        neighborList.findNeighbors(cell,"Spin0", "Spin0", 7.44, 7.49)
        self.assertEqual(len(neighborList),12)
        for pt in neighborList:
            self.assertAlmostEqual(np.linalg.norm(pt),7.46172134,places=6)

    def test_seventh_neighbors(self):
        cell = createMnBiCell()
        neighborList = SWG.Neighbors()
        neighborList.findNeighbors(cell,"Spin0", "Spin0", 7.9, 8.1)
        self.assertEqual(len(neighborList),12)
        for pt in neighborList:
            self.assertAlmostEqual(np.linalg.norm(pt),8.022375,places=6)
                
    def test_eight_neighbors(self):
        cell = createMnBiCell()
        neighborList = SWG.Neighbors()
        neighborList.findNeighbors(cell,"Spin0", "Spin0", 8.45, 8.65)
        self.assertEqual(len(neighborList),6)
        for pt in neighborList:
            self.assertAlmostEqual(np.linalg.norm(pt),8.5654,places=6)

    def test_nineth_neighbors(self):
        cell = createMnBiCell()
        neighborList = SWG.Neighbors()
        neighborList.findNeighbors(cell,"Spin0", "Spin0", 9.0, 9.1)
        self.assertEqual(len(neighborList),12)
        for pt in neighborList:
            self.assertAlmostEqual(np.linalg.norm(pt),9.093955,places=6)

    def test_tenth_neighbors(self):
        cell = createMnBiCell()
        neighborList = SWG.Neighbors()
        neighborList.findNeighbors(cell,"Spin0", "Spin0", 9.1, 9.2)
        self.assertEqual(len(neighborList),2)
        for pt in neighborList:
            self.assertAlmostEqual(np.linalg.norm(pt),9.16545,places=6)



if __name__ == '__main__':
    unittest.main()

import cmath, unittest
import SpinWaveGenie as swg

class TestThreeVectorsClass(unittest.TestCase):
    def test_insert(self):
        doubleTest = swg.ThreeVectors()
        complexTest = swg.ComplexThreeVectors()
        complexNumber = 1 + 1j
        doubleTest.insert(1.0, 1.0, 1.0)
        complexTest.insert(complexNumber, complexNumber, complexNumber)
        self.assertEqual(len(doubleTest), 1)
        self.assertEqual(len(complexTest), 1)

        it = iter(doubleTest)
        value = next(it)
        self.assertAlmostEqual(value[0], 1.0)
        self.assertAlmostEqual(value[1], 1.0)
        self.assertAlmostEqual(value[2], 1.0)

        it = iter(complexTest)
        value = next(it)
        self.assertAlmostEqual(value[0], complexNumber)
        self.assertAlmostEqual(value[1], complexNumber)
        self.assertAlmostEqual(value[2], complexNumber)
               
    def test_clear(self):
        doubleTest = swg.ThreeVectors()
        complexTest = swg.ComplexThreeVectors()
        complexNumber = 1 + 1j
        doubleTest.insert(1.0, 1.0, 1.0)
        complexTest.insert(complexNumber, complexNumber, complexNumber)
        self.assertEqual(len(doubleTest), 1)
        self.assertEqual(len(complexTest), 1)

        doubleTest.clear()
        complexTest.clear()

        self.assertEqual(len(doubleTest), 0)
        self.assertEqual(len(complexTest), 0)

    def test_iterator(self):
        doubleTest = swg.ThreeVectors()
        complexTest = swg.ComplexThreeVectors()
        complexNumber = 1 + 1j
        doubleTest.insert(1.0, 1.0, 1.0)
        doubleTest.insert(2.0, 2.0, 2.0)
        complexTest.insert(complexNumber, complexNumber, complexNumber)
        complexTest.insert(2.0*complexNumber, 2.0*complexNumber, 2.0*complexNumber)
        self.assertEqual(len(doubleTest), 2)
        self.assertEqual(len(complexTest), 2)

        it = iter(doubleTest)
        value = next(it)
        self.assertAlmostEqual(value[0], 1.0)
        self.assertAlmostEqual(value[1], 1.0)
        self.assertAlmostEqual(value[2], 1.0)

        value = next(it)
        self.assertAlmostEqual(value[0], 2.0)
        self.assertAlmostEqual(value[1], 2.0)
        self.assertAlmostEqual(value[2], 2.0)

        it = iter(complexTest)
        value = next(it)
        self.assertAlmostEqual(value[0], complexNumber)
        self.assertAlmostEqual(value[1], complexNumber)
        self.assertAlmostEqual(value[2], complexNumber)

        value = next(it)
        self.assertAlmostEqual(value[0], 2.0 * complexNumber)
        self.assertAlmostEqual(value[1], 2.0 * complexNumber)
        self.assertAlmostEqual(value[2], 2.0 * complexNumber)

if __name__ == '__main__':
    unittest.main()

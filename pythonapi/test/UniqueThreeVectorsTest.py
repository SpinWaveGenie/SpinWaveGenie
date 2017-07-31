import cmath, unittest
import SpinWaveGenie as swg

class TestThreeVectorsClass(unittest.TestCase):
    def test_insert(self):
        doubleTest = swg.UniqueThreeVectors()
        doubleTest.insert(1.0, 1.0, 1.0)
        self.assertEquals(len(doubleTest), 1)

        it = iter(doubleTest)
        value = next(it)
        self.assertAlmostEquals(value[0], 1.0)
        self.assertAlmostEquals(value[1], 1.0)
        self.assertAlmostEquals(value[2], 1.0)

    def test_clear(self):
        doubleTest = swg.ThreeVectors()
        doubleTest.insert(1.0, 1.0, 1.0)
        self.assertEquals(len(doubleTest), 1)

        doubleTest.clear()

        self.assertEquals(len(doubleTest), 0)

    def test_iterator(self):
        doubleTest = swg.ThreeVectors()
        doubleTest.insert(1.0, 1.0, 1.0)
        doubleTest.insert(2.0, 2.0, 2.0)
        self.assertEquals(len(doubleTest), 2)

        it = iter(doubleTest)
        value = next(it)
        self.assertAlmostEquals(value[0], 1.0)
        self.assertAlmostEquals(value[1], 1.0)
        self.assertAlmostEquals(value[2], 1.0)

        value = next(it)
        self.assertAlmostEquals(value[1], 2.0)
        self.assertAlmostEquals(value[2], 2.0)

if __name__ == '__main__':
    unittest.main()

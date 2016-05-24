import unittest
import SpinWaveGenie as swg
import numpy as np

class TestInteractionFactoryClass(unittest.TestCase):

    def test_getExchange(self):
        factory = swg.InteractionFactory()
        exchange = factory.getExchange("J",1.0,"SA","SB",2.0,2.1)
        self.assertEqual(exchange.getName(),"J")
        sl = exchange.sublattices()
        self.assertEqual(sl[0],"SA")
        self.assertEqual(sl[1],"SB")

    def test_getDzyaloshinskiiMoriya(self):
        factory = swg.InteractionFactory()
        dm = factory.getDzyaloshinskiiMoriya("D",1.0,[1.0,0.0,0.0],"SA","SB",2.0,2.1)
        self.assertEqual(exchange.getName(),"J")
        sl = dm.sublattices()
        self.assertEqual(sl[0],"SA")
        self.assertEqual(sl[1],"SB")

    def test_getAnisotropy(self):
        factory = swg.InteractionFactory()
        anisotropy = factory.getAnisotropy("K",1.0,[0.0,0.0,1.0],"SA")
        self.assertEqual(anisotropy.getName(),"K")
        sl = anisotropy.sublattices()
        self.assertEqual(sl[0],"SA")
        self.assertEqual(sl[1],"SA")

    def test_getMagneticField(self):
        factory = swg.InteractionFactory()
        magnetic_field = factory.getMagneticField("K",1.0,[0.0,0.0,1.0],"SA")
        self.assertEqual(magnetic_field.getName(),"K")
        sl = magnetic_fieldgit .sublattices()
        self.assertEqual(sl[0],"SA")
        self.assertEqual(sl[1],"SA")
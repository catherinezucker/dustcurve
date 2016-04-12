from unittest import TestCase

import dustcurve

class TestJoke(TestCase):
    def test_is_string(self):
        s = dustcurve.fillerfunc()
        self.assertTrue(isinstance(s, basestring))

from unittest import TestCase

import dustcurve

class TestFillerFunc(TestCase):
    def test_is_string(self):
        s = dustcurve.fillerfunc()
        self.assertTrue(isinstance(s, basestring))

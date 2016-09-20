from unittest import TestCase

import nanoraw

# not actually working, but provides a structure for testing
class TestCorrect(TestCase):
    def test_correct(self):
        corr_read = correct.correct_raw_data()

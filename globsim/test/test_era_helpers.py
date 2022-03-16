import unittest
from globsim.download.era_helpers import make_monthly_chunks
from datetime import datetime


class TestMakeMonthlyChunks(unittest.TestCase):

    def setUp(self) -> None:
        self.same_year_12_mo = make_monthly_chunks(datetime(2000, 1, 1), datetime(2000, 12, 31))
        self.partial_start = make_monthly_chunks(datetime(2000, 1, 29), datetime(2010, 12, 31))
        self.partial_end = make_monthly_chunks(datetime(2000, 1, 1), datetime(2004, 6, 15))

    def test_same_year(self):
        self.assertEqual(12, len(self.same_year_12_mo))

    def test_start_dates(self):
        self.assertEqual(3, len(self.partial_start[0]['day']))

    def test_end_dates(self):
        self.assertEqual(15, len(self.partial_end[-1]['day']))

    def test_correct_number_months(self):
        self.assertEqual(12, len(self.same_year_12_mo))
        self.assertEqual(10 * 12 + 12, len(self.partial_start))
        self.assertEqual(4 * 12 + 6, len(self.partial_end))

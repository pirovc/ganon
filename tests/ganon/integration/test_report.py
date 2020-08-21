import unittest, shlex, pickle, sys
from pathlib import Path
sys.path.append('src')
from ganon import ganon

class TestReport(unittest.TestCase):
    def test_report(self):
        """
        Test if report on sample data is working
        """
        path_data = "tests/ganon/integration/data/"
        output_file = "test_report.tre"

        ret = ganon.main("report", 
                        output_report=output_file, 
                        db_prefix=path_data+"sample_bacteria",
                        rep_file=path_data+"results.rep",
                        rank="all")

        # check if ran okay
        self.assertTrue(ret, "ganon report finish with an error")
       
        # check if files were created
        self.assertTrue(Path(output_file).is_file() , "File (" + output_file +") was not created")

if __name__ == '__main__':
    unittest.main()
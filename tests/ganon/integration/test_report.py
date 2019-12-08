import unittest
import shlex, pickle
from pathlib import Path
from src.ganon import ganon

class TestReport(unittest.TestCase):
    def test_report(self):
        """
        Test if report on sample data is working
        """
        path_data = "tests/ganon/integration/data/"
        output_file = "test_report.tre"
        ret = ganon.main(shlex.split("./ganon report --output-report "+output_file+" --db-prefix "+path_data+"sample_bacteria --rep-file "+path_data+"results.rep --rank all"))
        
        # check if ran okay
        self.assertFalse(ret, "ganon report finish with an error")
       
        # check if files were created
        self.assertTrue(Path(output_file).is_file() , "File (" + output_file +") was not created")

if __name__ == '__main__':
    unittest.main()
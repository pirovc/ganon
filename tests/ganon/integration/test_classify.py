import unittest
import shlex, pickle
from pathlib import Path
from src.ganon import ganon

class TestClassify(unittest.TestCase):
    def test_classify(self):
        """
        Test if classify on sample data is working
        """
        path_data = "tests/ganon/integration/data/"
        prefix = "test_classify"
        ret = ganon.main(shlex.split("ganon classify --db-prefix "+path_data+"sample_bacteria --reads "+path_data+"bacteria.simulated.1.fq -o "+prefix))
        
        # check if ran okay
        self.assertFalse(ret, "ganon classify finish with an error")
       
        # check if files were created
        for ext in ["lca", "out", "rep", "tre"]:
            self.assertTrue(Path(prefix+"."+ext).is_file() , "File (" + ext +") was not created") # TODO check file contents

if __name__ == '__main__':
    unittest.main()
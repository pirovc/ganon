import unittest, shlex, pickle, sys
from pathlib import Path
sys.path.append('src')
from ganon import ganon

class TestClassify(unittest.TestCase):
    def test_classify(self):
        """
        Test if classify on sample data is working
        """
        path_data = "tests/ganon/integration/data/"
        prefix = "test_classify"

        ret = ganon.main("classify", 
                        output_all=True, 
                        db_prefix=path_data+"sample_bacteria",
                        single_reads=path_data+"bacteria.simulated.1.fq",
                        output_prefix=prefix)

        # check if ran okay
        self.assertFalse(ret, "ganon classify finish with an error")
       
        # check if files were created
        for ext in ["lca", "all", "rep", "tre"]:
            self.assertTrue(Path(prefix+"."+ext).is_file() , "File (" + ext +") was not created") # TODO check file contents

if __name__ == '__main__':
    unittest.main()
import unittest
import sys
import os

sys.path.append('src')
from ganon import ganon
from ganon.config import Config

base_dir = "tests/ganon/"
sys.path.append(base_dir)
from utils import setup_dir
from utils import list_files_folder
from utils import list_sequences
from utils import build_sanity_check_and_parse
data_dir = base_dir + "data/"


class TestBuildCustomOffline(unittest.TestCase):

    results_dir = base_dir + "results/integration/build-custom/"

    default_params = {"input": data_dir + "build-custom/",
                      "ncbi_file_info": data_dir + "build-custom/assembly_summary.txt",
                      "taxonomy": "ncbi",
                      "taxonomy_files": data_dir + "build-custom/taxdump.tar.gz",
                      "threads": 1,
                      "write_info_file": True,
                      "keep_files": True,
                      "verbose": False,
                      "quiet": True}

    @classmethod
    def setUpClass(self):
        setup_dir(self.results_dir)

    def test_input_folder(self):
        """
        ganon build-custom with folder as --input with --extension
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_folder"
        params["input"] = data_dir + "build-custom/"
        params["input_extension"] = "fna.gz"
        cfg = Config("build-custom", **params)
        self.assertTrue(ganon.main(cfg=cfg), "ganon build-custom run failed")
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        files = list_files_folder(params["input"], params["input_extension"])
        self.assertTrue(res["target"]["file"].isin(files).all(), "Files missing from target")
        self.assertEqual(len(files), res["target"].shape[0], "Wrong number of files on target")
        self.assertTrue(res["info"]["file"].isin(files).all(), "Files missing from info")
        self.assertEqual(len(files), res["info"].shape[0], "Wrong number of files on info")

        # Wrong extension
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_folder_wrong_extension"
        params["input"] = data_dir + "build-custom/"
        params["input_extension"] = "xxx.gz"
        cfg = Config("build-custom", **params)
        self.assertFalse(ganon.main(cfg=cfg), "ganon build-custom ran but it should fail")

        # Wrong folder
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_folder_wrong_folder"
        params["input"] = data_dir + "wrong-place/"
        params["input_extension"] = "fna.gz"
        cfg = Config("build-custom", **params)
        self.assertFalse(ganon.main(cfg=cfg), "ganon build-custom ran but it should fail")

    def test_input_files(self):
        """
        ganon build-custom with files as --input
        """
        files = list_files_folder(data_dir + "build-custom/", "fna.gz")
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_files"
        params["input"] = files
        params["input_extension"] = ""
        cfg = Config("build-custom", **params)
        self.assertTrue(ganon.main(cfg=cfg), "ganon build-custom run failed")
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        self.assertTrue(res["target"]["file"].isin(files).all(), "Files missing from target")
        self.assertEqual(len(files), res["target"].shape[0], "Wrong number of files on target")
        self.assertTrue(res["info"]["file"].isin(files).all(), "Files missing from info")
        self.assertEqual(len(files), res["info"].shape[0], "Wrong number of files on info")

        # All files are invalid
        files = [f+".xxx" for f in files]
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_files_invalid"
        params["input"] = files
        params["input_extension"] = ""
        cfg = Config("build-custom", **params)
        self.assertFalse(ganon.main(cfg=cfg), "ganon build-custom ran but it should fail")

    def test_input_folders_files(self):
        """
        ganon build-custom with files and folders as --input with --extension
        """
        files = list_files_folder(data_dir + "build-custom/", "fna.gz")
        folder = data_dir + "build-custom/more/"
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_folders_files"
        params["input"] = files + [folder]
        params["input_extension"] = "fna.gz"
        cfg = Config("build-custom", **params)
        self.assertTrue(ganon.main(cfg=cfg), "ganon build-custom run failed")
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        files.extend(list_files_folder(folder, params["input_extension"]))
        self.assertTrue(res["target"]["file"].isin(files).all(), "Files missing from target")
        self.assertEqual(len(files), res["target"].shape[0], "Wrong number of files on target")
        self.assertTrue(res["info"]["file"].isin(files).all(), "Files missing from info")
        self.assertEqual(len(files), res["info"].shape[0], "Wrong number of files on info")

    def test_taxonomy(self):
        """
        ganon build-custom with --taxonomy ncbi,gtdb,skip
        """
        #ncbi
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_taxonomy_ncbi"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        cfg = Config("build-custom", **params)
        self.assertTrue(ganon.main(cfg=cfg), "ganon build-custom run failed")
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        #gtdb
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_taxonomy_gtdb"
        params["taxonomy"] = "gtdb"
        params["taxonomy_files"] = [data_dir + "build-custom/ar53_taxonomy.tsv.gz",
                                    data_dir + "build-custom/bac120_taxonomy.tsv.gz"]
        cfg = Config("build-custom", **params)
        self.assertTrue(ganon.main(cfg=cfg), "ganon build-custom run failed")
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        #skip
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_taxonomy_skip"
        params["taxonomy"] = "skip"
        params["taxonomy_files"] = ""
        cfg = Config("build-custom", **params)
        self.assertTrue(ganon.main(cfg=cfg), "ganon build-custom run failed")
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

    def test_input_target_file(self):
        """
        ganon build-custom with --input-target file
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_target_file"
        params["input_target"] = "file"
        cfg = Config("build-custom", **params)
        self.assertTrue(ganon.main(cfg=cfg), "ganon build-custom run failed")
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        files = list_files_folder(params["input"], "fna.gz")
        self.assertTrue(res["target"]["file"].isin(files).all(), "Files missing from target")
        self.assertEqual(len(files), res["target"].shape[0], "Wrong number of files on target")
        self.assertTrue(res["info"]["file"].isin(files).all(), "Files missing from info")
        self.assertEqual(len(files), res["info"].shape[0], "Wrong number of files on info")

    def test_input_target_sequence(self):
        """
        ganon build-custom with --input-target sequence
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_target_sequence"
        params["input_target"] = "sequence"
        params["taxonomy"] = "skip"
        cfg = Config("build-custom", **params)
        self.assertTrue(ganon.main(cfg=cfg), "ganon build-custom run failed")
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        sequences = list_sequences(list_files_folder(params["input"], "fna.gz"))
        self.assertTrue(res["target"]["sequence"].isin(sequences).all(), "Files missing from target")
        self.assertEqual(len(sequences), res["target"].shape[0], "Wrong number of files on target")
        self.assertTrue(res["info"]["target"].isin(sequences).all(), "Files missing from info")
        self.assertEqual(len(sequences), res["info"].shape[0], "Wrong number of files on info")

    def test_level_file(self):
        """
        ganon build-custom --input-target file and --level default, assembly or ranks
        """
        # --level default (file)
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_file_default"
        params["input_target"] = "file"
        cfg = Config("build-custom", **params)
        self.assertTrue(ganon.main(cfg=cfg), "ganon build-custom run failed")
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")
        # target (accession from file) should be final target
        self.assertTrue((res["target"]["target"] == res["info"]["target"]).all(), "Wrong target")
        # Should not have any specialization, target is the file
        self.assertTrue(res["info"]["specialization"].isna().all(), "Undefined specialization")
        self.assertTrue(res["info"]["specialization_name"].isna().all(), "Undefined specialization")
        # Tax must have "file" rank
        self.assertTrue("file" in res["tax"]._ranks.values(), "file rank not found")

        # --level genus (file)
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_file_genus"
        params["input_target"] = "file"
        params["level"] = "genus"
        cfg = Config("build-custom", **params)
        self.assertTrue(ganon.main(cfg=cfg), "ganon build-custom run failed")
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")
        # node should be final target
        self.assertTrue((res["target"]["target"] == res["info"]["node"]).all(), "Wrong target")
        # Should not have any specialization, target is the file
        self.assertTrue(res["info"]["specialization"].isna().all(), "Undefined specialization")
        self.assertTrue(res["info"]["specialization_name"].isna().all(), "Undefined specialization")
        # Tax must have "genus" rank and not "species"
        self.assertTrue("genus" in res["tax"]._ranks.values(), "file rank not found")
        self.assertFalse("species" in res["tax"]._ranks.values(), "file rank not found")

        # --level leaves (file)
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_file_leaves"
        params["input_target"] = "file"
        params["level"] = "leaves"
        cfg = Config("build-custom", **params)
        self.assertTrue(ganon.main(cfg=cfg), "ganon build-custom run failed")
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")
        # node should be final target
        self.assertTrue((res["target"]["target"] == res["info"]["node"]).all(), "Wrong target")
        # Should not have any specialization, target is the file
        self.assertTrue(res["info"]["specialization"].isna().all(), "Undefined specialization")
        self.assertTrue(res["info"]["specialization_name"].isna().all(), "Undefined specialization")
        # Tax must have "genus" and "species"
        self.assertTrue("genus" in res["tax"]._ranks.values(), "file rank not found")
        self.assertTrue("species" in res["tax"]._ranks.values(), "file rank not found")

        # --level assembly (file) - specialization
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_file_assembly"
        params["input_target"] = "file"
        params["level"] = "assembly"
        cfg = Config("build-custom", **params)
        self.assertTrue(ganon.main(cfg=cfg), "ganon build-custom run failed")
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")
        # specialization should be final target
        self.assertTrue((res["target"]["target"] == res["info"]["specialization"]).all(), "Wrong target")
        # Should have specialization
        self.assertFalse(res["info"]["specialization"].isna().any(), "Missing specialization")
        self.assertFalse(res["info"]["specialization_name"].isna().any(), "Missing specialization")
        # Tax must have "assembly" rank
        self.assertTrue("assembly" in res["tax"]._ranks.values(), "assembly rank not found")

    # def test_level_sequence(self):
    #     """
    #     ganon build-custom --input-target sequence and --level default, assembly or ranks
    #     """

    #     # --level default (sequence)
    #     params["db_prefix"] = self.results_dir + "test_level_sequence"
    #     params["input_target"] = "sequence"
    #     params["taxonomy"] = "skip"
    #     cfg = Config("build-custom", **params)
    #     self.assertTrue(ganon.main(cfg=cfg), "ganon build-custom run failed")
    #     res = build_sanity_check_and_parse(vars(cfg))
    #     self.assertIsNotNone(res, "ganon build-custom sanity check failed")
    #     # Should not have any specialization, target is the sequence
    #     self.assertTrue(res["info"]["specialization"].isna().all(), "Undefined specialization")

    #     print(res)

if __name__ == '__main__':
    unittest.main()

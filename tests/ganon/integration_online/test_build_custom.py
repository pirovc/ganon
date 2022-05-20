import unittest
import sys

sys.path.append('src')
from ganon import ganon
from ganon.config import Config

base_dir = "tests/ganon/"
sys.path.append(base_dir)
from utils import setup_dir
from utils import build_sanity_check_and_parse
data_dir = base_dir + "data/"


class TestBuildCustom(unittest.TestCase):

    results_dir = base_dir + "results/integration_online/build-custom/"

    default_params = {"input": data_dir + "build-custom/files/",
                      "taxonomy": "skip",
                      "threads": 1,
                      "write_info_file": True,
                      "keep_files": True,
                      "verbose": True,
                      "quiet": True}

    @classmethod
    def setUpClass(self):
        setup_dir(self.results_dir)

    def test_taxonomy(self):
        """
        ganon build-custom with --taxonomy ncbi,gtdb
        """
        #ncbi
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_taxonomy_ncbi"
        params["taxonomy"] = "ncbi"
        params["ncbi_file_info"] = data_dir + "build-custom/assembly_summary.txt"
        cfg = Config("build-custom", **params)
        self.assertTrue(ganon.main(cfg=cfg), "ganon build-custom run failed")
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        #gtdb
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_taxonomy_gtdb"
        params["taxonomy"] = "gtdb"
        cfg = Config("build-custom", **params)
        self.assertTrue(ganon.main(cfg=cfg), "ganon build-custom run failed")
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

    def test_level_sequence_default_gtdb(self):
        """
        ganon build-custom --input-target sequence and --level default
        needs eutils to link sequence -> accession -> gtdb tax
        """
        # --level default (sequence) - GTDB
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_sequence_default_gtdb"
        params["input_target"] = "sequence"
        params["taxonomy"] = "gtdb"
        params["taxonomy_files"] = [data_dir + "build-custom/ar53_taxonomy.tsv.gz",
                                    data_dir + "build-custom/bac120_taxonomy.tsv.gz"]
        cfg = Config("build-custom", **params)
        self.assertTrue(ganon.main(cfg=cfg), "ganon build-custom run failed")
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

    def test_level_sequence_taxrank_gtdb(self):
        """
        ganon build-custom --input-target sequence and --level {tax.rank}
        needs eutils to link sequence -> accession -> gtdb tax
        """
        # --level genus NCBI
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_sequence_taxrank_gtdb"
        params["input_target"] = "sequence"
        params["level"] = "genus"
        params["taxonomy"] = "gtdb"
        params["taxonomy_files"] = [data_dir + "build-custom/ar53_taxonomy.tsv.gz",
                                    data_dir + "build-custom/bac120_taxonomy.tsv.gz"]
        cfg = Config("build-custom", **params)
        self.assertTrue(ganon.main(cfg=cfg), "ganon build-custom run failed")
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")
        # Tax must not have "species" (filtered out)
        self.assertFalse("species" in res["tax"]._ranks.values(), "rank found")

    def test_level_sequence_leaves_gtdb(self):
        """
        ganon build-custom --input-target sequence and --level leaves
        needs eutils to link sequence -> accession -> gtdb tax
        """
        # --level leaves NCBI
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_sequence_leaves_ncbi"
        params["input_target"] = "sequence"
        params["level"] = "leaves"
        params["taxonomy"] = "gtdb"
        params["taxonomy_files"] = [data_dir + "build-custom/ar53_taxonomy.tsv.gz",
                                    data_dir + "build-custom/bac120_taxonomy.tsv.gz"]
        cfg = Config("build-custom", **params)
        self.assertTrue(ganon.main(cfg=cfg), "ganon build-custom run failed")
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

    def test_level_sequence_assembly(self):
        """
        ganon build-custom --input-target sequence and --level assembly (specialization)
        """
        # --level assembly no tax
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_file_assembly"
        params["input_target"] = "file"
        params["level"] = "assembly"
        params["ncbi_file_info"] = "refseq"
        cfg = Config("build-custom", **params)
        self.assertTrue(ganon.main(cfg=cfg), "ganon build-custom run failed")
        self.assertIsNotNone(build_sanity_check_and_parse(vars(cfg)), "ganon build-custom sanity check failed")

        # --level assembly NCBI
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_file_assembly_ncbi"
        params["input_target"] = "file"
        params["level"] = "assembly"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["ncbi_file_info"] = "refseq"
        cfg = Config("build-custom", **params)
        self.assertTrue(ganon.main(cfg=cfg), "ganon build-custom run failed")
        self.assertIsNotNone(build_sanity_check_and_parse(vars(cfg)), "ganon build-custom sanity check failed")

        # --level assembly GTDB
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_file_assembly_gtdb"
        params["input_target"] = "file"
        params["level"] = "assembly"
        params["taxonomy"] = "gtdb"
        cfg = Config("build-custom", **params)
        self.assertTrue(ganon.main(cfg=cfg), "ganon build-custom run failed")
        self.assertIsNotNone(build_sanity_check_and_parse(vars(cfg)), "ganon build-custom sanity check failed")

    def test_ncbi_sequence_info(self):
        """
        ganon build-custom --ncbi-sequence-info files
        """
        # Test with one, two or invalid file
        pass

    def test_ncbi_file_info(self):
        """
        ganon build-custom --ncbi-file-info files
        """
        # Test with one, two or invalid file
        pass


if __name__ == '__main__':
    unittest.main()

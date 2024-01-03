import unittest
import sys
import os

sys.path.append('src')
from ganon.config import Config
from ganon.util import download

base_dir = "tests/ganon/"
sys.path.append(base_dir)
from utils import run_ganon
from utils import setup_dir
from utils import build_sanity_check_and_parse
from utils import download_bulk_files
data_dir = base_dir + "data/"


class TestBuildCustom(unittest.TestCase):

    results_dir = base_dir + "results/integration_online/build-custom/"
    download_dir = base_dir + "downloads/"

    default_params = {"input": data_dir + "build-custom/files/",
                      "taxonomy": "skip",
                      "threads": 1,
                      "ncbi_url": "file://" + os.path.abspath(download_dir),
                      "gtdb_url": "file://" + os.path.abspath(download_dir),
                      "write_info_file": True,
                      "keep_files": True,
                      "verbose": True,
                      "quiet": True}

    @classmethod
    def setUpClass(self):
        setup_dir(self.results_dir)

        # Download full files once (and simulate download by setting ncbi_url and gtdb_url on params)
        download_bulk_files(self.download_dir)


    def test_taxonomy(self):
        """
        ganon build-custom with --taxonomy ncbi,gtdb (downloads taxonomy)
        """
        #ncbi
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_taxonomy_ncbi"
        params["taxonomy"] = "ncbi"
        params["ncbi_file_info"] = data_dir + "build-custom/assembly_summary.txt"
        cfg = Config("build-custom", **params)
        self.assertTrue(run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed")
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        #gtdb
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_taxonomy_gtdb"
        params["taxonomy"] = "gtdb"
        cfg = Config("build-custom", **params)
        self.assertTrue(run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed")
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
        params["filter_type"] = "ibf"
        params["input_target"] = "sequence"
        params["taxonomy"] = "gtdb"
        params["taxonomy_files"] = [data_dir + "build-custom/ar53_taxonomy.tsv.gz",
                                    data_dir + "build-custom/bac120_taxonomy.tsv.gz"]
        cfg = Config("build-custom", **params)
        self.assertTrue(run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed")
        res = build_sanity_check_and_parse(vars(cfg), skipped_targets=True)
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

    def test_level_sequence_taxrank_gtdb(self):
        """
        ganon build-custom --input-target sequence and --level {tax.rank}
        needs eutils to link sequence -> accession -> gtdb tax
        """
        # --level genus NCBI
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_sequence_taxrank_gtdb"
        params["filter_type"] = "ibf"
        params["input_target"] = "sequence"
        params["level"] = "genus"
        params["taxonomy"] = "gtdb"
        params["taxonomy_files"] = [data_dir + "build-custom/ar53_taxonomy.tsv.gz",
                                    data_dir + "build-custom/bac120_taxonomy.tsv.gz"]
        cfg = Config("build-custom", **params)
        self.assertTrue(run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed")
        res = build_sanity_check_and_parse(vars(cfg), skipped_targets=True)
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
        params["filter_type"] = "ibf"
        params["input_target"] = "sequence"
        params["level"] = "leaves"
        params["taxonomy"] = "gtdb"
        params["taxonomy_files"] = [data_dir + "build-custom/ar53_taxonomy.tsv.gz",
                                    data_dir + "build-custom/bac120_taxonomy.tsv.gz"]
        cfg = Config("build-custom", **params)
        self.assertTrue(run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed")
        res = build_sanity_check_and_parse(vars(cfg), skipped_targets=True)
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

    def test_level_sequence_assembly(self):
        """
        ganon build-custom --input-target sequence and --level assembly (specialization)
        """
        # --level assembly no tax
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_sequence_assembly"
        params["filter_type"] = "ibf"
        params["input_target"] = "sequence"
        params["level"] = "assembly"
        params["ncbi_sequence_info"] = ["nucl_gb"]
        cfg = Config("build-custom", **params)
        self.assertTrue(run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed")
        self.assertIsNotNone(build_sanity_check_and_parse(vars(cfg), skipped_targets=True), "ganon build-custom sanity check failed")

        # --level assembly NCBI
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_sequence_assembly_ncbi"
        params["filter_type"] = "ibf"
        params["input_target"] = "sequence"
        params["level"] = "assembly"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["ncbi_sequence_info"] = ["nucl_gb"]
        cfg = Config("build-custom", **params)
        self.assertTrue(run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed")
        self.assertIsNotNone(build_sanity_check_and_parse(vars(cfg), skipped_targets=True), "ganon build-custom sanity check failed")

        # --level assembly GTDB
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_sequence_assembly_gtdb"
        params["filter_type"] = "ibf"
        params["input_target"] = "sequence"
        params["level"] = "assembly"
        params["taxonomy"] = "gtdb"
        cfg = Config("build-custom", **params)
        self.assertTrue(run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed")
        self.assertIsNotNone(build_sanity_check_and_parse(vars(cfg), skipped_targets=True), "ganon build-custom sanity check failed")

    def test_ncbi_sequence_info(self):
        """
        ganon build-custom --ncbi-sequence-info files
        """
        # dead_nucl nucl_gb
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_ncbi_sequence_info"
        params["filter_type"] = "ibf"
        params["input_target"] = "sequence"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["ncbi_sequence_info"] = ["dead_nucl", "nucl_gb"]
        cfg = Config("build-custom", **params)
        self.assertTrue(run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed")
        res = build_sanity_check_and_parse(vars(cfg), skipped_targets=True)
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        # dead_nucl - none found
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_ncbi_sequence_info_wrong"
        params["filter_type"] = "ibf"
        params["input_target"] = "sequence"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["ncbi_sequence_info"] = ["dead_nucl"]
        cfg = Config("build-custom", **params)
        self.assertFalse(run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed")

        # eutils
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_ncbi_sequence_info_eutils"
        params["filter_type"] = "ibf"
        params["input_target"] = "sequence"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["ncbi_sequence_info"] = "eutils"
        cfg = Config("build-custom", **params)
        self.assertTrue(run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed")
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

    def test_ncbi_file_info(self):
        """
        ganon build-custom --ncbi-file-info files
        """
        # refseq
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_ncbi_sequence_info_refseq"
        params["input_target"] = "file"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["ncbi_file_info"] = "refseq"
        cfg = Config("build-custom", **params)
        self.assertTrue(run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed")
        res = build_sanity_check_and_parse(vars(cfg), skipped_targets=True)
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        # genbank
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_ncbi_sequence_info_refseq"
        params["input_target"] = "file"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["ncbi_file_info"] = "genbank"
        cfg = Config("build-custom", **params)
        self.assertTrue(run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed")
        res = build_sanity_check_and_parse(vars(cfg), skipped_targets=True)
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        # genbank refseq and historical
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_ncbi_sequence_info_genbank_refseq"
        params["input_target"] = "file"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["ncbi_file_info"] = ["refseq", "refseq_historical", "genbank", "genbank_historical"]
        cfg = Config("build-custom", **params)
        self.assertTrue(run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed")
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")


if __name__ == '__main__':
    unittest.main()

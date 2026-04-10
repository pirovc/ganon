import unittest
import os
from ganon.config import Config
from tests.ganon.utils import (
    run_ganon,
    setup_dir,
    list_files_folder,
    list_sequences,
    build_sanity_check_and_parse,
    write_input_file,
)
from parameterized import parameterized_class

base_dir = "tests/ganon/"
data_dir = base_dir + "data/"


@parameterized_class(
    [
        {"filter_type": "ibf"},
        {"filter_type": "hibf"},
    ]
)
class TestBuildCustom(unittest.TestCase):
    default_params = {
        "input": data_dir + "build-custom/files/",
        "taxonomy": "skip",
        "threads": 1,
        "write_info_file": True,
        "keep_files": True,
        "verbose": True,
        "quiet": False,
    }

    @classmethod
    def setUpClass(self):
        # parametrization for filter type
        self.default_params["filter_type"] = self.filter_type
        self.results_dir = (
            base_dir + "results/integration/build-custom_" + self.filter_type + "/"
        )
        setup_dir(self.results_dir)

    def test_input_folder(self):
        """
        ganon build-custom with folder as --input with --input-extension
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_folder"
        params["input"] = data_dir + "build-custom/files/"
        params["input_extension"] = "fna.gz"
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        files = list_files_folder(params["input"], params["input_extension"])
        self.assertTrue(
            res["target"]["file"].isin(files).all(), "Files missing from target"
        )
        self.assertEqual(
            len(files), res["target"].shape[0], "Wrong number of files on target"
        )
        self.assertTrue(
            res["info"]["file"].isin(files).all(), "Files missing from info"
        )
        self.assertEqual(
            len(files), res["info"].shape[0], "Wrong number of files on info"
        )

        # Wrong extension
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_folder_wrong_extension"
        params["input"] = data_dir + "build-custom/files/"
        params["input_extension"] = "xxx.gz"
        cfg = Config("build-custom", **params)
        self.assertFalse(
            run_ganon(cfg, params["db_prefix"]),
            "ganon build-custom ran but it should fail",
        )

        # Wrong folder
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_folder_wrong_folder"
        params["input"] = data_dir + "wrong-place/"
        params["input_extension"] = "fna.gz"
        cfg = Config("build-custom", **params)
        self.assertFalse(
            run_ganon(cfg, params["db_prefix"]),
            "ganon build-custom ran but it should fail",
        )

    def test_input_folder_recursive(self):
        """
        ganon build-custom with folder as --input with --input-extension and --input-recursive
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_folder_recursive"
        params["input"] = data_dir + "build-custom/files/"
        params["input_extension"] = "fna.gz"
        params["input_recursive"] = True

        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        # list files from base folder and "more" (got recursively)
        files = list_files_folder(
            params["input"], ext=params["input_extension"], recursive=True
        )

        self.assertTrue(
            res["target"]["file"].isin(files).all(), "Files missing from target"
        )
        self.assertEqual(
            len(files), res["target"].shape[0], "Wrong number of files on target"
        )
        self.assertTrue(
            res["info"]["file"].isin(files).all(), "Files missing from info"
        )
        self.assertEqual(
            len(files), res["info"].shape[0], "Wrong number of files on info"
        )

    def test_input_single_file(self):
        """
        ganon build-custom with one file as --input  and --input-target sequence
        """
        files = list_files_folder(data_dir + "build-custom/files/", ext="fna.gz")
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_single_file"
        params["input"] = files[0]
        params["input_extension"] = ""
        params["input_target"] = "sequence"
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        sequences = list_sequences([params["input"]])
        self.assertTrue(
            res["target"]["target"].isin(sequences).all(),
            "Sequence missing from target",
        )
        self.assertEqual(
            len(sequences),
            res["target"].shape[0],
            "Wrong number of sequences on target",
        )
        self.assertTrue(
            res["info"]["target"].isin(sequences).all(), "Sequences missing from info"
        )
        self.assertEqual(
            len(sequences), res["info"].shape[0], "Wrong number of sequences on info"
        )

    def test_input_files(self):
        """
        ganon build-custom with files as --input
        """
        files = list_files_folder(data_dir + "build-custom/files/", ext="fna.gz")
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_files"
        params["input"] = files
        params["input_extension"] = ""
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        self.assertTrue(
            res["target"]["file"].isin(files).all(), "Files missing from target"
        )
        self.assertEqual(
            len(files), res["target"].shape[0], "Wrong number of files on target"
        )
        self.assertTrue(
            res["info"]["file"].isin(files).all(), "Files missing from info"
        )
        self.assertEqual(
            len(files), res["info"].shape[0], "Wrong number of files on info"
        )

        # All files are invalid
        files = [f + ".xxx" for f in files]
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_files_invalid"
        params["input"] = files
        params["input_extension"] = ""
        cfg = Config("build-custom", **params)
        self.assertFalse(
            run_ganon(cfg, params["db_prefix"]),
            "ganon build-custom ran but it should fail",
        )

    def test_input_folders_files(self):
        """
        ganon build-custom with files and folders as --input with --input-extension
        """
        files = list_files_folder(data_dir + "build-custom/files/", ext="fna.gz")
        folder = data_dir + "build-custom/files/more/"
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_folders_files"
        params["input"] = files + [folder]
        params["input_extension"] = "fna.gz"
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        files.extend(list_files_folder(folder, ext=params["input_extension"]))
        self.assertTrue(
            res["target"]["file"].isin(files).all(), "Files missing from target"
        )
        self.assertEqual(
            len(files), res["target"].shape[0], "Wrong number of files on target"
        )
        self.assertTrue(
            res["info"]["file"].isin(files).all(), "Files missing from info"
        )
        self.assertEqual(
            len(files), res["info"].shape[0], "Wrong number of files on info"
        )

    def test_taxonomy_download(self):
        """
        ganon build-custom with --taxonomy ncbi,gtdb,skip simulating download with local files
        """
        # ncbi
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_taxonomy_download_ncbi"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"

        # Simulate download from local files (assembly_summary and species_genome_size)
        params["ncbi_url"] = (
            "file://" + os.path.abspath(data_dir) + "/build-custom/remote/"
        )
        params["ncbi_file_info"] = ["refseq", "genbank"]

        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        # gtdb
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_taxonomy_download_gtdb"
        params["taxonomy"] = "gtdb"
        params["taxonomy_files"] = [
            data_dir + "build-custom/ar53_taxonomy.tsv.gz",
            data_dir + "build-custom/bac120_taxonomy.tsv.gz",
        ]

        # Simulate download from local files (ar and bac metadata)
        params["gtdb_url"] = "file://" + os.path.abspath(data_dir) + "/build-custom/"

        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

    def test_taxonomy(self):
        """
        ganon build-custom with --taxonomy ncbi,gtdb,skip
        """
        # ncbi
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_taxonomy_ncbi"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["ncbi_file_info"] = data_dir + "build-custom/assembly_summary.txt"
        params["genome_size_files"] = (
            data_dir + "build-custom/species_genome_size.txt.gz"
        )
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        # gtdb
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_taxonomy_gtdb"
        params["taxonomy"] = "gtdb"
        params["taxonomy_files"] = [
            data_dir + "build-custom/ar53_taxonomy.tsv.gz",
            data_dir + "build-custom/bac120_taxonomy.tsv.gz",
        ]
        params["genome_size_files"] = [
            data_dir + "build-custom/ar53_metadata.tsv.gz",
            data_dir + "build-custom/bac120_metadata.tsv.gz",
        ]
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        # skip
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_taxonomy_skip"
        params["taxonomy"] = "skip"
        params["taxonomy_files"] = ""
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

    def test_taxonomy_keep_invalid_taxa(self):
        """
        ganon build-custom with --taxonomy ncbi --keep-invalid-taxa
        """

        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_taxonomy_keep_invalid_taxa"
        params["taxonomy"] = "ncbi"
        params["level"] = "strain"
        params["keep_invalid_taxa"] = True
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["ncbi_file_info"] = data_dir + "build-custom/assembly_summary.txt"
        params["genome_size_files"] = (
            data_dir + "build-custom/species_genome_size.txt.gz"
        )
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        # Only one entry at strain level, others set to root
        self.assertCountEqual(res["info"]["node"].unique(), ["1", "871271"])

    def test_convert_taxonomy_gtdb_gtdb(self):
        """
        ganon build-custom with --taxonomy gtdb-95 --convert-taxonomy gtdb-226

        same RS_GCF_900200805.1
            95: s__Neisseria meningitidis -> 226: s__Neisseria meningitidis
        missing GB_GCA_003520315.1
            95: s__Bact-08 sp003520315 -> 226:
        one-to-one RS_GCF_003473685.1
            95: s__Ruminococcus_A sp003011855 -> 226: s__Oliverpabstia intestinalis
        one-to-many RS_GCF_002198735.1
            95: g__JOSHI-001 -> 226: g__Aquabacterium_A, g__AHLZ01 (lca: f__Burkholderiaceae)

        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_convert_taxonomy_gtdb_gtdb"
        params["input"] = None
        params["input_file"] = data_dir + "build-custom/convert/convert_gtdb.tsv"
        params["input_target"] = "sequence"
        params["skip_genome_size"] = True
        params["level"] = "leaves"
        params["taxonomy"] = "gtdb-95"
        params["taxonomy_files"] = [
            data_dir + "build-custom/convert/bac120_taxonomy_r95.tsv.gz",
        ]
        params["convert_taxonomy"] = "gtdb-226"
        params["convert_taxonomy_files"] = [
            data_dir + "build-custom/convert/bac120_taxonomy_r226.tsv.gz",
        ]
        params["convert_gtdb_files"] = [
            data_dir + "build-custom/convert/95_acc_rep_lin_ncbi.tsv.gz",
            data_dir + "build-custom/convert/226_acc_rep_lin_ncbi.tsv.gz",
        ]
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg), skipped_targets=True)
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        # Lost one entry
        # Check conversion
        self.assertCountEqual(
            res["target"]["target"],
            [
                "s__Neisseria meningitidis",
                "s__Oliverpabstia intestinalis",
                "f__Burkholderiaceae",
            ],
        )

    def test_convert_taxonomy_gtdb_ncbi(self):
        """
        ganon build-custom with
            --taxonomy gtdb-95
            --convert-taxonomy ncbi-latest
            --level family

        RS_GCF_900200805.1
            95: s__Neisseria meningitidis -> ncbi: 481
        GB_GCA_003520315.1
            95: s__Bact-08 sp003520315 -> ncbi: 171550
        RS_GCF_003473685.1
            95: s__Ruminococcus_A sp003011855 -> ncbi: 186803
        RS_GCF_002198735.1
            95: g__JOSHI-001 -> ncbi: 2975441

        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_convert_taxonomy_gtdb_ncbi"
        params["input"] = None
        params["input_file"] = data_dir + "build-custom/convert/convert_gtdb.tsv"
        params["input_target"] = "sequence"
        params["skip_genome_size"] = True
        params["level"] = "family"
        params["taxonomy"] = "gtdb-95"
        params["taxonomy_files"] = [
            data_dir + "build-custom/convert/bac120_taxonomy_r95.tsv.gz",
        ]
        params["convert_taxonomy"] = "ncbi-latest"
        params["convert_taxonomy_files"] = [
            data_dir + "build-custom/convert/convert_nodes.dmp",
        ]
        params["convert_gtdb_files"] = [
            data_dir + "build-custom/convert/95_acc_rep_lin_ncbi.tsv.gz",
        ]
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        # Kept all entries
        # Check conversion to family nodes
        self.assertCountEqual(
            res["target"]["target"], ["481", "171550", "186803", "2975441"]
        )

    def test_convert_taxonomy_ncbi_gtdb(self):
        """
        ganon build-custom with
            --taxonomy ncbi
            --convert-taxonomy gtdb-226
            --level species

        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_convert_taxonomy_ncbi_gtdb"
        params["input"] = None
        params["input_file"] = data_dir + "build-custom/convert/convert_ncbi.tsv"
        params["input_target"] = "sequence"
        params["skip_genome_size"] = True
        params["level"] = "species"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = [
            data_dir + "build-custom/convert/convert_nodes.dmp",
        ]
        params["convert_taxonomy"] = "gtdb-226"
        params["convert_taxonomy_files"] = [
            data_dir + "build-custom/convert/bac120_taxonomy_r226.tsv.gz",
        ]
        params["convert_gtdb_files"] = [
            data_dir + "build-custom/convert/226_acc_rep_lin_ncbi.tsv.gz",
        ]
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg), skipped_targets=True)
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        # one entry cannot be translated (GCF_003473685.1 not present in 226)
        # Check conversion to family nodes
        self.assertCountEqual(
            res["target"]["target"],
            [
                "s__Neisseria meningitidis",
                "s__Aquabacterium_A sp001770815",
                "s__Aquabacterium_A sp002198735",
            ],
        )

    def test_convert_taxonomy_ncbi_ncbi(self):
        """
        ganon build-custom with
            --taxonomy ncbi
            --convert-taxonomy ncbi-latest
            --level class
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_convert_taxonomy_ncbi_ncbi"
        params["input"] = None
        params["input_file"] = data_dir + "build-custom/convert/convert_ncbi.tsv"
        params["input_target"] = "sequence"
        params["skip_genome_size"] = True
        params["level"] = "class"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = [
            data_dir + "build-custom/convert/convert_nodes.dmp",
        ]
        params["convert_taxonomy"] = "ncbi-latest"
        params["convert_taxonomy_files"] = [
            data_dir + "build-custom/convert/convert_nodes.dmp",
        ]

        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        # Check conversion to family nodes
        self.assertCountEqual(
            res["target"]["target"], ["28216", "28216", "186801", "28216"]
        )

    def test_input_target_file(self):
        """
        ganon build-custom with --input-target file
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_target_file"
        params["input_target"] = "file"
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        files = list_files_folder(params["input"], ext="fna.gz")
        self.assertTrue(
            res["target"]["file"].isin(files).all(), "Files missing from target"
        )
        self.assertEqual(
            len(files), res["target"].shape[0], "Wrong number of files on target"
        )
        self.assertTrue(
            res["info"]["file"].isin(files).all(), "Files missing from info"
        )
        self.assertEqual(
            len(files), res["info"].shape[0], "Wrong number of files on info"
        )

    def test_input_target_sequence(self):
        """
        ganon build-custom with --input-target sequence
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_target_sequence"
        params["input_target"] = "sequence"
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        sequences = list_sequences(list_files_folder(params["input"], ext="fna.gz"))
        self.assertTrue(
            res["target"]["target"].isin(sequences).all(), "Files missing from target"
        )
        self.assertEqual(
            len(sequences), res["target"].shape[0], "Wrong number of files on target"
        )
        self.assertTrue(
            res["info"]["target"].isin(sequences).all(), "Files missing from info"
        )
        self.assertEqual(
            len(sequences), res["info"].shape[0], "Wrong number of files on info"
        )

    def test_level_file_default(self):
        """
        ganon build-custom --input-target file and --level default
        """

        # --level default (file) - no tax
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_file_default"
        params["input_target"] = "file"
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        # --level default (file) NCBI
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_file_default_ncbi"
        params["input_target"] = "file"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["ncbi_file_info"] = data_dir + "build-custom/assembly_summary.txt"
        params["genome_size_files"] = (
            data_dir + "build-custom/species_genome_size.txt.gz"
        )
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        # --level default (file) GTDB
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_file_default_gtdb"
        params["input_target"] = "file"
        params["taxonomy"] = "gtdb"
        params["taxonomy_files"] = [
            data_dir + "build-custom/ar53_taxonomy.tsv.gz",
            data_dir + "build-custom/bac120_taxonomy.tsv.gz",
        ]
        params["genome_size_files"] = [
            data_dir + "build-custom/ar53_metadata.tsv.gz",
            data_dir + "build-custom/bac120_metadata.tsv.gz",
        ]
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

    def test_level_file_taxrank(self):
        """
        ganon build-custom --input-target file and --level {taxonomic rank}
        """

        # --level genus NCBI
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_file_genus_ncbi"
        params["input_target"] = "file"
        params["level"] = "genus"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["ncbi_file_info"] = data_dir + "build-custom/assembly_summary.txt"
        params["genome_size_files"] = (
            data_dir + "build-custom/species_genome_size.txt.gz"
        )
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(
            vars(cfg), skipped_targets=True
        )  # may skip targets without genus entry
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")
        # Tax must not have "species" (filtered out)
        self.assertFalse("species" in res["tax"]._ranks.values(), "rank found")

        # --level genus GTDB
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_file_genus_gtdb"
        params["input_target"] = "file"
        params["level"] = "genus"
        params["taxonomy"] = "gtdb"
        params["taxonomy_files"] = [
            data_dir + "build-custom/ar53_taxonomy.tsv.gz",
            data_dir + "build-custom/bac120_taxonomy.tsv.gz",
        ]
        params["genome_size_files"] = [
            data_dir + "build-custom/ar53_metadata.tsv.gz",
            data_dir + "build-custom/bac120_metadata.tsv.gz",
        ]
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")
        # Tax must not have "species" (filtered out)
        self.assertFalse("species" in res["tax"]._ranks.values(), "rank found")

    def test_level_file_leaves(self):
        """
        ganon build-custom --input-target file and --level leaves
        """

        # --level leaves NCBI
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_file_leaves_ncbi"
        params["input_target"] = "file"
        params["level"] = "leaves"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["ncbi_file_info"] = data_dir + "build-custom/assembly_summary.txt"
        params["genome_size_files"] = (
            data_dir + "build-custom/species_genome_size.txt.gz"
        )
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        self.assertIsNotNone(
            build_sanity_check_and_parse(vars(cfg)),
            "ganon build-custom sanity check failed",
        )

        # --level leaves GTDB
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_file_leaves_gtdb"
        params["input_target"] = "file"
        params["level"] = "leaves"
        params["taxonomy"] = "gtdb"
        params["taxonomy_files"] = [
            data_dir + "build-custom/ar53_taxonomy.tsv.gz",
            data_dir + "build-custom/bac120_taxonomy.tsv.gz",
        ]
        params["genome_size_files"] = [
            data_dir + "build-custom/ar53_metadata.tsv.gz",
            data_dir + "build-custom/bac120_metadata.tsv.gz",
        ]
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        self.assertIsNotNone(
            build_sanity_check_and_parse(vars(cfg)),
            "ganon build-custom sanity check failed",
        )

    def test_level_file_specialization(self):
        """
        ganon build-custom --input-target file and --level assembly/custom (specialization)
        """

        # --level assembly no tax
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_file_assembly"
        params["input_target"] = "file"
        params["level"] = "assembly"
        params["ncbi_file_info"] = data_dir + "build-custom/assembly_summary.txt"
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        self.assertIsNotNone(
            build_sanity_check_and_parse(vars(cfg)),
            "ganon build-custom sanity check failed",
        )

        # --level assembly NCBI
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_file_assembly_ncbi"
        params["input_target"] = "file"
        params["level"] = "assembly"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["ncbi_file_info"] = data_dir + "build-custom/assembly_summary.txt"
        params["genome_size_files"] = (
            data_dir + "build-custom/species_genome_size.txt.gz"
        )
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        self.assertIsNotNone(
            build_sanity_check_and_parse(vars(cfg)),
            "ganon build-custom sanity check failed",
        )

        # --level assembly GTDB
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_file_assembly_gtdb"
        params["input_target"] = "file"
        params["level"] = "assembly"
        params["taxonomy"] = "gtdb"
        params["taxonomy_files"] = [
            data_dir + "build-custom/ar53_taxonomy.tsv.gz",
            data_dir + "build-custom/bac120_taxonomy.tsv.gz",
        ]
        params["genome_size_files"] = [
            data_dir + "build-custom/ar53_metadata.tsv.gz",
            data_dir + "build-custom/bac120_metadata.tsv.gz",
        ]
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        self.assertIsNotNone(
            build_sanity_check_and_parse(vars(cfg)),
            "ganon build-custom sanity check failed",
        )

        # --level custom NCBI
        # uses info.tsv from last tests
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_file_custom_ncbi"
        params["input"] = []
        params["input_target"] = "file"
        params["input_file"] = (
            self.results_dir + "test_level_file_assembly_ncbi.info.tsv"
        )
        params["level"] = "custom"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["genome_size_files"] = (
            data_dir + "build-custom/species_genome_size.txt.gz"
        )
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        self.assertIsNotNone(
            build_sanity_check_and_parse(vars(cfg)),
            "ganon build-custom sanity check failed",
        )

        # --level custom GTDB
        # uses info.tsv from last tests
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_file_custom_gtdb"
        params["input"] = []
        params["input_target"] = "file"
        params["input_file"] = (
            self.results_dir + "test_level_file_assembly_gtdb.info.tsv"
        )
        params["level"] = "custom"
        params["taxonomy"] = "gtdb"
        params["taxonomy_files"] = [
            data_dir + "build-custom/ar53_taxonomy.tsv.gz",
            data_dir + "build-custom/bac120_taxonomy.tsv.gz",
        ]
        params["genome_size_files"] = [
            data_dir + "build-custom/ar53_metadata.tsv.gz",
            data_dir + "build-custom/bac120_metadata.tsv.gz",
        ]
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        self.assertIsNotNone(
            build_sanity_check_and_parse(vars(cfg)),
            "ganon build-custom sanity check failed",
        )

    def test_level_sequence_default(self):
        """
        ganon build-custom --input-target sequence and --level default
        """
        # --level default (sequence) - no tax
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_sequence_default"
        params["input_target"] = "sequence"
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        # --level default (sequence) - NCBI
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_sequence_default_ncbi"
        params["input_target"] = "sequence"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["ncbi_sequence_info"] = (
            data_dir + "build-custom/nucl_gb.accession2taxid.gz"
        )
        params["genome_size_files"] = (
            data_dir + "build-custom/species_genome_size.txt.gz"
        )
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

    def test_level_sequence_taxrank(self):
        """
        ganon build-custom --input-target sequence and --level {tax.rank}
        """
        # --level genus NCBI
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_sequence_genus_ncbi"
        params["input_target"] = "sequence"
        params["level"] = "genus"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["ncbi_sequence_info"] = (
            data_dir + "build-custom/nucl_gb.accession2taxid.gz"
        )
        params["genome_size_files"] = (
            data_dir + "build-custom/species_genome_size.txt.gz"
        )
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(
            vars(cfg), skipped_targets=True
        )  # may skip targets without genus entry
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")
        # Tax must not have "species" (filtered out)
        self.assertFalse("species" in res["tax"]._ranks.values(), "rank found")

    def test_level_sequence_leaves_ncbi(self):
        """
        ganon build-custom --input-target sequence and --level leaves
        """
        # --level leaves NCBI
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_level_sequence_leaves_ncbi"
        params["input_target"] = "sequence"
        params["level"] = "leaves"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["ncbi_sequence_info"] = (
            data_dir + "build-custom/nucl_gb.accession2taxid.gz"
        )
        params["genome_size_files"] = (
            data_dir + "build-custom/species_genome_size.txt.gz"
        )
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

    def test_ncbi_sequence_info(self):
        """
        ganon build-custom --ncbi-sequence-info local files
        """
        # one accession2taxid
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_ncbi_sequence_info"
        params["input_target"] = "sequence"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["ncbi_sequence_info"] = (
            data_dir + "build-custom/nucl_gb.accession2taxid.gz"
        )
        params["genome_size_files"] = (
            data_dir + "build-custom/species_genome_size.txt.gz"
        )
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        # two accession2taxid, finds all on the first, skips second
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_ncbi_sequence_info_multi"
        params["input_target"] = "sequence"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["ncbi_sequence_info"] = [
            data_dir + "build-custom/nucl_gb.accession2taxid.gz",
            data_dir + "build-custom/nucl_gb.accession2taxid.gz",
        ]
        params["genome_size_files"] = (
            data_dir + "build-custom/species_genome_size.txt.gz"
        )
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        # wrong accession2taxid
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_ncbi_sequence_info_wrong"
        params["input_target"] = "sequence"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["ncbi_sequence_info"] = [
            data_dir + "build-custom/assembly_summary.txt",
            data_dir + "build-custom/nucl_gb.accession2taxid.gz",
        ]
        params["genome_size_files"] = (
            data_dir + "build-custom/species_genome_size.txt.gz"
        )
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        # fail accession2taxid
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_ncbi_sequence_info_fail"
        params["input_target"] = "sequence"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["ncbi_sequence_info"] = (
            data_dir + "build-custom/assembly_summary.txt"
        )  # wrong, should fail
        params["genome_size_files"] = (
            data_dir + "build-custom/species_genome_size.txt.gz"
        )
        cfg = Config("build-custom", **params)
        self.assertFalse(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )

    def test_ncbi_sequence_info_download(self):
        """
        ganon build-custom --ncbi-sequence-info files simulating download
        """

        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_ncbi_sequence_info_download"
        params["input_target"] = "sequence"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"

        # Simulate download from local files (nucl_gb and species_genome_size)
        params["ncbi_url"] = (
            "file://" + os.path.abspath(data_dir) + "/build-custom/remote/"
        )
        params["ncbi_sequence_info"] = ["nucl_gb"]

        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

    def test_ncbi_file_info(self):
        """
        ganon build-custom --ncbi-file-info files
        """
        # one assembly_summary.txt
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_ncbi_file_info"
        params["input_target"] = "file"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["ncbi_file_info"] = data_dir + "build-custom/assembly_summary.txt"
        params["genome_size_files"] = (
            data_dir + "build-custom/species_genome_size.txt.gz"
        )
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        # two assembly_summary.txt finds all on the first, skips second
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_ncbi_file_info_multi"
        params["input_target"] = "file"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["ncbi_file_info"] = [
            data_dir + "build-custom/assembly_summary.txt",
            data_dir + "build-custom/assembly_summary.txt",
        ]
        params["genome_size_files"] = (
            data_dir + "build-custom/species_genome_size.txt.gz"
        )
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        # wrong assembly_summary.txt
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_ncbi_file_info_wrong"
        params["input_target"] = "file"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["ncbi_file_info"] = [
            data_dir + "build-custom/assembly_summary_empty.txt",
            data_dir + "build-custom/assembly_summary.txt",
        ]
        params["genome_size_files"] = (
            data_dir + "build-custom/species_genome_size.txt.gz"
        )
        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        # fail assembly_summary.txt
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_ncbi_file_info_fail"
        params["input_target"] = "file"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["ncbi_file_info"] = (
            data_dir + "build-custom/assembly_summary_empty.txt"
        )  # wrong, should fail
        params["genome_size_files"] = (
            data_dir + "build-custom/species_genome_size.txt.gz"
        )
        cfg = Config("build-custom", **params)
        self.assertFalse(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )

    def test_input_file_1col(self):
        """
        ganon build-custom --input-file with 1 col
        """
        files = list_files_folder(data_dir + "build-custom/files/", ext="fna.gz")

        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_file_1col"
        params["input"] = ""
        params["input_file"] = params["db_prefix"] + ".input_file.tsv"

        write_input_file(
            files,
            data_dir + "build-custom/assembly_summary.txt",
            params["input_file"],
            cols=["file"],
        )

        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        self.assertTrue(
            res["target"]["file"].isin(files).all(), "Files missing from target"
        )
        self.assertEqual(
            len(files), res["target"].shape[0], "Wrong number of files on target"
        )
        self.assertTrue(
            res["info"]["file"].isin(files).all(), "Files missing from info"
        )
        self.assertEqual(
            len(files), res["info"].shape[0], "Wrong number of files on info"
        )

    def test_input_file_2col(self):
        """
        ganon build-custom --input-file with 2 cols
        """
        files = list_files_folder(data_dir + "build-custom/files/", ext="fna.gz")

        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_file_2col"
        params["input"] = ""
        params["input_file"] = params["db_prefix"] + ".input_file.tsv"

        write_input_file(
            files,
            data_dir + "build-custom/assembly_summary.txt",
            params["input_file"],
            cols=["file", "target"],
        )

        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        self.assertTrue(
            res["target"]["file"].isin(files).all(), "Files missing from target"
        )
        self.assertEqual(
            len(files), res["target"].shape[0], "Wrong number of files on target"
        )
        self.assertTrue(
            res["info"]["file"].isin(files).all(), "Files missing from info"
        )
        self.assertEqual(
            len(files), res["info"].shape[0], "Wrong number of files on info"
        )

    def test_input_file_3col(self):
        """
        ganon build-custom --input-file with one 3 cols
        """
        files = list_files_folder(data_dir + "build-custom/files/", ext="fna.gz")

        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_file_3col"
        params["input"] = ""
        params["input_file"] = params["db_prefix"] + ".input_file.tsv"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["genome_size_files"] = (
            data_dir + "build-custom/species_genome_size.txt.gz"
        )

        write_input_file(
            files,
            data_dir + "build-custom/assembly_summary.txt",
            params["input_file"],
            cols=["file", "target", "node"],
        )

        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        # Check files
        self.assertTrue(
            res["target"]["file"].isin(files).all(), "Files missing from target"
        )
        self.assertEqual(
            len(files), res["target"].shape[0], "Wrong number of files on target"
        )
        self.assertTrue(
            res["info"]["file"].isin(files).all(), "Files missing from info"
        )
        self.assertEqual(
            len(files), res["info"].shape[0], "Wrong number of files on info"
        )
        # Check tax. level (default file)
        leaves_rank = set(map(res["tax"].rank, res["tax"].leaves()))
        self.assertEqual(len(leaves_rank), 1, "Wrong taxonomy")
        self.assertEqual(leaves_rank.pop(), "file", "Wrong taxonomy")
        self.assertEqual(
            sorted(res["tax"].leaves()), sorted(res["info"]["target"]), "Wrong taxonomy"
        )

    def test_input_file_3col_level_leaves(self):
        """
        ganon build-custom --input-file with one 3 cols with tax. --level leaves
        """
        files = list_files_folder(data_dir + "build-custom/files/", ext="fna.gz")

        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_file_3col_level_tax"
        params["input"] = ""
        params["input_file"] = params["db_prefix"] + ".input_file.tsv"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["genome_size_files"] = (
            data_dir + "build-custom/species_genome_size.txt.gz"
        )
        params["level"] = "leaves"

        write_input_file(
            files,
            data_dir + "build-custom/assembly_summary.txt",
            params["input_file"],
            cols=["file", "target", "node"],
        )

        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        self.assertTrue(
            res["target"]["file"].isin(files).all(), "Files missing from target"
        )
        self.assertEqual(
            len(files), res["target"].shape[0], "Wrong number of files on target"
        )
        self.assertTrue(
            res["info"]["file"].isin(files).all(), "Files missing from info"
        )
        self.assertEqual(
            len(files), res["info"].shape[0], "Wrong number of files on info"
        )

        # Check tax. level (node)
        self.assertEqual(
            sorted(res["tax"].leaves()), sorted(res["info"]["node"]), "Wrong taxonomy"
        )

    def test_input_file_3col_level_species(self):
        """
        ganon build-custom --input-file with one 3 cols with tax. --level species
        """
        files = list_files_folder(data_dir + "build-custom/files/", ext="fna.gz")

        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_file_3col_level_species"
        params["input"] = ""
        params["input_file"] = params["db_prefix"] + ".input_file.tsv"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["genome_size_files"] = (
            data_dir + "build-custom/species_genome_size.txt.gz"
        )
        params["level"] = "species"

        write_input_file(
            files,
            data_dir + "build-custom/assembly_summary.txt",
            params["input_file"],
            cols=["file", "target", "node"],
        )

        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        self.assertTrue(
            res["target"]["file"].isin(files).all(), "Files missing from target"
        )
        self.assertEqual(
            len(files), res["target"].shape[0], "Wrong number of files on target"
        )
        self.assertTrue(
            res["info"]["file"].isin(files).all(), "Files missing from info"
        )
        self.assertEqual(
            len(files), res["info"].shape[0], "Wrong number of files on info"
        )

        # Check tax. level (node)
        leaves_rank = set(map(res["tax"].rank, res["tax"].leaves()))
        self.assertEqual(len(leaves_rank), 1, "Wrong taxonomy")
        self.assertEqual(leaves_rank.pop(), "species", "Wrong taxonomy")

    def test_input_file_4col(self):
        """
        ganon build-custom --input-file with 4 cols
        """
        files = list_files_folder(data_dir + "build-custom/files/", ext="fna.gz")

        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_file_4col"
        params["input"] = ""
        params["input_file"] = params["db_prefix"] + ".input_file.tsv"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["genome_size_files"] = (
            data_dir + "build-custom/species_genome_size.txt.gz"
        )
        params["level"] = "custom"

        write_input_file(
            files,
            data_dir + "build-custom/assembly_summary.txt",
            params["input_file"],
            cols=["file", "target", "node", "specialization"],
        )

        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        self.assertTrue(
            res["target"]["file"].isin(files).all(), "Files missing from target"
        )
        self.assertEqual(
            len(files), res["target"].shape[0], "Wrong number of files on target"
        )
        self.assertTrue(
            res["info"]["file"].isin(files).all(), "Files missing from info"
        )
        self.assertEqual(
            len(files), res["info"].shape[0], "Wrong number of files on info"
        )
        # Check tax. level (default file)
        leaves_rank = set(map(res["tax"].rank, res["tax"].leaves()))
        self.assertEqual(len(leaves_rank), 1, "Wrong taxonomy")
        self.assertEqual(leaves_rank.pop(), "custom", "Wrong taxonomy")
        self.assertEqual(
            sorted(res["tax"].leaves()),
            sorted(res["info"]["specialization"]),
            "Wrong taxonomy",
        )
        # Specialization name is the same as specialization
        self.assertEqual(
            sorted(res["tax"].leaves()),
            sorted(map(res["tax"].name, res["tax"].leaves())),
            "Wrong taxonomy",
        )

    def test_input_file_5col(self):
        """
        ganon build-custom --input-file with 5 cols
        """
        files = list_files_folder(data_dir + "build-custom/files/", ext="fna.gz")

        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_file_5col"
        params["input"] = ""
        params["input_file"] = params["db_prefix"] + ".input_file.tsv"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["genome_size_files"] = (
            data_dir + "build-custom/species_genome_size.txt.gz"
        )
        params["level"] = "custom"

        write_input_file(
            files,
            data_dir + "build-custom/assembly_summary.txt",
            params["input_file"],
            cols=["file", "target", "node", "specialization", "specialization_name"],
        )

        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        self.assertTrue(
            res["target"]["file"].isin(files).all(), "Files missing from target"
        )
        self.assertEqual(
            len(files), res["target"].shape[0], "Wrong number of files on target"
        )
        self.assertTrue(
            res["info"]["file"].isin(files).all(), "Files missing from info"
        )
        self.assertEqual(
            len(files), res["info"].shape[0], "Wrong number of files on info"
        )
        # Check tax. level (default file)
        leaves_rank = set(map(res["tax"].rank, res["tax"].leaves()))
        self.assertEqual(len(leaves_rank), 1, "Wrong taxonomy")
        self.assertEqual(leaves_rank.pop(), "custom", "Wrong taxonomy")
        self.assertEqual(
            sorted(res["tax"].leaves()),
            sorted(res["info"]["specialization"]),
            "Wrong taxonomy",
        )
        self.assertEqual(
            sorted(map(res["tax"].name, res["tax"].leaves())),
            sorted(res["info"]["specialization_name"]),
            "Wrong taxonomy",
        )

    def test_input_file_1col_sequence(self):
        """
        ganon build-custom --input-file with one 1 col and --input-target sequence
        """
        files = list_files_folder(data_dir + "build-custom/files/", ext="fna.gz")
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_file_1col_sequence"
        params["input"] = ""
        params["input_file"] = params["db_prefix"] + ".input_file.tsv"
        params["input_target"] = "sequence"

        write_input_file(
            files,
            data_dir + "build-custom/assembly_summary.txt",
            params["input_file"],
            cols=["file"],
            sequence=True,
        )

        cfg = Config("build-custom", **params)
        # Should fail, since no sequence information was given to match the files
        self.assertFalse(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )

    def test_input_file_2col_sequence(self):
        """
        ganon build-custom --input-file with one 2 cols and --input-target sequence
        """
        files = list_files_folder(data_dir + "build-custom/files/", ext="fna.gz")
        sequences = list_sequences(files)

        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_file_2col_sequence"
        params["input"] = ""
        params["input_file"] = params["db_prefix"] + ".input_file.tsv"
        params["input_target"] = "sequence"

        write_input_file(
            files,
            data_dir + "build-custom/assembly_summary.txt",
            params["input_file"],
            cols=["file", "target"],
            sequence=True,
        )

        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        # Check files
        self.assertEqual(
            sorted(res["info"]["file"].unique()),
            sorted(files),
            "Files missing from target",
        )
        # Check sequences
        self.assertTrue(
            res["target"]["target"].isin(sequences).all(),
            "Sequences missing from target",
        )
        self.assertEqual(
            len(sequences),
            res["target"].shape[0],
            "Wrong number of sequences on target",
        )
        self.assertTrue(
            res["info"]["target"].isin(sequences).all(), "Sequences missing from info"
        )
        self.assertEqual(
            len(sequences), res["info"].shape[0], "Wrong number of sequences on info"
        )

    def test_input_file_3col_sequence(self):
        """
        ganon build-custom --input-file with one 3 cols and --input-target sequence
        """
        files = list_files_folder(data_dir + "build-custom/files/", ext="fna.gz")
        sequences = list_sequences(files)

        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_file_3col_sequence"
        params["input"] = ""
        params["input_file"] = params["db_prefix"] + ".input_file.tsv"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["genome_size_files"] = (
            data_dir + "build-custom/species_genome_size.txt.gz"
        )
        params["input_target"] = "sequence"

        write_input_file(
            files,
            data_dir + "build-custom/assembly_summary.txt",
            params["input_file"],
            cols=["file", "target", "node"],
            sequence=True,
        )

        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        # Check files
        self.assertEqual(
            sorted(res["info"]["file"].unique()),
            sorted(files),
            "Files missing from target",
        )
        # Check sequences
        self.assertTrue(
            res["target"]["target"].isin(sequences).all(),
            "Sequences missing from target",
        )
        self.assertEqual(
            len(sequences),
            res["target"].shape[0],
            "Wrong number of sequences on target",
        )
        self.assertTrue(
            res["info"]["target"].isin(sequences).all(), "Sequences missing from info"
        )
        self.assertEqual(
            len(sequences), res["info"].shape[0], "Wrong number of sequences on info"
        )
        # Check tax. level (default sequence)
        leaves_rank = set(map(res["tax"].rank, res["tax"].leaves()))
        self.assertEqual(len(leaves_rank), 1, "Wrong taxonomy")
        self.assertEqual(leaves_rank.pop(), "sequence", "Wrong taxonomy")
        self.assertEqual(
            sorted(res["tax"].leaves()), sorted(res["info"]["target"]), "Wrong taxonomy"
        )

    def test_input_file_3col_sequence_level_leaves(self):
        """
        ganon build-custom --input-file with one 3 cols and --input-target sequence and --level leaves
        """
        files = list_files_folder(data_dir + "build-custom/files/", ext="fna.gz")
        sequences = list_sequences(files)

        params = self.default_params.copy()
        params["db_prefix"] = (
            self.results_dir + "test_input_file_3col_sequence_level_leaves"
        )
        params["input"] = ""
        params["input_file"] = params["db_prefix"] + ".input_file.tsv"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["genome_size_files"] = (
            data_dir + "build-custom/species_genome_size.txt.gz"
        )
        params["input_target"] = "sequence"
        params["level"] = "leaves"

        write_input_file(
            files,
            data_dir + "build-custom/assembly_summary.txt",
            params["input_file"],
            cols=["file", "target", "node"],
            sequence=True,
        )

        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        # Check files
        self.assertEqual(
            sorted(res["info"]["file"].unique()),
            sorted(files),
            "Files missing from target",
        )
        # Check sequences
        self.assertEqual(
            sorted(res["info"]["target"]),
            sorted(sequences),
            "Sequences missing from target",
        )
        # Check tax. level (default sequence)
        self.assertEqual(
            sorted(res["tax"].leaves()),
            sorted(res["info"]["node"].unique()),
            "Wrong taxonomy",
        )

    def test_input_file_3col_sequence_level_species(self):
        """
        ganon build-custom --input-file with one 3 cols and --input-target sequence and --level species
        """
        files = list_files_folder(data_dir + "build-custom/files/", ext="fna.gz")
        sequences = list_sequences(files)

        params = self.default_params.copy()
        params["db_prefix"] = (
            self.results_dir + "test_input_file_3col_sequence_level_species"
        )
        params["input"] = ""
        params["input_file"] = params["db_prefix"] + ".input_file.tsv"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["genome_size_files"] = (
            data_dir + "build-custom/species_genome_size.txt.gz"
        )
        params["input_target"] = "sequence"
        params["level"] = "species"

        write_input_file(
            files,
            data_dir + "build-custom/assembly_summary.txt",
            params["input_file"],
            cols=["file", "target", "node"],
            sequence=True,
        )

        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        # Check files
        self.assertEqual(
            sorted(res["info"]["file"].unique()),
            sorted(files),
            "Files missing from target",
        )
        # Check sequences
        self.assertEqual(
            sorted(res["info"]["target"]),
            sorted(sequences),
            "Sequences missing from target",
        )
        # Check tax. level (species)
        leaves_rank = set(map(res["tax"].rank, res["tax"].leaves()))
        self.assertEqual(len(leaves_rank), 1, "Wrong taxonomy")
        self.assertEqual(leaves_rank.pop(), "species", "Wrong taxonomy")

    def test_input_file_4col_sequence(self):
        """
        ganon build-custom --input-file with 4 cols and --input-target sequence
        """
        files = list_files_folder(data_dir + "build-custom/files/", ext="fna.gz")
        sequences = list_sequences(files)

        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_file_4col"
        params["input"] = ""
        params["input_file"] = params["db_prefix"] + ".input_file.tsv"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["genome_size_files"] = (
            data_dir + "build-custom/species_genome_size.txt.gz"
        )
        params["input_target"] = "sequence"
        params["level"] = "custom"

        write_input_file(
            files,
            data_dir + "build-custom/assembly_summary.txt",
            params["input_file"],
            cols=["file", "target", "node", "specialization"],
            sequence=True,
        )

        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        # Check files
        self.assertEqual(
            sorted(res["info"]["file"].unique()),
            sorted(files),
            "Files missing from target",
        )
        # Check sequences
        self.assertTrue(
            res["target"]["target"].isin(res["info"]["specialization"]).all(),
            "Sequences missing from target",
        )
        self.assertTrue(
            res["info"]["target"].isin(sequences).all(), "Sequences missing from info"
        )
        self.assertEqual(
            len(sequences), res["info"].shape[0], "Wrong number of sequences on info"
        )

        # Check tax. level (default file)
        leaves_rank = set(map(res["tax"].rank, res["tax"].leaves()))
        self.assertEqual(len(leaves_rank), 1, "Wrong taxonomy")
        self.assertEqual(leaves_rank.pop(), "custom", "Wrong taxonomy")
        self.assertEqual(
            sorted(res["tax"].leaves()),
            sorted(res["info"]["specialization"].unique()),
            "Wrong taxonomy",
        )
        # Specialization name is the same as specialization
        self.assertEqual(
            sorted(res["tax"].leaves()),
            sorted(map(res["tax"].name, res["tax"].leaves())),
            "Wrong taxonomy",
        )

    def test_input_file_5col_sequence(self):
        """
        ganon build-custom --input-file with 5 cols and --input-target sequence
        """
        files = list_files_folder(data_dir + "build-custom/files/", ext="fna.gz")
        sequences = list_sequences(files)

        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_file_4col"
        params["input"] = ""
        params["input_file"] = params["db_prefix"] + ".input_file.tsv"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["genome_size_files"] = (
            data_dir + "build-custom/species_genome_size.txt.gz"
        )
        params["input_target"] = "sequence"
        params["level"] = "custom"

        write_input_file(
            files,
            data_dir + "build-custom/assembly_summary.txt",
            params["input_file"],
            cols=["file", "target", "node", "specialization", "specialization_name"],
            sequence=True,
        )

        cfg = Config("build-custom", **params)
        self.assertTrue(
            run_ganon(cfg, params["db_prefix"]), "ganon build-custom run failed"
        )
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

        # Check files
        self.assertEqual(
            sorted(res["info"]["file"].unique()),
            sorted(files),
            "Files missing from target",
        )
        # Check sequences
        self.assertTrue(
            res["target"]["target"].isin(res["info"]["specialization"]).all(),
            "Sequences missing from target",
        )
        self.assertTrue(
            res["info"]["target"].isin(sequences).all(), "Sequences missing from info"
        )
        self.assertEqual(
            len(sequences), res["info"].shape[0], "Wrong number of sequences on info"
        )

        # Check tax. level (default file)
        leaves_rank = set(map(res["tax"].rank, res["tax"].leaves()))
        self.assertEqual(len(leaves_rank), 1, "Wrong taxonomy")
        self.assertEqual(leaves_rank.pop(), "custom", "Wrong taxonomy")
        self.assertEqual(
            sorted(res["tax"].leaves()),
            sorted(res["info"]["specialization"].unique()),
            "Wrong taxonomy",
        )
        self.assertEqual(
            sorted(map(res["tax"].name, res["tax"].leaves())),
            sorted(res["info"]["specialization_name"].unique()),
            "Wrong taxonomy",
        )


if __name__ == "__main__":
    unittest.main()

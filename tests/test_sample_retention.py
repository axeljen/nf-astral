import os
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class SampleRetentionTest(unittest.TestCase):
    def test_multiple_single_and_no_retained_alignments(self):
        with tempfile.TemporaryDirectory(prefix="nf-astral-retention-") as directory:
            root = Path(directory)
            (root / "one.fasta").write_text(">ind1\nACGT\n>ind2\nACGT\n>reference\nACGT\n")
            (root / "two.fasta").write_text(">ind2\nACGT\n")
            (root / "failed.fasta").write_text("")
            (root / "one.filter_summary.txt").write_text(
                "one\t3\t4\t2\t4\tPASS\n"
            )
            (root / "two.filter_summary.txt").write_text(
                "two\t3\t4\t0\t0\tFAIL\n"
            )
            (root / "nextflow.config").write_text("process.executor = 'local'\n")
            module = (ROOT / "modules/summarize_sample_retention").as_posix()
            (root / "main.nf").write_text(f"""
nextflow.enable.dsl = 2
include {{ SUMMARIZE_SAMPLE_SURVIVAL as MULTIPLE }} from '{module}'
include {{ SUMMARIZE_SAMPLE_SURVIVAL as SINGLE }} from '{module}'
include {{ SUMMARIZE_SAMPLE_SURVIVAL as NONE }} from '{module}'
include {{ SUMMARIZE_FILTERING }} from '{ROOT / "modules/summarize_filtering"}'

workflow {{
    sample_names = Channel.of(
        tuple('ind1', 'mask1'), tuple('ind2', 'mask2'), tuple('ind3', 'mask3')
    ).toList().map {{ rows -> rows.collect {{ it[0] }} }}
    retained = Channel.of(
        tuple('one', file('one.fasta')),
        tuple('two', file('two.fasta')),
        tuple('failed', file('failed.fasta'))
    ).filter {{ window_id, fasta -> fasta.size() > 0 }}
     .map {{ window_id, fasta -> fasta }}
     .toList()
    single = Channel.of(tuple('one', file('one.fasta')))
        .filter {{ window_id, fasta -> fasta.size() > 0 }}
        .map {{ window_id, fasta -> fasta }}
        .toList()
    none = Channel.of(tuple('failed', file('failed.fasta')))
        .filter {{ window_id, fasta -> fasta.size() > 0 }}
        .map {{ window_id, fasta -> fasta }}
        .toList()

    MULTIPLE(sample_names, retained)
    SINGLE(sample_names, single)
    NONE(sample_names, none)
    MULTIPLE.out.summary.subscribe {{
        assert it.name == 'sample_survival_summary.tsv'
        it.copyTo(file('multiple.tsv'))
    }}
    SINGLE.out.summary.subscribe {{
        assert it.name == 'sample_survival_summary.tsv'
        it.copyTo(file('single.tsv'))
    }}
    NONE.out.summary.subscribe {{
        assert it.name == 'sample_survival_summary.tsv'
        it.copyTo(file('none.tsv'))
    }}
    SUMMARIZE_FILTERING(Channel.of(
        file('two.filter_summary.txt'), file('one.filter_summary.txt')
    ).collect())
    SUMMARIZE_FILTERING.out.summary.subscribe {{
        assert it.name == 'alignment_filtering_summary.tsv'
        it.copyTo(file('filtering.tsv'))
    }}
}}
""")
            environment = dict(os.environ, NXF_OFFLINE="true", NXF_ANSI_LOG="false")
            result = subprocess.run(
                ["nextflow", "-C", str(root / "nextflow.config"), "run", "main.nf"],
                cwd=root,
                env=environment,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                timeout=180,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            header = "sample\tnumber_of_trees\n"
            for name, counts in [
                ("multiple", (1, 2, 0)),
                ("single", (1, 1, 0)),
                ("none", (0, 0, 0)),
            ]:
                with self.subTest(case=name):
                    expected = header + "".join(
                        f"ind{index}\t{count}\n"
                        for index, count in enumerate(counts, start=1)
                    )
                    self.assertEqual((root / f"{name}.tsv").read_text(), expected)
            self.assertEqual(
                (root / "filtering.tsv").read_text(),
                "window_id\tseq_count_orig\tsite_count_orig\t"
                "seq_count_retained\tsite_count_retained\ttree_generated\n"
                "one\t3\t4\t2\t4\tPASS\n"
                "two\t3\t4\t0\t0\tFAIL\n",
            )

    def test_workflow_imports(self):
        with tempfile.TemporaryDirectory(prefix="nf-astral-imports-") as directory:
            root = Path(directory)
            config = root / "nextflow.config"
            config.write_text("params.reference = null\n")
            result = subprocess.run(
                ["nextflow", "-C", str(config), "run", str(ROOT / "main.nf")],
                cwd=root,
                env=dict(os.environ, NXF_OFFLINE="true", NXF_ANSI_LOG="false"),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                timeout=180,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(
                "Please provide --reference <path to reference fasta>",
                result.stdout + result.stderr,
            )


if __name__ == "__main__":
    unittest.main()

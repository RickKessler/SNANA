"""Run with: python -m unittest discover -s util/submit_batch/tests -v.

Uses temporary worker logs and mocked SLURM commands; requires PyYAML.
No batch jobs are submitted or cancelled.
"""

import builtins
import logging
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch


SNANA_ROOT = Path(__file__).resolve().parents[3]
os.environ.setdefault("SNANA_DIR", str(SNANA_ROOT))
os.environ.setdefault("SHELL", "/bin/bash")
sys.path.insert(0, str(SNANA_ROOT / "util" / "submit_batch"))

import submit_util as util
from submit_prog_base import Program


def accounting_result(stdout="", stderr="", returncode=0):
    return subprocess.CompletedProcess(["sacct"], returncode, stdout, stderr)


class GrepEvidenceTests(unittest.TestCase):
    def test_counts_and_original_line_numbers_without_duplicate_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            path = str(Path(directory) / "CPU0000.LOG")
            Path(path).write_bytes(
                b"started\n_@_ copied ERROR\nTMPDIR error ignored\n"
                b"slurmstepd: error: KILL memory limit\nerror: bad byte \xff\n")
            evidence = {}
            result = util.grep([path], ["error", "KILL", "ERROR"], ["TMPDIR"],
                               False, evidence)
            self.assertEqual(result, (3, {path: 3}, {"error": 2, "KILL": 1, "ERROR": 0}))
            self.assertEqual(evidence[path], [
                (4, "slurmstepd: error: KILL memory limit"),
                (5, "error: bad byte \ufffd")])
            # Existing callers still receive the same three return values.
            self.assertEqual(util.grep([path], ["error", "KILL", "ERROR"],
                                       ["TMPDIR"], False), result)

    def test_evidence_is_bounded_but_all_matches_are_counted(self):
        with tempfile.TemporaryDirectory() as directory:
            path = str(Path(directory) / "CPU0000.LOG")
            Path(path).write_text(("error " + "x" * 3000 + "\n") * 25)
            evidence = {}
            total, _, _ = util.grep([path], ["error"], [], False, evidence)
            self.assertEqual(total, 25)
            self.assertEqual(len(evidence[path]), 20)
            self.assertIn("line truncated", evidence[path][0][1])
            self.assertLess(len(evidence[path][0][1]), 2100)


class AccountingTests(unittest.TestCase):
    @patch("submit_util.subprocess.run")
    def test_single_query_keeps_batch_step_failure_and_full_state(self, run):
        run.return_value = accounting_result(
            "101|RUNNING|0:0|00:01:00|01:00:00||4Gn|node1\n"
            "101.batch|OUT_OF_MEMORY|0:9|00:00:59||4096M||node1\n"
            "102|CANCELLED by 1234|0:15|00:00:10|01:00:00||4Gn|node2\n"
            "999|FAILED|1:0|||||\n"
            "malformed row\n")
        records, note = util.get_slurm_accounting([101, "101", 102])
        self.assertEqual(note, "")
        self.assertEqual(len(records["101"]), 2)
        self.assertEqual(records["101"][1]["State"], "OUT_OF_MEMORY")
        self.assertEqual(records["102"][0]["State"], "CANCELLED by 1234")
        self.assertNotIn("999", records)
        run.assert_called_once()
        command = run.call_args.args[0]
        self.assertIn("--jobs=101,102", command)
        self.assertNotIn("--allocations", command)
        self.assertNotIn("Reason", " ".join(command))
        self.assertEqual(run.call_args.kwargs["timeout"], 10)
        self.assertFalse(run.call_args.kwargs.get("shell", False))

    @patch("submit_util.subprocess.run")
    def test_unavailable_accounting_returns_explanation(self, run):
        for exception, expected in [
            (FileNotFoundError("sacct missing"), "sacct missing"),
            (subprocess.TimeoutExpired("sacct", 10), "timed out after 10 seconds"),
            (PermissionError("sacct denied"), "sacct denied"),
        ]:
            with self.subTest(exception=exception):
                run.side_effect = exception
                records, note = util.get_slurm_accounting([101])
                self.assertEqual(records, {"101": []})
                self.assertIn(expected, note)
        run.side_effect = None
        run.return_value = accounting_result(stderr="database\nnot available", returncode=1)
        _, note = util.get_slurm_accounting([101])
        self.assertIn("sacct exit 1; database not available", note)

    @patch("submit_util.subprocess.run")
    def test_missing_or_invalid_ids_do_not_launch_command(self, run):
        self.assertIn("ID unavailable", util.get_slurm_accounting([])[1])
        self.assertIn("invalid job ID", util.get_slurm_accounting(["101; echo hello"])[1])
        run.assert_not_called()


class FailureReportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="snana-slurm-")
        self.addCleanup(self.temp.cleanup)
        self.output = Path(self.temp.name) / "output.with.dots"
        self.scripts = self.output / "SPLIT JOBS"
        self.scripts.mkdir(parents=True)
        self.failed = self.scripts / "CPU0000_JOBLIST_node.LOG"
        self.healthy = self.scripts / "CPU0001_JOBLIST_node.LOG"
        self.done = self.scripts / "CPU0002_JOBLIST_node.LOG"
        self.failed.write_text("started\nslurmstepd: error: Detected 1 oom_kill event\n")
        self.healthy.write_text("started\n")
        self.done.write_text("ERROR from a finished worker is not inspected\n")
        self.done.with_suffix(".DONE").touch()

        self.program = object.__new__(Program)
        self.program.config_yaml = {
            "args": SimpleNamespace(devel_flag=0, kill=False),
            "CONFIG": {"BATCH_INFO": "sbatch template 3"},
        }
        self.program.config_prep = {
            "output_dir": str(self.output),
            "submit_info_yaml": {
                "SCRIPT_DIR": str(self.scripts),
                "SUBMIT_MODE": "BATCH",
                "SBATCH_LIST": [[0, 101, "test-CPU0000"],
                                [1, 102, "test-CPU0001"],
                                [2, 103, "test-CPU0002"]],
                "DONE_STAMP_LIST": ["ALL.DONE"],
            },
        }
        self.program.create_output_dir = Mock()
        # Exercise the real kill_jobs report/stamp writes, but never scancel.
        self.program.kill_sbatch_jobs = Mock()
        self.program.kill_ssh_jobs = Mock()
        self.sleep = patch("submit_prog_base.time.sleep").start()
        self.addCleanup(patch.stopall)

    @patch("submit_util.subprocess.run")
    def test_report_persisted_before_cancellation_and_peers_point_to_culprit(self, run):
        run.return_value = accounting_result(
            "101|RUNNING|0:0|00:01:00|01:00:00||4Gn|node1\n"
            "101.batch|OUT_OF_MEMORY|0:9|00:00:59||4096M||node1\n")

        def at_cancellation():
            report = (self.output / "MERGE.LOG").read_text()
            self.assertIn("CPU0000_JOBLIST_node.LOG:2: slurmstepd: error: Detected 1 oom_kill event", report)
            self.assertIn("SLURM job ID: 101", report)
            self.assertIn("JobIDRaw=101.batch State=OUT_OF_MEMORY ExitCode=0:9", report)
            self.assertNotIn(self.done.name, report)
            self.assertNotIn("slurm failure identified in " + self.healthy.name, report)
            self.assertEqual((self.output / "ALL.DONE").read_text().strip(), "FAIL")
            self.assertIn("MERGE.LOG", self.healthy.read_text())
            self.assertIn(self.failed.name, self.healthy.read_text())

        self.program.kill_sbatch_jobs.side_effect = at_cancellation
        with self.assertLogs(level=logging.INFO) as captured:
            self.program.check_for_slurm_failure()
        self.program.kill_sbatch_jobs.assert_called_once()
        self.assertTrue(self.program.config_yaml["args"].kill)
        self.assertIn("--jobs=101", run.call_args.args[0])
        # Report text copied into the merge worker's log is excluded on rescan.
        report_lines = [entry for entry in captured.output if "_@_" in entry]
        with self.healthy.open("a") as stream:
            stream.write("\n".join(report_lines) + "\n")
        self.assertEqual(util.grep([str(self.healthy)], ["FAIL", "KILL", "ERROR",
                                      "CANCEL", "fail", "kill", "error"], [], False)[0], 0)

    @patch("submit_util.subprocess.run")
    def test_accounting_failure_does_not_prevent_report_or_shutdown(self, run):
        for exception in [FileNotFoundError("no sacct"), subprocess.TimeoutExpired("sacct", 10)]:
            with self.subTest(exception=exception):
                run.side_effect = exception
                self.program.check_for_slurm_failure()
                report = (self.output / "MERGE.LOG").read_text()
                self.assertIn("Detected 1 oom_kill event", report)
                self.assertIn("SLURM accounting unavailable", report)
                self.assertEqual((self.output / "ALL.DONE").read_text().strip(), "FAIL")
        self.assertEqual(self.program.kill_sbatch_jobs.call_count, 2)

    @patch("submit_util.subprocess.run")
    def test_empty_accounting_is_not_misrepresented_as_a_diagnosis(self, run):
        run.return_value = accounting_result()
        self.program.check_for_slurm_failure()
        report = (self.output / "MERGE.LOG").read_text()
        self.assertIn("No SLURM accounting record available yet", report)
        self.assertNotIn("State=OUT_OF_MEMORY", report)
        self.program.kill_sbatch_jobs.assert_called_once()

    @patch("submit_util.subprocess.run")
    def test_multiple_failing_cpus_keep_their_own_accounting(self, run):
        self.healthy.write_text("slurmstepd: error: TIME LIMIT reached\n")
        run.return_value = accounting_result(
            "101.batch|OUT_OF_MEMORY|0:9|00:00:59||4096M||node1\n"
            "102|TIMEOUT|0:15|01:00:00|01:00:00||4Gn|node2\n")
        self.program.check_for_slurm_failure()
        run.assert_called_once()
        self.assertIn("--jobs=101,102", run.call_args.args[0])
        report = (self.output / "MERGE.LOG").read_text()
        first, second = report.split("slurm failure identified in " + self.healthy.name)
        self.assertIn("State=OUT_OF_MEMORY", first)
        self.assertNotIn("State=TIMEOUT", first)
        self.assertIn("SLURM job ID: 102", second)
        self.assertIn("State=TIMEOUT", second)
        self.assertNotIn("State=OUT_OF_MEMORY", second)

    @patch("submit_util.subprocess.run")
    def test_unwritable_peer_log_does_not_prevent_shutdown(self, run):
        run.return_value = accounting_result()
        real_open = builtins.open

        def open_with_unwritable_peer(path, mode="r", *args, **kwargs):
            if str(path) == str(self.healthy) and mode == "at":
                raise PermissionError("peer log is read-only")
            return real_open(path, mode, *args, **kwargs)

        with patch("builtins.open", side_effect=open_with_unwritable_peer):
            with self.assertLogs(level=logging.WARNING) as captured:
                self.program.check_for_slurm_failure()
        self.assertIn("_@_ Could not append cancellation notice", "\n".join(captured.output))
        self.assertIn("Detected 1 oom_kill event", (self.output / "MERGE.LOG").read_text())
        self.program.kill_sbatch_jobs.assert_called_once()

    @patch("submit_util.subprocess.run")
    def test_no_failure_never_queries_accounting_or_cancels(self, run):
        self.failed.write_text("running normally\nTMPDIR error ignored\n_@_ copied ERROR\n")
        self.program.check_for_slurm_failure()
        run.assert_not_called()
        self.program.kill_sbatch_jobs.assert_not_called()
        self.assertFalse((self.output / "MERGE.LOG").exists())

    @patch("submit_util.subprocess.run")
    def test_missing_submit_ids_still_reports_log_evidence(self, run):
        del self.program.config_prep["submit_info_yaml"]["SBATCH_LIST"]
        self.program.check_for_slurm_failure()
        report = (self.output / "MERGE.LOG").read_text()
        self.assertIn("SLURM job ID unavailable in SUBMIT.INFO", report)
        self.assertIn("Detected 1 oom_kill event", report)
        run.assert_not_called()
        self.program.kill_sbatch_jobs.assert_called_once()

    @patch("submit_util.subprocess.run")
    def test_non_slurm_report_does_not_query_slurm(self, run):
        for submit_mode, batch_info in [("SSH", ""), ("BATCH", "qsub template 3")]:
            with self.subTest(submit_mode=submit_mode, batch_info=batch_info):
                self.program.config_prep["submit_info_yaml"]["SUBMIT_MODE"] = submit_mode
                self.program.config_yaml["CONFIG"]["BATCH_INFO"] = batch_info
                report = self.program.slurm_failure_report(
                    [str(self.failed)], {str(self.failed): [(2, "original error")]})
                self.assertIn("original error", report)
                self.assertIn("not applicable", report)
        run.assert_not_called()


if __name__ == "__main__":
    unittest.main()

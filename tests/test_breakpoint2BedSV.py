import subprocess
from pathlib import Path
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "bin" / "breakpoint2BedSV.py"
BASE = Path("tests/data")
CASES = [
    "test_01_angle-bracketed_notation",
    "test_02_square-bracketed_notation",
    "test_03_sequence_notation",
    "test_04_TRA_without_chr",
    "test_05_bcf",
    "test_06_compressedInputVCF",
    "test_07_single_breakend_with_redundancy",
    "test_08_SNV_with_redundancy",
]
ERROR_CASES = [
    "test_09_no_vcf_header_line",
]


def find_input(case_dir):
    for ext in ["vcf", "vcf.gz", "bcf"]:
        f = case_dir / "input" / f"test.{ext}"
        if f.exists():
            return f
    raise FileNotFoundError("[ERROR] No test.vcf or test.bcf found")


def run_case(tmp_path, case):
    case_dir = BASE / case

    input_vcf = find_input(case_dir)
    expected = case_dir / "expected" / "test.bed"

    out_file = tmp_path / "out.bed"

    result=subprocess.run(
        [
            "python3",
            SCRIPT,
            "-i", str(input_vcf),
            "-o", str(out_file),
        ],
		capture_output=True,
        text=True,
    )

    return result, out_file, expected


@pytest.mark.parametrize("case", CASES)
def test_breakpoint2BedSV(tmp_path, case):
    result, out, exp = run_case(tmp_path, case)
    assert result.returncode == 0
    assert out.read_text() == exp.read_text()


@pytest.mark.parametrize("case", ERROR_CASES)
def test_breakpoint2BedSV_errors(tmp_path, case):
    result, _, _ = run_case(tmp_path, case)
    assert result.returncode != 0
    assert "ERROR" in result.stderr




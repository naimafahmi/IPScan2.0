from pathlib import Path

import pytest

from ipscan import build_parser, discover_samples, load_config, validate_args


def test_discover_samples_prefers_bams_over_coverage_directories(tmp_path):
    (tmp_path / "sample1.bam").touch()
    coverage_dir = tmp_path / "sample2"
    coverage_dir.mkdir()
    (coverage_dir / "chr1.txt").touch()

    assert discover_samples(tmp_path) == [tmp_path / "sample1.bam"]


def test_discover_samples_finds_coverage_directories(tmp_path):
    coverage_dir = tmp_path / "sample1"
    coverage_dir.mkdir()
    (coverage_dir / "chr1.txt").touch()

    assert discover_samples(tmp_path) == [coverage_dir]


def test_single_mode_rejects_second_input_directory(tmp_path):
    annotation = tmp_path / "annotation.tsv"
    annotation.touch()
    group1 = tmp_path / "group1"
    group2 = tmp_path / "group2"
    group1.mkdir()
    group2.mkdir()
    (group1 / "sample.bam").touch()
    (group2 / "sample.bam").touch()

    args = build_parser().parse_args([
        "--mode", "single",
        "--species", "hg38",
        "--input-dir1", str(group1),
        "--input-dir2", str(group2),
        "--annotation", str(annotation),
        "--output-dir", str(tmp_path / "output"),
        "--log-dir", str(tmp_path / "logs"),
    ])

    with pytest.raises(ValueError, match="Single mode accepts only one"):
        validate_args(args)


def test_differential_mode_requires_second_input_directory(tmp_path):
    annotation = tmp_path / "annotation.tsv"
    annotation.touch()
    group1 = tmp_path / "group1"
    group1.mkdir()
    (group1 / "sample.bam").touch()

    args = build_parser().parse_args([
        "--mode", "differential",
        "--species", "hg38",
        "--input-dir1", str(group1),
        "--annotation", str(annotation),
        "--output-dir", str(tmp_path / "output"),
        "--log-dir", str(tmp_path / "logs"),
    ])

    with pytest.raises(ValueError, match="requires --input-dir2"):
        validate_args(args)


def test_three_prime_directory_is_optional_but_must_exist(tmp_path):
    annotation = tmp_path / "annotation.tsv"
    annotation.touch()
    group1 = tmp_path / "group1"
    group1.mkdir()
    (group1 / "sample.bam").touch()

    args = build_parser().parse_args([
        "--mode", "single",
        "--species", "hg38",
        "--input-dir1", str(group1),
        "--annotation", str(annotation),
        "--output-dir", str(tmp_path / "output"),
        "--log-dir", str(tmp_path / "logs"),
    ])
    validate_args(args)
    assert args.three_prime_dir is None


def test_peptides_flag_is_optional():
    parser = build_parser()
    args = parser.parse_args([
        "--mode", "single",
        "--species", "hg38",
        "--input-dir1", "group1",
        "--annotation", "annotation.tsv",
        "--output-dir", "output",
        "--log-dir", "logs",
        "--peptides",
    ])

    assert args.peptides is True


def test_toml_config_populates_parser_and_cli_overrides(tmp_path):
    config_path = tmp_path / "run.toml"
    config_path.write_text(
        'mode = "single"\n'
        'species = "hg38"\n'
        'input_dir1 = "from-config"\n'
        'peptides = true\n'
    )

    config = load_config(config_path)
    args = build_parser(config).parse_args([
        "--config", str(config_path),
        "--input-dir1", "from-cli",
        "--no-peptides",
    ])

    assert args.input_dir1 == "from-cli"
    assert args.peptides is False

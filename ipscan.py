import argparse
import logging
import os
import time
from pathlib import Path
import sys

try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10 and earlier
    import tomli as tomllib

from methods.print_and_log import print_and_log
from methods.methods import generate_ipscan_events, get_chromosomes, prepare_annotation_context
from methods.peaks import generate_three_prime_peaks
from methods.peptides import generate_peptide_sequences


def load_config(config_path):
    with open(config_path, "rb") as config_file:
        config = tomllib.load(config_file)
    if not isinstance(config, dict):
        raise ValueError("Configuration file must contain a TOML table")
    return config


def build_parser(config=None):
    config = config or {}
    parser = argparse.ArgumentParser(description="Detect intronic polyadenylation events.")
    parser.add_argument("--config", help="Optional TOML configuration file.")
    parser.add_argument("--mode", choices=("single", "differential"), default=config.get("mode"))
    parser.add_argument(
        "-species", "--species",
        choices=("hg19", "hg38", "mm10", "mm39"),
        default=config.get("species"),
    )
    parser.add_argument(
        "-s1", "--sample1", "--input-dir1",
        dest="input_dir1", default=config.get("input_dir1"),
        help="Directory containing BAM samples or sample coverage directories.",
    )
    parser.add_argument(
        "-s2", "--sample2", "--input-dir2",
        dest="input_dir2", default=config.get("input_dir2"),
        help="Second group directory; required for differential mode.",
    )
    parser.add_argument(
        "-anno", "--annotation", default=config.get("annotation"),
        help="Tab-separated gene annotation file.",
    )
    parser.add_argument("-out", "--output-dir", dest="output_dir", default=config.get("output_dir"))
    parser.add_argument("-log", "--log-dir", dest="log_dir", default=config.get("log_dir"))
    parser.add_argument(
        "--three-prime-dir", default=config.get("three_prime_dir"),
        help="Optional directory containing 3'-end-seq BAM files for peak generation.",
    )
    peptide_options = parser.add_mutually_exclusive_group()
    peptide_options.add_argument(
        "--peptides", dest="peptides", action="store_true",
        help="Generate peptide sequences when supported.",
    )
    peptide_options.add_argument(
        "--no-peptides", dest="peptides", action="store_false",
        help="Disable peptide sequence generation.",
    )
    parser.set_defaults(peptides=config.get("peptides", False))
    return parser


def configure_logging(log_dir):
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, "ipscan.log")
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    return logging.getLogger("ipscan")


def discover_samples(input_dir):
    """Return BAM files or coverage directories contained in a group directory."""
    group_dir = Path(input_dir)
    if not group_dir.is_dir():
        raise NotADirectoryError("Input group directory not found: " + str(group_dir))

    if any(group_dir.glob("chr*.txt")):
        return [group_dir]

    samples = sorted(group_dir.glob("*.bam"))
    if not samples:
        samples = sorted(
            sample_dir
            for sample_dir in group_dir.iterdir()
            if sample_dir.is_dir() and any(sample_dir.glob("chr*.txt"))
        )
    if not samples:
        raise FileNotFoundError(
            "No BAM files or chromosome coverage directories found in: " + str(group_dir)
        )
    return samples


def validate_args(args):
    required_values = {
        "--mode": args.mode,
        "--species": args.species,
        "--input-dir1": args.input_dir1,
        "--annotation": args.annotation,
        "--output-dir": args.output_dir,
        "--log-dir": args.log_dir,
    }
    missing_values = [name for name, value in required_values.items() if not value]
    if missing_values:
        raise ValueError("Missing required arguments: " + ", ".join(missing_values))
    if not os.path.isfile(args.annotation):
        raise FileNotFoundError("Annotation file not found: " + args.annotation)
    if not os.path.isdir(args.input_dir1):
        raise NotADirectoryError("Input group 1 directory not found: " + args.input_dir1)
    if args.mode == "differential" and not args.input_dir2:
        raise ValueError("Differential mode requires --input-dir2")
    if args.input_dir2 and not os.path.isdir(args.input_dir2):
        raise NotADirectoryError("Input group 2 directory not found: " + args.input_dir2)
    if args.mode == "single" and args.input_dir2:
        raise ValueError("Single mode accepts only one input directory")
    if args.three_prime_dir and not os.path.isdir(args.three_prime_dir):
        raise NotADirectoryError("3'-end-seq directory not found: " + args.three_prime_dir)
    discover_samples(args.input_dir1)
    if args.mode == "differential":
        discover_samples(args.input_dir2)
    os.makedirs(args.output_dir, exist_ok=True)


def run_group(
    species,
    annotation,
    input_dir,
    output_dir,
    logger,
    group_name,
    peak_file=None,
    generate_peptides=False,
):
    samples = discover_samples(input_dir)
    print_and_log(f"{group_name}: {len(samples)} sample(s)", logger)
    annotation_context = prepare_annotation_context(species, annotation, logger, peak_file=peak_file)
    results = []
    for sample in samples:
        print_and_log("Processing sample: " + str(sample), logger)
        outputs = generate_ipscan_events(
            species,
            annotation,
            str(sample),
            output_dir,
            logger,
            annotation_context=annotation_context,
        )
        if generate_peptides:
            outputs["peptides"] = generate_peptide_sequences(
                outputs["filtered"],
                annotation,
                output_dir,
                logger,
            )
        results.append({"sample": str(sample), "outputs": outputs})
    return results


def compare_groups(group1_results, group2_results, output_dir, logger):
    """Comparison boundary for differential analysis; metric is intentionally pending."""
    print_and_log("Differential comparison is not implemented yet; group results are available.", logger)
    return None


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    config_path = None
    for index, value in enumerate(argv):
        if value == "--config" and index + 1 < len(argv):
            config_path = argv[index + 1]
            break
        if value.startswith("--config="):
            config_path = value.split("=", 1)[1]
            break
    config = {}
    if config_path:
        config = load_config(config_path)
    args = build_parser(config).parse_args(argv)
    validate_args(args)
    logger = configure_logging(args.log_dir)
    start = time.perf_counter()

    print_and_log("Beginning IPScan run", logger)
    print_and_log("Mode: " + args.mode, logger)
    peak_file = None
    if args.three_prime_dir:
        peak_file = os.path.join(
            Path(__file__).resolve().parent,
            "anno",
            f"{args.species}_polyA_sites_with_3end.csv",
        )
        generate_three_prime_peaks(
            args.species,
            args.annotation,
            args.three_prime_dir,
            peak_file,
            os.path.join(
                Path(__file__).resolve().parent,
                "anno",
                f"{args.species}_polyA_sites.csv",
            ),
            get_chromosomes(args.species),
        )
        print_and_log("Using merged PolyA_DB and 3'-end-seq peaks: " + peak_file, logger)
    group1_results = run_group(
        args.species,
        args.annotation,
        args.input_dir1,
        args.output_dir,
        logger,
        "Group 1",
        peak_file,
        args.peptides,
    )
    if args.mode == "differential":
        group2_output_dir = os.path.join(args.output_dir, "group2")
        group2_results = run_group(
            args.species,
            args.annotation,
            args.input_dir2,
            group2_output_dir,
            logger,
            "Group 2",
            peak_file,
            args.peptides,
        )
        compare_groups(group1_results, group2_results, args.output_dir, logger)
    print_and_log("IPScan completed in %.2f seconds" % (time.perf_counter() - start), logger)


if __name__ == "__main__":
    main()


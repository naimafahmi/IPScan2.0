"""Coverage table utilities used by the IPScan event detector."""

import os
from bisect import bisect_left

import pandas as pd
import pysam


def bam_to_coverage_tables(bam_path, output_dir, chromosomes):
    """Write legacy chromosome coverage tables from an indexed BAM file.

    Each output row contains 1-based genomic position and read depth. The
    chromosome is represented by the output filename.
    """
    bam_path = os.fspath(bam_path)
    os.makedirs(output_dir, exist_ok=True)
    expected_files = [os.path.join(output_dir, chromosome + ".txt") for chromosome in chromosomes]
    if all(os.path.isfile(path) for path in expected_files):
        return False

    bai_paths = (bam_path + ".bai", os.path.splitext(bam_path)[0] + ".bai")
    if not any(os.path.isfile(path) for path in bai_paths):
        pysam.index(bam_path)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        references = set(bam.references)
        for chromosome in chromosomes:
            output_path = os.path.join(output_dir, chromosome + ".txt")
            with open(output_path, "w", encoding="ascii") as output_file:
                if chromosome not in references:
                    continue
                for column in bam.pileup(
                    chromosome,
                    truncate=True,
                    stepper="all",
                    min_base_quality=0,
                ):
                    output_file.write(f"{column.reference_pos + 1}\t{column.nsegments}\n")

    return True


def load_coverage_table(coverage_path):
    """Load and validate a two-column chromosome coverage table."""
    try:
        coverage = pd.read_csv(coverage_path, sep="\t", header=None)
    except pd.errors.EmptyDataError:
        return None

    if coverage.shape[1] != 2:
        raise ValueError(
            f"Coverage file must contain exactly 2 tab-separated columns: {coverage_path}"
        )

    coverage = coverage.apply(pd.to_numeric, errors="coerce")
    if coverage.isna().any().any():
        raise ValueError(f"Coverage file contains non-numeric values: {coverage_path}")

    positions = coverage.iloc[:, 0]
    read_counts = coverage.iloc[:, 1]
    if not (positions == positions.astype("int64")).all():
        raise ValueError(f"Coverage positions must be integers: {coverage_path}")
    if not (read_counts == read_counts.astype("int64")).all():
        raise ValueError(f"Coverage values must be integers: {coverage_path}")
    if (read_counts < 0).any():
        raise ValueError(f"Coverage values must be non-negative: {coverage_path}")
    if not positions.is_monotonic_increasing or positions.duplicated().any():
        raise ValueError(
            f"Coverage positions must be sorted and unique: {coverage_path}"
        )

    return coverage.astype("int64")


class CoverageIndex:
    def __init__(self, coverage_rows, positions):
        self.coverage_rows = coverage_rows
        self.positions = positions
        self.prefix_totals = [0]
        for row in coverage_rows:
            self.prefix_totals.append(self.prefix_totals[-1] + int(row[1]))

    def sum(self, start, end):
        if end < start or not self.coverage_rows:
            return 0, 0

        first = bisect_left(self.positions, start)
        last = bisect_left(self.positions, end)
        if last == len(self.coverage_rows):
            last -= 1
        elif int(self.coverage_rows[last][0]) != end:
            last -= 1

        total = 0
        if first <= last:
            total = self.prefix_totals[last + 1] - self.prefix_totals[first]
        return total, end - start + 1


def count_total_read_count(start, end, coverage_rows, positions):
    """Sum coverage over an inclusive interval.

    ``coverage_rows`` must contain position in column zero and coverage in
    column one. Positions must be sorted in ascending order.
    """
    if end < start or not coverage_rows:
        return 0, 0

    first = bisect_left(positions, start)
    last = bisect_left(positions, end)
    if last < len(coverage_rows) and int(coverage_rows[last][0]) != end:
        last -= 1

    total = sum(int(coverage_rows[index][1]) for index in range(first, last + 1))
    return total, end - start + 1

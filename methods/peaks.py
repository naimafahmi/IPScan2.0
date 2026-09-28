"""Polyadenylation peak generation and database merging."""

from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
import pysam

from methods.annotation import get_intron_list, make_full_dictionary, merge_chromosome_dictionary


def discover_three_prime_bams(input_dir):
    paths = sorted(Path(input_dir).glob("*.bam"))
    if not paths:
        raise FileNotFoundError("No 3'-end-seq BAM files found in: " + str(input_dir))
    return paths


def _read_three_prime_counts(bam_paths):
    counts = Counter()
    for bam_path in bam_paths:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if read.reference_end is None:
                    continue
                if read.is_reverse:
                    position = read.reference_start + 1
                else:
                    position = read.reference_end
                counts[(read.reference_name, position)] += 1
    return counts


def _cluster_positions(counts, merge_distance, min_support):
    selected = sorted(position for position, count in counts.items() if count >= min_support)
    clusters = []
    for position in selected:
        if not clusters or position - clusters[-1][-1] > merge_distance:
            clusters.append([position])
        else:
            clusters[-1].append(position)
    return [max(cluster, key=lambda position: counts[position]) for cluster in clusters]


def generate_three_prime_peaks(
    species,
    annotation_file,
    three_prime_dir,
    output_path,
    polyadb_path,
    chromosomes,
    min_support=3,
    merge_distance=25,
):
    """Generate intronic 3'-end peaks and merge them with PolyA_DB peaks."""
    bam_paths = discover_three_prime_bams(three_prime_dir)
    annotation = pd.read_csv(annotation_file, sep="\t")
    chromosome_dict = make_full_dictionary(annotation, chromosomes)
    merged_dict = merge_chromosome_dictionary(chromosome_dict, chromosomes)

    introns_by_chromosome = defaultdict(list)
    for chromosome in chromosomes:
        for gene, (exons, strand) in merged_dict[chromosome].items():
            for start, end in get_intron_list(exons):
                introns_by_chromosome[chromosome].append((start, end, gene, strand))

    generated = defaultdict(Counter)
    for (chromosome, position), count in _read_three_prime_counts(bam_paths).items():
        if chromosome not in introns_by_chromosome:
            continue
        for start, end, gene, strand in introns_by_chromosome[chromosome]:
            if start <= position <= end:
                generated[(chromosome, gene, strand)][position] += count

    base = pd.read_csv(polyadb_path, sep="\t")
    required_columns = {"Chrom", "Gene Name", "Strand", "Positions"}
    missing_columns = required_columns - set(base.columns)
    if missing_columns:
        raise ValueError(
            "PolyA_DB file is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    merged = defaultdict(set)
    for _, row in base.iterrows():
        key = (str(row["Chrom"]).strip(), str(row["Gene Name"]).strip().upper(), str(row["Strand"]))
        positions = row["Positions"].strip("[]").split(",") if isinstance(row["Positions"], str) else []
        merged[key].update(int(position) for position in positions if position.strip())

    for (chromosome, gene, strand), position_counts in generated.items():
        merged[(chromosome, gene, strand)].update(
            _cluster_positions(position_counts, merge_distance, min_support)
        )

    output = pd.DataFrame(
        [
            {"Chrom": chrom, "Gene Name": gene, "Strand": strand, "Positions": sorted(positions)}
            for (chrom, gene, strand), positions in sorted(merged.items())
        ],
        columns=["Chrom", "Gene Name", "Strand", "Positions"],
    )
    output.to_csv(output_path, sep="\t", index=False)
    return output_path

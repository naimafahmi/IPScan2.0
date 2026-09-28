"""Annotation parsing and interval construction utilities."""


def make_full_dictionary(annotation_df, chromosomes):
    """Build a chromosome -> gene -> exon/strand lookup."""
    chromosome_dict = {}
    for chromosome in chromosomes:
        gene_dict = {}
        chromosome_rows = annotation_df[annotation_df["chrom"] == chromosome]
        for gene in set(chromosome_rows["name2"]):
            gene_rows = chromosome_rows[chromosome_rows["name2"] == gene]
            exon_list = []
            for _, row in gene_rows.iterrows():
                starts = [value for value in row["exonStarts"].split(",") if value]
                ends = [value for value in row["exonEnds"].split(",") if value]
                for start, end in zip(starts, ends):
                    interval = (int(start), int(end))
                    if interval not in exon_list:
                        exon_list.append(interval)
            gene_dict[gene.strip().upper()] = (exon_list, row["strand"])
        chromosome_dict[chromosome] = gene_dict
    return chromosome_dict


def merge_intervals(intervals):
    """Merge overlapping intervals and return them in ascending order."""
    if not intervals:
        return []
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def merge_chromosome_dictionary(chromosome_dict, chromosomes):
    """Merge exon intervals in every chromosome/gene record."""
    return {
        chromosome: {
            gene: (merge_intervals(exons), strand)
            for gene, (exons, strand) in chromosome_dict[chromosome].items()
        }
        for chromosome in chromosomes
    }


def get_intron_list(exons):
    """Return the genomic gaps between ordered exon intervals."""
    if len(exons) < 2:
        return []
    return [
        (left[1] + 1, right[0] - 1)
        for left, right in zip(exons, exons[1:])
    ]

import ast
import glob, os, csv, sys, time
from pathlib import Path
from operator import attrgetter
import numpy as np
import pandas as pd
import bisect
import warnings
from bisect import bisect_left, bisect_right
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
csv.field_size_limit(sys.maxsize)

from methods.print_and_log import print_and_log
from methods.annotation import (
	get_intron_list as build_intron_list,
	make_full_dictionary as build_annotation_dictionary,
	merge_chromosome_dictionary,
)
from methods.coverage import CoverageIndex, bam_to_coverage_tables, load_coverage_table

warnings.simplefilter(action="ignore", category=FutureWarning)

HUMAN_CHROMOSOMES = tuple(f"chr{index}" for index in range(1, 23)) + ("chrX", "chrY")
MOUSE_CHROMOSOMES = tuple(f"chr{index}" for index in range(1, 20)) + ("chrX", "chrY")
C2_CUTOFF = 10
MAX_C2_C1_RATIO = 0.2

def parse_peak_positions(value):
	positions = ast.literal_eval(str(value))
	if not isinstance(positions, (list, tuple)):
		raise ValueError("Peak Positions values must be lists of genomic positions")
	return sorted(int(position) for position in positions)


def get_chromosomes(species):
	if species in {"hg38", "hg19"}:
		return HUMAN_CHROMOSOMES
	if species in {"mm10", "mm39"}:
		return MOUSE_CHROMOSOMES
	raise ValueError("Unsupported species: choose hg38, hg19, mm10 or mm39")


def prepare_annotation_context(species, annotation_file, logger, peak_file=None):
	chromosomes = get_chromosomes(species)
	ann_df = pd.read_csv(annotation_file, delimiter="\t")
	required_annotation_columns = {"chrom", "name2", "strand", "exonStarts", "exonEnds"}
	missing_annotation_columns = required_annotation_columns - set(ann_df.columns)
	if missing_annotation_columns:
		raise ValueError(
			"Annotation file is missing required columns: "
			+ ", ".join(sorted(missing_annotation_columns))
		)

	print_and_log("Creating the chromosome dictionary...", logger)
	start = time.perf_counter()
	chrom_dict = build_annotation_dictionary(ann_df, chromosomes)
	print_and_log(f"Chromosome dictionary created in {time.perf_counter() - start:.2f} seconds", logger)

	print_and_log("Merging the Exon intervals...", logger)
	start = time.perf_counter()
	merged_chrom_dict = merge_chromosome_dictionary(chrom_dict, chromosomes)
	print_and_log(f"Exon intervals merged in {time.perf_counter() - start:.2f} seconds", logger)

	if peak_file is None:
		peak_file = Path(__file__).resolve().parent.parent / "anno" / f"{species}_polyA_sites.csv"
	peak_df = pd.read_csv(peak_file, delimiter="\t")
	required_peak_columns = {"Chrom", "Gene Name", "Positions"}
	missing_peak_columns = required_peak_columns - set(peak_df.columns)
	if missing_peak_columns:
		raise ValueError(
			"Peak file is missing required columns: "
			+ ", ".join(sorted(missing_peak_columns))
		)
	peak_df["Chrom"] = peak_df["Chrom"].str.strip()
	peak_df["Gene Name"] = peak_df["Gene Name"].str.strip().str.upper()
	peaks_by_gene = {
		(chrom, gene): [parse_peak_positions(value) for value in rows["Positions"]]
		for (chrom, gene), rows in peak_df.groupby(["Chrom", "Gene Name"])
	}

	return {
		"chromosomes": chromosomes,
		"merged_chrom_dict": merged_chrom_dict,
		"peaks_by_gene": peaks_by_gene,
	}


def calculate_event_metrics(strand, exon_intervals, intron_start, intron_end, peak_pos, coverage_index):
	if strand not in {"+", "-"}:
		raise ValueError(f"Unsupported strand value: {strand!r}")

	R1, L1 = 0, 0
	for exon_start, exon_end in exon_intervals:
		read_count, length = coverage_index.sum(exon_start, exon_end)
		R1 += read_count
		L1 += length

	if strand == "+":
		R2, L2 = coverage_index.sum(intron_start, peak_pos)
		R3, L3 = coverage_index.sum(peak_pos + 1, intron_end)
	else:
		R2, L2 = coverage_index.sum(peak_pos, intron_end)
		R3, L3 = coverage_index.sum(intron_start, peak_pos - 1)

	if 0 in (L1, L2, L3) or R1 <= 0 or R3 <= 0:
		return None

	C1 = float(R1) / L1
	C2 = float(R2) / L2
	C3 = float(R3) / L3
	if C1 <= 0:
		return None

	ratio_1 = C2 / C1
	ratio_2 = C3 / C1
	if ratio_2 <= 0:
		return None

	return {
		"ratio_1": ratio_1,
		"ratio_2": ratio_2,
		"R1": R1,
		"R2": R2,
		"R3": R3,
		"L1": L1,
		"L2": L2,
		"L3": L3,
		"C1": C1,
		"C2": C2,
		"C3": C3,
	}


def filter_ipscan_events(events, c2_cutoff=C2_CUTOFF, max_c2_c1_ratio=MAX_C2_C1_RATIO):
	filtered_events = events.copy()
	filtered_events["C2/C1"] = filtered_events["C2"] / filtered_events["C1"]
	return filtered_events.loc[
		(filtered_events["C2"] > c2_cutoff)
		& (filtered_events["C2/C1"] <= max_c2_c1_ratio)
		& (filtered_events["C1"] > filtered_events["C2"])
		& (filtered_events["C1"] > filtered_events["C3"])
		& (filtered_events["C2"] > filtered_events["C3"])
	]


def generate_ipscan_events(species, annotation_file, sample, output_dir, logger, annotation_context=None):
	os.makedirs(output_dir, exist_ok=True)
	bam_path = None
	if os.path.isfile(sample) and sample.lower().endswith(".bam"):
		g1_name = os.path.splitext(os.path.basename(sample))[0]
		bam_path = os.path.abspath(sample)
		sample = os.path.splitext(sample)[0]
	else:
		g1_name = os.path.basename(sample)
	print_and_log("----------------------------------------------------------------", logger)
	print_and_log("| STEP 1: Detecting IPA events for sample: "+g1_name+"...                                  |", logger)
	print_and_log("----------------------------------------------------------------", logger)

	start = time.perf_counter()
	context = annotation_context or prepare_annotation_context(species, annotation_file, logger)
	chromosomes = context["chromosomes"]
	merged_chrom_dict = context["merged_chrom_dict"]
	peaks_by_gene = context["peaks_by_gene"]

	if bam_path is not None:
		if bam_to_coverage_tables(bam_path, sample, chromosomes):
			print_and_log("Generated chromosome coverage tables", logger)
		else:
			print_and_log("Using existing chromosome coverage tables", logger)


	writer_list = []
	output_columns = [
			'Chrom', 
			'Gene Name', 
			'Strand', 
			'Intron Start', 
			'Intron End', 
			'Peak position', 
			'Ratio_1',
			'Ratio_2',
			'Ratio_1/Ratio_2', 
			'R1',
			'R2',
			'R3',
			'L1',
			'L2',
			'L3',
			'C1',
			'C2',
			'C3'
	]

	for chrom in chromosomes:
		#print("Starting:",chrom)

		coverage_file = os.path.join(sample, chrom + ".txt")
		if not os.path.isfile(coverage_file):
			print_and_log(f"Coverage file not found for {chrom}; skipping chromosome", logger)
			continue

		if os.path.getsize(coverage_file) > 0:
			bam_df = load_coverage_table(coverage_file)
			if bam_df is None:
				continue
			position_row = bam_df.iloc[:, 0].tolist()
			bam_list = bam_df.values.tolist()
			coverage_index = CoverageIndex(bam_list, position_row)

			GeneDict_merged = merged_chrom_dict[chrom]
			for geneID in GeneDict_merged.keys():
				peak_positions_rows = peaks_by_gene.get((chrom, geneID), ())

				if peak_positions_rows:
					(mergedExonList, strand) = GeneDict_merged[geneID]
					range_start, range_end = min(mergedExonList)[0], max(mergedExonList)[1]
					if (range_end - range_start) <= 2500000:     ### maximum gene length of human. to stop genes occurring in two regions in the annotation, which makes this range very big
						intronList = build_intron_list(mergedExonList)
						
						for peak_positions in peak_positions_rows:
							for ind in range(len(intronList)):
								(intron_start, intron_end) = intronList[ind]
								peak_start = bisect_left(peak_positions, intron_start)
								peak_end = bisect_right(peak_positions, intron_end + 1)
								for peak_pos in peak_positions[peak_start:peak_end]:
									metrics = calculate_event_metrics(
										strand,
										mergedExonList,
										intron_start,
										intron_end,
										peak_pos,
										coverage_index,
									)
									if metrics is not None:
												writer_list.append(
													(
														chrom, 
														geneID, 
														strand, 
														intron_start, 
														intron_end, 
														peak_pos, 
													metrics["ratio_1"],
													metrics["ratio_2"],
													round(metrics["ratio_1"] / metrics["ratio_2"], 3),
													metrics["R1"],
													metrics["R2"],
													metrics["R3"],
													metrics["L1"],
													metrics["L2"],
													metrics["L3"],
													metrics["C1"],
													metrics["C2"],
													metrics["C3"]
													)
												)
												#print("Found", chrom, geneID, strand, intron_start, intron_end, peak_pos, ratio_1, ratio_2, round(ratio_1/ratio_2, 3), R1, R2, R3, L1, L2, L3, C1, C2, C3)


	df_out = pd.DataFrame(writer_list, columns=output_columns)
	full_output_path = os.path.join(output_dir, "IPScan_"+g1_name+"_full.csv")
	df_out.to_csv(full_output_path, sep="\t", index=False)

	newdf = filter_ipscan_events(df_out)
	filtered_output_path = os.path.join(output_dir, "IPScan_"+g1_name+".csv")
	newdf.to_csv(filtered_output_path, sep='\t', header=True, index=False)
	print_and_log(f"IPA sites detected in {time.perf_counter() - start:.2f} seconds", logger)
	return {"full": full_output_path, "filtered": filtered_output_path}
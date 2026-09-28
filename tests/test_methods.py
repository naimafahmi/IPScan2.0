import pytest
import pandas as pd

from methods.annotation import get_intron_list, merge_intervals
from methods.methods import calculate_event_metrics, filter_ipscan_events
from methods.coverage import (
    CoverageIndex,
    bam_to_coverage_tables,
    count_total_read_count,
    load_coverage_table,
)


def test_merge_intervals_merges_overlaps_and_is_ordered():
    assert merge_intervals([(10, 20), (1, 5), (4, 12)]) == [(1, 20)]


def test_get_intron_list_returns_gaps_between_exons():
    assert get_intron_list([(1, 10), (21, 30), (41, 50)]) == [(11, 20), (31, 40)]


def test_count_total_read_count_handles_inclusive_bounds():
    rows = [[1, 2], [2, 3], [3, 5], [5, 7]]
    assert count_total_read_count(2, 3, rows, [1, 2, 3, 5]) == (8, 2)


def test_coverage_index_matches_inclusive_interval_semantics():
    rows = [[1, 2], [2, 3], [3, 5], [5, 7]]
    index = CoverageIndex(rows, [1, 2, 3, 5])
    assert index.sum(2, 3) == (8, 2)


def test_bam_conversion_reuses_complete_coverage_tables(tmp_path):
    output_dir = tmp_path / "sample"
    output_dir.mkdir()
    chromosomes = ("chr1", "chr2")
    for chromosome in chromosomes:
        (output_dir / f"{chromosome}.txt").write_text("1\t1\n")

    assert bam_to_coverage_tables("missing.bam", output_dir, chromosomes) is False


def test_load_coverage_table_validates_two_column_sorted_input(tmp_path):
    coverage_path = tmp_path / "chr1.txt"
    coverage_path.write_text("2\t3\n5\t7\n")

    coverage = load_coverage_table(coverage_path)

    assert coverage.values.tolist() == [[2, 3], [5, 7]]


def test_load_coverage_table_rejects_unsorted_positions(tmp_path):
    coverage_path = tmp_path / "chr1.txt"
    coverage_path.write_text("5\t7\n2\t3\n")

    with pytest.raises(ValueError, match="sorted and unique"):
        load_coverage_table(coverage_path)


def test_calculate_event_metrics_handles_both_strands():
    rows = [[1, 2], [2, 2], [3, 1], [4, 1], [5, 2], [6, 2]]
    index = CoverageIndex(rows, [row[0] for row in rows])

    plus_metrics = calculate_event_metrics("+", [(1, 2), (5, 6)], 3, 4, 3, index)
    minus_metrics = calculate_event_metrics("-", [(1, 2), (5, 6)], 3, 4, 4, index)

    assert plus_metrics["ratio_1"] == 0.5
    assert minus_metrics["ratio_2"] == 0.5


def test_calculate_event_metrics_rejects_invalid_strand():
    index = CoverageIndex([[1, 1]], [1])

    with pytest.raises(ValueError, match="Unsupported strand"):
        calculate_event_metrics("?", [(1, 1)], 2, 3, 2, index)


def test_filter_ipscan_events_applies_named_thresholds():
    events = pd.DataFrame(
        [
            {"C1": 100, "C2": 10, "C3": 1},
            {"C1": 100, "C2": 20, "C3": 2},
            {"C1": 50, "C2": 20, "C3": 2},
        ]
    )

    filtered = filter_ipscan_events(events)

    assert len(filtered) == 1
    assert filtered.iloc[0]["C2"] == 20
    assert filtered.iloc[0]["C2/C1"] == 0.2
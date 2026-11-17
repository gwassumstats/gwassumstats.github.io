"""h² on the liability scale: build/liability.py."""
import math

import pytest

from build.liability import factor, prevalences, scaled
from conftest import write


def test_h2_is_converted_with_the_prevalence_and_the_sample_s_share_of_cases():
    assert factor(0.01, 500, 500) == pytest.approx(0.5519, abs=1e-4)  # K² (1 - K)² / (φ(t)² v (1 - v)), Lee et al. 2011: K 1%, half cases
    assert factor(0.01, 500, 500, n=1000) == factor(0.01, 500, 500, n=600) == factor(0.01, 500, 500)  # a file of a part of the sample has its share of cases
    assert factor(0.01, 500, 500, n=1150) == pytest.approx(1.15 * factor(0.01, 500, 500))  # an N that counts some individuals twice lowered h² as many times
    assert factor(0.2, 100, 900) < factor(0.2, 500, 500) * 4  # the same h² in a sample of fewer cases is more on the liability scale, by 0.25 / 0.09 here


def test_a_study_has_none_without_cases_a_prevalence_or_with_proxy_cases():
    assert scaled("binary", False, 500, 500, 1000, 0.1, 0.01, 0.01) == pytest.approx((0.05519, 0.005519), abs=1e-5)
    assert scaled("binary", True, 500, 500, 1000, 0.1, 0.01, 0.01) == (None, None)  # proxy cases: relatives of a case
    assert scaled("binary", False, math.nan, math.nan, 1000, 0.1, 0.01, 0.01) == scaled("binary", False, 500, 500, 1000, 0.1, 0.01, None) == scaled("quantitative", False, math.nan, math.nan, 1000, 0.1, 0.01, 0.01) == (None, None)


def test_one_prevalence_a_trait_and_ancestry(tmp_path):
    head = ["trait", "anc", "prevalence", "source", "quote"]
    write(tmp_path / "prevalence.tsv", head, [["# a comment", "", "", "", ""], ["T", "EUR", "0.01", "PMID 1", "1%"], ["T", "EAS", "0.02", "PMID 2", "2%"]])
    assert prevalences(tmp_path) == {("T", "EUR"): {"k": 0.01, "source": "PMID 1", "quote": "1%"}, ("T", "EAS"): {"k": 0.02, "source": "PMID 2", "quote": "2%"}}
    write(tmp_path / "prevalence.tsv", head, [["T", "EUR", "0.01", "PMID 1", "1%"], ["T", "EUR", "0.02", "PMID 2", "2%"]])
    with pytest.raises(ValueError, match="T \\(EUR\\) has two rows"):
        prevalences(tmp_path)
    write(tmp_path / "prevalence.tsv", head, [["T", "EUR", "1.5", "PMID 1", "150%"]])
    with pytest.raises(ValueError, match="not between 0 and 1"):
        prevalences(tmp_path)

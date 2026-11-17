"""The effective N of a study, and the choice of a study for each trait and ancestry: build/recommend.py."""
import pandas as pd
import pytest

from build.recommend import chosen, effective, older, sex


def test_effective_n_is_that_of_the_cases_and_controls():
    assert effective("quantitative", 1000) == 1000
    assert effective("binary", 300, 100, 200) == pytest.approx(4 / (1 / 100 + 1 / 200))  # 4 v (1 - v) N is the same number
    assert effective("binary", 300) is None and effective("binary", 300, float("nan"), float("nan")) is None and effective("", 300, 100, 200) is None  # cases, or the kind, not known
    # the cases and controls are those of the file's sample, whose N is the largest of its SNPs: the inventory's is their median
    assert effective("binary", 280, 100, 200, largest=300) == pytest.approx(4 / (1 / 100 + 1 / 200))
    # a larger sample is recorded than the file holds: the file has its share of cases
    assert effective("binary", 150, 100, 200, largest=150) == pytest.approx(4 * (1 / 3) * (2 / 3) * 150)
    # the study's N is above its sample, as where it counts another stage: the sample's cases and controls stand
    assert effective("binary", 900, 100, 200, largest=900) == pytest.approx(4 / (1 / 100 + 1 / 200))


def test_a_short_name_says_the_sex_of_a_study_and_a_date_how_far_it_is_known():
    assert [sex(name) for name in ("Lee_2020_BMIFemale", "Lee_2020_BMIMale", "Lee_2020_BMI", "Male_2020_BMI", "Lee_2020_FemaleInfertility")] == ["Female", "Male", "", "", "Female"]
    assert older("2018-07", "2018-10-15") and older("2018-07-01", "2018-07-02") and not older("2018-07", "2018-07-15") and not older("", "2018") and not older("2019", "2018-07")


def test_the_largest_effective_n_among_the_studies_that_pass_is_recommended():
    d = pd.DataFrame([
        ("big_warns", "T", "EUR", "warn", 950, 900.0, 30.0, "A_2022_T", "2022-01-01"),
        ("large", "T", "EUR", "pass", 600, 500.0, 10.0, "B_2020_T", "2020-06-01"),
        ("many", "T", "EUR", "pass", 5000, 400.0, 40.0, "M_2021_T", "2021-01-01"),  # many controls, few cases: more individuals, less power
        ("old", "T", "EUR", "pass", 100, 100.0, 20.0, "C_2010_T", "2010-01-01"),
        ("old_fails", "T", "EUR", "fail", 60, 50.0, 2.0, "D_2012_T", "2012-01-01"),
        ("women", "T", "EUR", "pass", 300, 300.0, 9.0, "E_2015_TFemale", "2015-01-01"),
        ("newer", "T", "EUR", "pass", 200, 200.0, 9.0, "F_2024_T", "2024-01-01"),
        ("same_paper", "T", "EUR", "pass", 450, 450.0, 9.0, "B_2020_TSubtype", "2020-06-01"),
        ("unknown_warns", "T", "EUR", "warn", 9000, None, 99.0, "G_2019_T", "2019-01-01"),  # it does not pass: its unknown cases change nothing
        ("east", "T", "EAS", "pass", 80, 80.0, 5.0, "H_2018_TEAS", "2018-01-01"),
        ("tie_a", "U", "EUR", "pass", 700, 700.0, 5.0, "I_2018_U", "2018-01"),
        ("tie_b", "U", "EUR", "pass", 700, 700.0, 6.0, "J_2018_U", "2018-01-20"),
        ("alone", "V", "EUR", "pass", 40, None, 3.0, "K_2018_V", ""),
        ("none", "W", "EUR", "warn", 10, 10.0, 1.0, "L_2018_W", "2018-01-01"),
        ("known", "X", "EUR", "pass", 900, 300.0, 8.0, "N_2021_X", "2021-01-01"),
        ("not_public", "X", "EUR", "pass", 9000, None, 4.0, "O_2018_X", "2018-01-01"),  # a biobank's study of few cases: the most individuals, and not the most power
        ("strong", "X", "EUR", "pass", 500, 480.0, 9.0, "P_2022_X", "2022-01-01"),
        ("equal_small", "Y", "EUR", "pass", 100, None, 7.0, "Q_2018_Y", "2018-01-01"),
        ("equal_large", "Y", "EUR", "pass", 200, 150.0, 7.0, "R_2019_Y", "2019-01-01"),
    ], columns=["id", "of", "anc", "qc", "n", "neff", "h2_z", "name", "published"]).set_index("id")
    recommended, superseded = chosen(d)
    # one a trait and ancestry, by effective N; h² z settles a tie; a trait of which nothing passes has none
    # where the cases of a study that passes are not known, h² z ranks them all, whatever their N, and N settles a tie
    assert recommended.dropna().to_dict() == {"large": "neff", "east": "neff", "tie_b": "neff", "alone": "z", "strong": "z", "equal_large": "z"}
    # older studies of the same sex are superseded, whatever the checks say of them; a study of one sex, a newer one and one of the same paper are other analyses
    assert superseded.dropna().to_dict() == {"old": "large", "old_fails": "large", "unknown_warns": "large", "known": "strong", "not_public": "strong", "equal_small": "equal_large"}

from __future__ import annotations

import hashlib
import json
import math
import os
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import chi2, norm
import statsmodels
import statsmodels.api as sm
import statsmodels.formula.api as smf

REPO_ROOT = Path(__file__).resolve().parents[2]


def required_path(name: str) -> Path:
    value = os.environ.get(name)
    if not value:
        raise RuntimeError(
            f"Required environment variable {name} is not set. "
            "See config/paths.example.env and docs/DATA_ACCESS.md."
        )
    return Path(value).expanduser().resolve()


RAW = required_path("CHARLS_DATA_ROOT")
MAP = required_path("CHARLS_ITEM_MAP_FILE")
OUTPUT_ROOT = Path(os.environ.get("CHARLS_PRIMARY_OUTPUT", REPO_ROOT / "outputs" / "charls_v2")).expanduser().resolve()
OUT_PRIMARY = OUTPUT_ROOT / "primary"
OUT_AUDIT = OUTPUT_ROOT / "audit"
OUT_SENS = OUTPUT_ROOT / "sensitivity"
OUT_QC = OUTPUT_ROOT / "qc"

FILES = {
    "grip2011": RAW / "2011" / "biomarkers.dta",
    "grip2015": RAW / "2015" / "Biomarker.dta",
    "hsf2011": RAW / "2011" / "health_status_and_functioning.dta",
    "hsf2015": RAW / "2015" / "Health_Status_and_Functioning.dta",
    "hsf2018": RAW / "2018" / "Health_Status_and_Functioning.dta",
    "hsf2020": RAW / "2020" / "Health_Status_and_Functioning.dta",
    "demo2011": RAW / "2011" / "demographic_background.dta",
    "demo2015": RAW / "2015" / "Demographic_Background.dta",
}


def write_csv(path: Path, rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False, encoding="utf-8-sig")


def native(s: pd.Series) -> pd.Series:
    def f(v):
        if pd.isna(v):
            return pd.NA
        if isinstance(v, bytes):
            v = v.decode("utf-8")
        v = str(v).strip()
        return pd.NA if v == "" else v
    return s.map(f).astype("string")


def xid(s: pd.Series, wave: int) -> pd.Series:
    x = native(s)
    z = x.dropna()
    if wave == 2011:
        if not z.str.fullmatch(r"\d{11}").all():
            raise ValueError("Unexpected 2011 ID format")
        x.loc[z.index] = z.str[:9] + "0" + z.str[-2:]
    elif not z.str.fullmatch(r"\d{12}").all():
        raise ValueError(f"Unexpected {wave} ID format")
    return x


def read(path: Path, columns=None) -> pd.DataFrame:
    return pd.read_stata(path, columns=columns, convert_categoricals=False, preserve_dtypes=True)


def num(x):
    return pd.to_numeric(x, errors="coerce")


def merge1(a, b):
    return a.merge(b, on="pid", how="left", validate="one_to_one")


def grip(wave: int) -> pd.DataFrame:
    d = read(FILES[f"grip{wave}"], ["ID", "qc003", "qc004", "qc005", "qc006"])
    d["pid"] = xid(d.ID, wave)
    vals = d[["qc003", "qc004", "qc005", "qc006"]].apply(num)
    vals = vals.where((vals > 0) & (vals < 100))
    d[f"grip_{wave}"] = vals.max(axis=1)
    return d[["pid", f"grip_{wave}"]].drop_duplicates("pid")


amap = pd.read_csv(MAP, dtype=str).fillna("")


def outcomes(wave: int) -> pd.DataFrame:
    mp = amap[amap.wave.eq(str(wave))]
    extra = ["db001", "db004", "db005", "db006", "db007", "db008", "db009", "db015", "db016", "db034"] if wave == 2018 else []
    variables = list(dict.fromkeys(list(mp.raw_variable_name) + extra))
    d = read(FILES[f"hsf{wave}"], ["ID"] + variables)
    d["pid"] = xid(d.ID, wave)
    adl, iadl = [], []
    for (domain, activity), z in mp.groupby(["domain", "activity"]):
        variable = z.raw_variable_name.iloc[0]
        value = num(d[variable])
        valid = value.isin([1, 2, 3, 4])
        dependency_codes = z.loc[z.need_help.eq("True") | z.unable_to_do.eq("True"), "raw_code"].astype(float)
        state = value.isin(dependency_codes).where(valid)
        (adl if domain == "ADL" else iadl).append(state)

    def combine(parts):
        q = pd.concat(parts, axis=1)
        return q.any(axis=1).where(q.notna().all(axis=1))

    adl_state, iadl_state = combine(adl), combine(iadl)
    structural = pd.Series(False, index=d.index)
    if wave == 2018:
        structural = (
            d[["db001", "db004", "db005", "db006", "db007", "db008", "db009"]].apply(num).eq(1).all(axis=1)
            & d[["db010", "db011", "db012", "db013", "db014"]].isna().all(axis=1)
            & d["db015"].isna()
            & num(d["db016"]).isin([1, 2, 3, 4])
        )
        adl_state = adl_state.mask(structural, False)
    return pd.DataFrame({
        "pid": d.pid,
        f"strict_{wave}": combine([adl_state, iadl_state]),
        f"adl_{wave}": adl_state,
        f"iadl_{wave}": iadl_state,
        f"adl_structural_skip_{wave}": structural,
    }).drop_duplicates("pid")


def covariates() -> pd.DataFrame:
    d11 = read(FILES["demo2011"], ["ID", "rgender", "bd001", "be001"])
    d11["pid"] = xid(d11.ID, 2011)
    d15 = read(FILES["demo2015"], ["ID", "ba004_w3_1", "bb001_w3_2"])
    d15["pid"] = xid(d15.ID, 2015)
    h11cols = ["ID", "da059", "da061", "da067", "da069"] + [f"da007_{i}_" for i in range(1, 15)]
    h11 = read(FILES["hsf2011"], h11cols)
    h11["pid"] = xid(h11.ID, 2011)
    cov = d11[["pid", "rgender", "bd001", "be001"]].merge(
        d15[["pid", "ba004_w3_1", "bb001_w3_2"]], on="pid", how="left", validate="one_to_one"
    ).merge(h11.drop(columns="ID"), on="pid", how="left", validate="one_to_one")
    cov["age_2015"] = 2015 - num(cov.ba004_w3_1)
    cov["sex"] = num(cov.rgender).map({1: "male", 2: "female"})
    ed = num(cov.bd001)
    cov["education"] = pd.cut(ed, [0, 3, 4, 5, 7, 99], labels=["no_formal", "primary", "middle", "high_vocational", "college_plus"])
    cov["residence"] = np.where(num(cov.bb001_w3_2).isin([1, 2, 3]), "urban", np.where(num(cov.bb001_w3_2).isin([4, 5, 6, 7]), "rural", pd.NA))
    cov["marital_status"] = np.where(num(cov.be001).isin([1, 2]), "married", np.where(num(cov.be001).isin([3, 4, 5, 6]), "not_married", pd.NA))
    cov["smoking"] = np.select(
        [num(cov.da059).eq(2), num(cov.da059).eq(1) & num(cov.da061).eq(2), num(cov.da059).eq(1) & num(cov.da061).eq(1)],
        ["never", "former", "current"], default=None,
    )
    cov["alcohol"] = np.select(
        [num(cov.da067).isin([1, 2]), num(cov.da067).eq(3) & num(cov.da069).isin([2, 3]), num(cov.da067).eq(3) & num(cov.da069).eq(1)],
        ["current", "former", "never"], default=None,
    )
    chronic_items = pd.concat([num(cov[f"da007_{i}_"]).map({1: 1, 2: 0}) for i in range(1, 15)], axis=1)
    cov["chronic_condition_count"] = chronic_items.sum(axis=1).where(chronic_items.notna().all(axis=1))
    return cov[["pid", "age_2015", "sex", "education", "residence", "marital_status", "smoking", "alcohol", "chronic_condition_count"]]


def fit(formula: str, data: pd.DataFrame):
    model = smf.glm(formula, data=data, family=sm.families.Binomial(link=sm.families.links.CLogLog()))
    ordinary = model.fit()
    robust = model.fit(cov_type="cluster", cov_kwds={"groups": data["study_id"]})
    return ordinary, robust


def effect(model, term: str, scale: float = 1.0):
    b = float(model.params[term]) * scale
    se = float(model.bse[term]) * abs(scale)
    return {
        "term": term,
        "log_HR": b,
        "HR": math.exp(b),
        "CI_low": math.exp(b - 1.96 * se),
        "CI_high": math.exp(b + 1.96 * se),
        "p_value": 2 * norm.sf(abs(b / se)),
    }


def descriptive(group: pd.DataFrame, label: str):
    out = []
    for v in ["age_2015", "grip_2011", "grip_decline", "chronic_condition_count"]:
        x = num(group[v]).dropna()
        out.append({"group": label, "variable": v, "available_n": len(x), "mean": x.mean(), "sd": x.std(), "median": x.median()})
    for v in ["sex", "education"]:
        n = group[v].notna().sum()
        for level, count in group[v].value_counts(dropna=True).items():
            out.append({"group": label, "variable": f"{v}:{level}", "available_n": n, "mean": count / n if n else np.nan, "sd": np.nan, "median": np.nan})
    return out


def main():
    for p in FILES.values():
        if not p.exists():
            raise FileNotFoundError(p)
    if not MAP.exists():
        raise FileNotFoundError(MAP)

    g11, g15 = grip(2011), grip(2015)
    cohort = g11.merge(g15, on="pid", how="outer", validate="one_to_one")
    cohort["grip_decline"] = cohort.grip_2011 - cohort.grip_2015
    for wave in [2015, 2018, 2020]:
        cohort = merge1(cohort, outcomes(wave))
    cohort = merge1(cohort, covariates())

    valid11 = cohort.grip_2011.notna()
    valid15 = valid11 & cohort.grip_2015.notna()
    change = valid15 & cohort.grip_decline.notna()
    risk = change & cohort.strict_2015.eq(False)
    any_followup = cohort.strict_2018.notna() | cohort.strict_2020.notna()
    ascertained = risk & any_followup

    pp = []
    for _, r in cohort[ascertained].iterrows():
        if pd.notna(r.strict_2018):
            q = r.to_dict(); q.update(interval="2015-2018", event=int(bool(r.strict_2018))); pp.append(q)
            if not bool(r.strict_2018) and pd.notna(r.strict_2020):
                q = r.to_dict(); q.update(interval="2018-2020", event=int(bool(r.strict_2020))); pp.append(q)
    pp = pd.DataFrame(pp)
    pp["study_id"] = pp.pid.map(lambda x: hashlib.sha256(("v2-charls:" + str(x)).encode()).hexdigest()[:16])

    modelvars = ["event", "study_id", "interval", "grip_decline", "grip_2011", "age_2015", "sex", "education", "residence", "marital_status", "smoking", "alcohol", "chronic_condition_count"]
    analytic = pp.dropna(subset=modelvars).copy()
    included_ids = set(analytic.study_id.unique())
    participant = pp.sort_values("interval").drop_duplicates("study_id").copy()
    participant["included"] = participant.study_id.isin(included_ids)

    complete_person_mask = cohort.pid.map(lambda x: hashlib.sha256(("v2-charls:" + str(x)).encode()).hexdigest()[:16]).isin(included_ids)
    flow = [
        {"step": "valid grip 2011", "N": int(valid11.sum())},
        {"step": "valid grip 2015 among valid 2011", "N": int(valid15.sum())},
        {"step": "calculable grip change", "N": int(change.sum())},
        {"step": "2015 dependency-free risk set", "N": int(risk.sum())},
        {"step": "at least one valid person-period", "N": int(pp.study_id.nunique())},
        {"step": "primary covariates complete", "N": int(analytic.study_id.nunique())},
    ]
    write_csv(OUT_PRIMARY / "CHARLS_V2_SELECTION_FLOW.csv", flow)
    (OUT_PRIMARY / "CHARLS_V2_SELECTION_FLOW.md").write_text(
        "# CHARLS V2 selection flow\n\n" + "\n".join(f"- {r['step']}: N={r['N']:,}" for r in flow)
        + f"\n\nFinal: {analytic.study_id.nunique():,} participants, {len(analytic):,} person-periods, {int(analytic.event.sum()):,} events. Old biomarker-complete primary: 2,499 participants, 4,449 person-periods, 674 events.\n",
        encoding="utf-8",
    )

    # Missingness reasons among the ascertained risk-set person records.
    missing_rows = []
    for v in modelvars[3:]:
        missing_rows.append({"variable": v, "missing_participants": int(participant[v].isna().sum()), "denominator": len(participant), "percent": 100 * participant[v].isna().mean()})
    missing_rows += [
        {"variable": "outcome_not_ascertainable_in_risk_set", "missing_participants": int((risk & ~any_followup).sum()), "denominator": int(risk.sum()), "percent": 100 * int((risk & ~any_followup).sum()) / int(risk.sum())},
        {"variable": "no_valid_person_period_among_outcome_ascertained", "missing_participants": int(ascertained.sum() - pp.study_id.nunique()), "denominator": int(ascertained.sum()), "percent": 100 * (int(ascertained.sum() - pp.study_id.nunique())) / int(ascertained.sum())},
    ]
    write_csv(OUT_AUDIT / "CHARLS_MISSINGNESS_AUDIT.csv", missing_rows)
    write_csv(OUT_AUDIT / "CHARLS_INCLUDED_EXCLUDED_DESCRIPTIVE.csv", descriptive(participant[participant.included], "included") + descriptive(participant[~participant.included], "excluded"))

    sd_change = analytic.drop_duplicates("study_id").grip_decline.std()
    analytic["grip_decline_sd"] = analytic.grip_decline / sd_change
    formula_b = "event ~ C(interval) + grip_decline_sd + grip_2011 + age_2015 + C(sex) + C(education) + C(residence) + C(marital_status) + C(smoking) + C(alcohol) + chronic_condition_count"
    formula_a = formula_b.replace(" + grip_decline_sd", "")
    ordinary_a, robust_a = fit(formula_a, analytic)
    ordinary_b, robust_b = fit(formula_b, analytic)

    res_sd = effect(robust_b, "grip_decline_sd")
    res_5kg = effect(robust_b, "grip_decline_sd", 5.0 / sd_change)
    primary_rows = [
        {"analysis": "per 1 SD greater grip decline", "scale_value": sd_change, **res_sd},
        {"analysis": "per 5 kg greater grip decline", "scale_value": 5.0, **res_5kg},
    ]
    for r in primary_rows:
        r.update(participants=analytic.study_id.nunique(), person_periods=len(analytic), events=int(analytic.event.sum()), covariance="participant-clustered robust")
    write_csv(OUT_PRIMARY / "CHARLS_PRIMARY_GRIP_MODEL.csv", primary_rows)
    (OUT_PRIMARY / "CHARLS_PRIMARY_GRIP_MODEL.md").write_text(
        f"# CHARLS V2 primary grip model\n\nN={analytic.study_id.nunique():,}; person-periods={len(analytic):,}; events={int(analytic.event.sum()):,}; SD of change={sd_change:.6f} kg.\n\n"
        f"- Per SD: HR {res_sd['HR']:.6f} (95% CI {res_sd['CI_low']:.6f}-{res_sd['CI_high']:.6f}), p={res_sd['p_value']:.8g}.\n"
        f"- Per 5 kg: HR {res_5kg['HR']:.6f} (95% CI {res_5kg['CI_low']:.6f}-{res_5kg['CI_high']:.6f}), p={res_5kg['p_value']:.8g}.\n",
        encoding="utf-8",
    )

    lr = 2 * (ordinary_b.llf - ordinary_a.llf)
    comparison = [
        {"model": "baseline_only", "N": analytic.study_id.nunique(), "person_periods": len(analytic), "events": int(analytic.event.sum()), "log_likelihood": ordinary_a.llf, "AIC": ordinary_a.aic, "BIC": ordinary_a.bic_llf, "LR_statistic_vs_baseline": np.nan, "LR_df": np.nan, "LR_p": np.nan, "discrimination": "NOT_IMPLEMENTED_DUE_TO_MODEL_STRUCTURE", "Brier": "NOT_IMPLEMENTED_DUE_TO_MODEL_STRUCTURE"},
        {"model": "baseline_plus_change", "N": analytic.study_id.nunique(), "person_periods": len(analytic), "events": int(analytic.event.sum()), "log_likelihood": ordinary_b.llf, "AIC": ordinary_b.aic, "BIC": ordinary_b.bic_llf, "LR_statistic_vs_baseline": lr, "LR_df": 1, "LR_p": chi2.sf(lr, 1), "discrimination": "NOT_IMPLEMENTED_DUE_TO_MODEL_STRUCTURE", "Brier": "NOT_IMPLEMENTED_DUE_TO_MODEL_STRUCTURE"},
    ]
    write_csv(OUT_PRIMARY / "BASELINE_VS_CHANGE_MODEL_COMPARISON.csv", comparison)
    (OUT_PRIMARY / "BASELINE_VS_CHANGE_MODEL_COMPARISON.md").write_text(
        f"# Baseline-only versus baseline plus change\n\nSame analytic sample in both models. LR={lr:.6f}, df=1, p={chi2.sf(lr,1):.8g}. AIC: {ordinary_a.aic:.3f} versus {ordinary_b.aic:.3f}; BIC: {ordinary_a.bic_llf:.3f} versus {ordinary_b.bic_llf:.3f}. These statistics assess additional model information, not improved clinical prediction.\n",
        encoding="utf-8",
    )

    # Equivalent baseline/follow-up parameterization.
    analytic["grip_2015"] = analytic.grip_2011 - analytic.grip_decline
    formula_c = formula_a.replace("grip_2011", "grip_2011 + grip_2015")
    ordinary_c, robust_c = fit(formula_c, analytic)
    fitted_diff = float(np.max(np.abs(ordinary_b.fittedvalues - ordinary_c.fittedvalues)))
    (OUT_SENS / "CHANGE_SCORE_PARAMETERIZATION_AUDIT.md").write_text(
        "# Change-score parameterization audit\n\n"
        f"Model B log likelihood: {ordinary_b.llf:.12f}. Model C log likelihood: {ordinary_c.llf:.12f}. Maximum absolute fitted-probability difference: {fitted_diff:.3e}. "
        "Because decline = baseline - follow-up, the two formulations span the same design space (up to numerical precision). The decline coefficient conditional on baseline is the negative of the follow-up coefficient in the baseline-plus-follow-up formulation. It is a conditional association and not a pure causal effect of dynamic change; baseline-change mathematical coupling and measurement error remain relevant.\n",
        encoding="utf-8",
    )

    # Interval interaction.
    formula_i = formula_b.replace("C(interval) + grip_decline_sd", "C(interval) * grip_decline_sd")
    ordinary_i, robust_i = fit(formula_i, analytic)
    iterms = [x for x in robust_i.params.index if ":grip_decline_sd" in x]
    if len(iterms) != 1:
        raise RuntimeError(f"Unexpected interval interaction terms: {iterms}")
    interval_effect = effect(robust_i, iterms[0])
    interval_rows = [{"component": "interaction", **interval_effect}]
    interval_rows += [{"component": "events", "interval": k, "events": int(v.event.sum()), "person_periods": len(v)} for k, v in analytic.groupby("interval")]
    write_csv(OUT_SENS / "INTERVAL_INTERACTION.csv", interval_rows)

    # Extended adjustment deliberately not improvised without a verified frozen construction.
    write_csv(OUT_SENS / "EXTENDED_ADJUSTMENT_RESULTS.csv", [{"model": "BMI_CESD_extended", "status": "NOT_RUN_RELIABLE_FROZEN_BMI_CESD_CONSTRUCTION_NOT_LOCATED_IN_RELEASE_SCRIPT", "N": np.nan, "person_periods": np.nan, "events": np.nan, "grip_HR": np.nan, "lower_CI": np.nan, "upper_CI": np.nan, "p": np.nan, "covariate_set": "primary + BMI + CES-D10"}])

    # Pre-outcome IPCW specification and feasibility status.
    (OUT_SENS / "IPCW_SPECIFICATION.md").write_text(
        "# IPCW specification (frozen before any weighted outcome model)\n\nTarget: participants in the 2015 dependency-free grip-change risk set. Denominator observation model candidates: age, sex, education, marital status, residence, smoking, alcohol, chronic conditions, 2011 grip, grip change, and prior observed functional state. Numerator: intercept and interval. Stabilized weights; diagnostics include min, max, mean, median, p1, p99, effective sample size; truncation sensitivity at p1/p99.\n\nStatus: `IPCW_NOT_IDENTIFIABLE_FROM_AVAILABLE_DATA`. The released construction does not retain a validated wave-specific observation/attrition state that distinguishes death, non-response, structural missingness and item non-response for every eligible interval. No weighted outcome model was fitted and no substitute definition was created.\n",
        encoding="utf-8",
    )
    write_csv(OUT_SENS / "IPCW_RESULTS.csv", [{"status": "IPCW_NOT_IDENTIFIABLE_FROM_AVAILABLE_DATA", "reason": "validated wave-specific observation/attrition states unavailable", "weighted_HR": np.nan, "CI_low": np.nan, "CI_high": np.nan}])

    # QC gates.
    key_dupes = int(analytic.duplicated(["study_id", "interval"]).sum())
    multi_events = int((analytic.groupby("study_id").event.sum() > 1).sum())
    post_event = 0
    for _, z in analytic.sort_values(["study_id", "interval"]).groupby("study_id"):
        ev = np.where(z.event.to_numpy() == 1)[0]
        if len(ev) and ev[0] < len(z) - 1:
            post_event += 1
    bridge = int(((cohort.strict_2018.isna()) & cohort.strict_2020.notna() & cohort.pid.isin(pp.pid)).sum())
    X = ordinary_b.model.exog
    rank = np.linalg.matrix_rank(X)
    qc = [
        ("participant IDs unique at landmark", int(cohort.pid.duplicated().sum()) == 0, int(cohort.pid.duplicated().sum())),
        ("person-period duplicates", key_dupes == 0, key_dupes),
        ("event only once per participant", multi_events == 0, multi_events),
        ("stop after first event", post_event == 0, post_event),
        ("no invalid bridge interval", bridge == 0, bridge),
        ("baseline dependency excluded", not analytic.strict_2015.eq(True).any(), int(analytic.strict_2015.eq(True).sum())),
        ("exposure precedes outcome", True, "2011-2015 change; 2018/2020 outcome"),
        ("design matrix full rank", rank == X.shape[1], f"{rank}/{X.shape[1]}"),
        ("finite coefficients", np.isfinite(robust_b.params).all(), bool(np.isfinite(robust_b.params).all())),
        ("finite standard errors", np.isfinite(robust_b.bse).all(), bool(np.isfinite(robust_b.bse).all())),
        ("cluster robust covariance available", robust_b.cov_type == "cluster", robust_b.cov_type),
        ("event counts reconcile", int(analytic.event.sum()) == int(analytic.groupby("study_id").event.max().sum()), int(analytic.event.sum())),
        ("displayed numbers machine-readable", True, "CHARLS_PRIMARY_GRIP_MODEL.csv"),
    ]
    qc_rows = [{"check": a, "status": "PASS" if b else "FAIL", "detail": c} for a, b, c in qc]
    write_csv(OUT_QC / "CHARLS_PRIMARY_QC.csv", qc_rows)
    if any(r["status"] == "FAIL" for r in qc_rows):
        raise RuntimeError("Core CHARLS QC failure")

    summary = {
        "status": "CHARLS_V2_PRIMARY_COMPLETE",
        "participants": int(analytic.study_id.nunique()),
        "person_periods": int(len(analytic)),
        "events": int(analytic.event.sum()),
        "sd_change_kg": float(sd_change),
        "HR_per_SD": res_sd,
        "HR_per_5kg": res_5kg,
        "LR": float(lr),
        "LR_p": float(chi2.sf(lr, 1)),
        "interval_interaction": interval_effect,
        "software": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__, "statsmodels": statsmodels.__version__},
    }
    (OUT_PRIMARY / "CHARLS_V2_SUMMARY.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()

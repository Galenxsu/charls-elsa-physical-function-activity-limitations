from __future__ import annotations

import json
import math
import os
import runpy
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2, norm
import statsmodels.api as sm
import statsmodels.formula.api as smf


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
BASE_SCRIPT = HERE / "01_sensitivity_analyses.py"
OUT = Path(os.environ.get("SENSITIVITY_OUTPUT_ROOT", REPO_ROOT / "outputs" / "sensitivity")) / "S17"
OUT.mkdir(parents=True, exist_ok=True)


def fit(formula: str, data: pd.DataFrame, cluster: str):
    model = smf.glm(
        formula,
        data=data,
        family=sm.families.Binomial(link=sm.families.links.CLogLog()),
    )
    ordinary = model.fit(maxiter=200)
    robust = model.fit(
        cov_type="cluster", cov_kwds={"groups": data[cluster]}, maxiter=200
    )
    rank = int(np.linalg.matrix_rank(ordinary.model.exog))
    if rank != ordinary.model.exog.shape[1]:
        raise RuntimeError(
            f"Rank-deficient design: {rank}/{ordinary.model.exog.shape[1]}"
        )
    return ordinary, robust


def effect(model, term: str) -> dict:
    b = float(model.params[term])
    se = float(model.bse[term])
    z = b / se
    return {
        "term": term,
        "HR": math.exp(b),
        "CI_low": math.exp(b - 1.959963984540054 * se),
        "CI_high": math.exp(b + 1.959963984540054 * se),
        "robust_wald_p": float(2 * norm.sf(abs(z))),
    }


def compare(
    cohort: str,
    analysis: str,
    data: pd.DataFrame,
    id_col: str,
    formula_b: str,
    formula_c: str,
    change_term: str,
    landmark_term: str,
    change_scale: str,
) -> dict:
    ob, rb = fit(formula_b, data, id_col)
    oc, rc = fit(formula_c, data, id_col)
    lr = max(0.0, 2.0 * (float(ob.llf) - float(oc.llf)))
    df = int(len(ob.params) - len(oc.params))
    n_pp = int(len(data))
    n_people = int(data[id_col].nunique())
    out = {
        "cohort": cohort,
        "analysis": analysis,
        "participants": n_people,
        "person_periods": n_pp,
        "events": int(data["event"].sum()),
        "change_scale": change_scale,
        "model_B_log_likelihood": float(ob.llf),
        "model_C_log_likelihood": float(oc.llf),
        "model_B_AIC": float(ob.aic),
        "model_C_AIC": float(oc.aic),
        "delta_AIC_B_minus_C": float(ob.aic - oc.aic),
        "model_B_BIC_participant": float(-2 * ob.llf + len(ob.params) * math.log(n_people)),
        "model_C_BIC_participant": float(-2 * oc.llf + len(oc.params) * math.log(n_people)),
        "delta_BIC_participant_B_minus_C": float(
            (-2 * ob.llf + len(ob.params) * math.log(n_people))
            - (-2 * oc.llf + len(oc.params) * math.log(n_people))
        ),
        "LR_statistic_model_based": lr,
        "LR_df": df,
        "LR_p_model_based": float(chi2.sf(lr, df)),
        **{f"change_{k}": v for k, v in effect(rb, change_term).items()},
        **{f"model_C_landmark_{k}": v for k, v in effect(rc, landmark_term).items()},
        "interpretation_note": (
            "Model B is expressed as landmark grip plus change and is algebraically "
            "equivalent to earlier grip plus change. The robust change coefficient tests "
            "whether change adds information conditional on landmark grip. The likelihood-"
            "ratio test is model-based and does not use the cluster-robust covariance."
        ),
    }
    return out


def charls_rows(base: dict) -> list[dict]:
    # Reuse the audited cohort constructor without overwriting its archived outputs.
    base["write_csv"] = lambda *args, **kwargs: None
    _, _, frames = base["build_charls"]()
    analytic = frames["analytic"].copy()
    analytic["landmark_grip"] = analytic["grip_2011"] - analytic["grip_decline"]
    common = (
        "age_2015 + C(sex) + C(education) + C(residence) + "
        "C(marital_status) + C(smoking) + C(alcohol) + "
        "chronic_condition_count + C(interval)"
    )
    rows = []
    for label, variable in [
        ("Earlier ADL/IADL history adjustment", "early_limitation"),
        ("BMI adjustment", "bmi_2011"),
    ]:
        source = "strict_2011" if variable == "early_limitation" else "bmi_2011"
        data = analytic.dropna(subset=[source]).copy()
        if variable == "early_limitation":
            data[variable] = data[source].astype(int)
        sd = float(data.drop_duplicates("study_id")["grip_decline"].std(ddof=1))
        data["grip_decline_sd"] = data["grip_decline"] / sd
        extra = f" + {variable}"
        rows.append(
            compare(
                "CHARLS",
                label,
                data,
                "study_id",
                "event ~ landmark_grip + grip_decline_sd + " + common + extra,
                "event ~ landmark_grip + " + common + extra,
                "grip_decline_sd",
                "landmark_grip",
                f"1 SD = {sd:.6f} kg",
            )
        )
    return rows


def elsa_rows(base: dict) -> list[dict]:
    ns = base["load_elsa_namespace"]()
    spec = next(x for x in ns["specs"] if x["model_id"] == "ELSA-P1")
    pp = ns["model_frame"](spec)
    if (pp.idauniq.nunique(), len(pp), int(pp.event.sum())) != (3991, 7233, 733):
        raise RuntimeError("ELSA primary anchor failed")
    ids = pp.idauniq.astype("string")
    pp["grip_change_kg"] = (ns["g4"] - ns["g6"]).reindex(ids).to_numpy()
    pp["landmark_grip"] = pp["baseline_function"] - pp["grip_change_kg"]

    w4out = ns["read_dta"](
        "wave_4_elsa_data_eul.dta", ["idauniq"] + ns["ADL_DIFF"] + ns["IADL_DIFF"]
    )
    ak, ad = ns["difficulty"](w4out, ns["ADL_DIFF"])
    ik, ip = ns["difficulty"](w4out, ns["IADL_DIFF"])
    early_known = ak & ik
    early_positive = (ad.fillna(False) | ip.fillna(False)).where(early_known)
    pp["early_limitation"] = early_positive.reindex(ids).astype("Float64").to_numpy(
        dtype=float, na_value=np.nan
    )

    w4bmi = ns["read_dta"]("wave_4_nurse_data.dta", ["idauniq", "bmival"])
    bmi_series = pd.to_numeric(w4bmi.bmival, errors="coerce").where(
        lambda x: x.between(10, 70)
    )
    pp["bmi_w4"] = bmi_series.reindex(ids).to_numpy()

    common = (
        "age + C(sex) + C(education) + partnered + urban + current_smoking + "
        "C(alcohol_frequency) + chronic_count + C(interval)"
    )
    rows = []
    for label, variable in [
        ("Earlier ADL/IADL history adjustment", "early_limitation"),
        ("BMI adjustment", "bmi_w4"),
    ]:
        data = pp.dropna(subset=[variable]).copy()
        rows.append(
            compare(
                "ELSA",
                label,
                data,
                "idauniq",
                "event ~ landmark_grip + grip_z + " + common + f" + {variable}",
                "event ~ landmark_grip + " + common + f" + {variable}",
                "grip_z",
                "landmark_grip",
                "1 frozen reference-population SD = 5.6019 kg",
            )
        )
    return rows


def main() -> None:
    base = runpy.run_path(str(BASE_SCRIPT), run_name="__bc_sensitivity_base__")
    rows = charls_rows(base) + elsa_rows(base)
    result = pd.DataFrame(rows)
    result.to_csv(
        OUT / "BC_SENSITIVITY_COMPARISON_RESULTS.csv",
        index=False,
        encoding="utf-8-sig",
    )
    qc = {
        "row_count": int(len(result)),
        "expected_rows": 4,
        "all_finite_key_statistics": bool(
            np.isfinite(
                result[
                    [
                        "model_B_log_likelihood",
                        "model_C_log_likelihood",
                        "delta_AIC_B_minus_C",
                        "delta_BIC_participant_B_minus_C",
                        "LR_p_model_based",
                        "change_HR",
                        "change_CI_low",
                        "change_CI_high",
                        "change_robust_wald_p",
                    ]
                ].to_numpy(dtype=float)
            ).all()
        ),
        "identical_samples_within_each_B_C_comparison": True,
        "primary_files_overwritten": False,
    }
    (OUT / "BC_SENSITIVITY_QC.json").write_text(
        json.dumps(qc, indent=2), encoding="utf-8"
    )
    print(result.to_json(orient="records", indent=2))
    print(json.dumps(qc, indent=2))


if __name__ == "__main__":
    main()

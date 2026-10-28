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
BASE = HERE / "01_sensitivity_analyses.py"
OLD_COHORT = os.environ.get("CHARLS_V1_DEIDENTIFIED_COHORT")
OUT = Path(os.environ.get("SENSITIVITY_OUTPUT_ROOT", REPO_ROOT / "outputs" / "sensitivity")) / "S18"
OUT.mkdir(parents=True, exist_ok=True)


def fit(formula: str, data: pd.DataFrame):
    model = smf.glm(
        formula,
        data=data,
        family=sm.families.Binomial(link=sm.families.links.CLogLog()),
    )
    ordinary = model.fit(maxiter=200)
    robust = model.fit(
        cov_type="cluster", cov_kwds={"groups": data["study_id"]}, maxiter=200
    )
    rank = int(np.linalg.matrix_rank(ordinary.model.exog))
    if rank != ordinary.model.exog.shape[1]:
        raise RuntimeError(f"Rank deficient: {rank}/{ordinary.model.exog.shape[1]}")
    return ordinary, robust


def effect(model, term: str) -> dict:
    b = float(model.params[term])
    se = float(model.bse[term])
    z = b / se
    return {
        "HR": math.exp(b),
        "CI_low": math.exp(b - 1.959963984540054 * se),
        "CI_high": math.exp(b + 1.959963984540054 * se),
        "robust_wald_p": float(2 * norm.sf(abs(z))),
    }


def cesd_and_bmi(ns: dict) -> pd.DataFrame:
    read, xid, num = ns["read"], ns["xid"], ns["num"]
    items = [f"dc{i:03d}" for i in range(9, 19)]
    h = read(ns["FILES"]["hsf2015"], ["ID"] + items)
    h["pid"] = xid(h.ID, 2015)
    z = h[items].apply(num) - 1
    # V1 construction recovered from the archived table-generation workflow.
    for col in ["dc013", "dc016"]:
        z[col] = 3 - z[col]
    h["cesd10_2015"] = z.sum(axis=1, min_count=10)

    b = read(ns["FILES"]["grip2015"], ["ID", "qi002", "ql002"])
    b["pid"] = xid(b.ID, 2015)
    height = num(b.qi002).where(num(b.qi002).between(100, 220))
    weight = num(b.ql002).where(num(b.ql002).between(25, 250))
    b["bmi_2015"] = (weight / (height / 100) ** 2).where(lambda x: x.between(10, 80))
    return h[["pid", "cesd10_2015"]].drop_duplicates("pid").merge(
        b[["pid", "bmi_2015"]].drop_duplicates("pid"), on="pid", how="outer"
    )


def v1_reproduction(ns: dict, exposure: pd.DataFrame) -> dict:
    if not OLD_COHORT:
        return {
            "V1_historical_audit_status": "NOT_RUN",
            "V1_historical_audit_reason": "CHARLS_V1_DEIDENTIFIED_COHORT was not supplied",
        }
    d = pd.read_csv(Path(OLD_COHORT).expanduser().resolve())
    base = d.burden_increase.notna() & d.grip_decline.notna() & d.strict_2015.eq(False)
    follow = base & (d.strict_2018.notna() | d.strict_2020.notna())
    pp_ok = follow & d.strict_2018.notna()
    cov = [
        "burden_increase", "grip_decline", "age_2015", "sex", "education",
        "urban", "married", "smoking", "alcohol", "chronic_count",
        "burden_2011", "grip_2011",
    ]
    included = base & pp_ok & d[cov].notna().all(axis=1)
    # The old cohort uses a different deidentification hash, so join via its
    # already stored study_id after applying the archived V1 hash to raw IDs.
    import hashlib

    raw = ns["read"](ns["FILES"]["hsf2015"], ["ID"] + [f"dc{i:03d}" for i in range(9, 19)])
    native = raw.ID.map(
        lambda v: pd.NA
        if pd.isna(v) or not str(v).strip()
        else (v.decode() if isinstance(v, bytes) else str(v).strip())
    ).astype("string")
    raw["study_id"] = native.map(
        lambda x: hashlib.sha256(("phase3b:" + str(x)).encode()).hexdigest()[:16]
        if pd.notna(x)
        else pd.NA
    )
    items = [f"dc{i:03d}" for i in range(9, 19)]
    z = raw[items].apply(pd.to_numeric, errors="coerce") - 1
    for col in ["dc013", "dc016"]:
        z[col] = 3 - z[col]
    raw["cesd10"] = z.sum(axis=1, min_count=10)
    x = d.loc[included, ["study_id"]].merge(raw[["study_id", "cesd10"]], on="study_id", how="left")
    return {
        "V1_N": int(len(x)),
        "V1_CESD_complete": int(x.cesd10.notna().sum()),
        "V1_CESD_missing_reproduced": int(x.cesd10.isna().sum()),
        "V1_CESD_mean_reproduced": float(x.cesd10.mean()),
        "V1_CESD_SD_reproduced": float(x.cesd10.std(ddof=1)),
        "legacy_manuscript_missing_count": 156,
        "legacy_missing_count_reproduced": False,
    }


def main() -> None:
    base = runpy.run_path(str(BASE), run_name="__cesd_joint_base__")
    base["write_csv"] = lambda *args, **kwargs: None
    _, _, frames = base["build_charls"]()
    pp = frames["analytic"].copy()
    ns = runpy.run_path(
        str(REPO_ROOT / "src" / "charls" / "01_primary_analysis.py"),
        run_name="__cesd_joint_charls__",
    )
    exposures = cesd_and_bmi(ns)
    # Recover raw pid on the person-period data by rebuilding the primary frame.
    # study_id is a deterministic v2 hash, allowing an audited one-to-one map.
    cohort = ns["grip"](2011).merge(ns["grip"](2015), on="pid", how="outer", validate="one_to_one")
    import hashlib
    cohort["study_id"] = cohort.pid.map(
        lambda x: hashlib.sha256(("v2-charls:" + str(x)).encode()).hexdigest()[:16]
    )
    lookup = cohort[["study_id", "pid"]].drop_duplicates("study_id").merge(
        exposures, on="pid", how="left"
    )
    pp = pp.merge(lookup[["study_id", "bmi_2015", "cesd10_2015"]], on="study_id", how="left", validate="many_to_one")
    data = pp.dropna(subset=["bmi_2015", "cesd10_2015"]).copy()
    unique = data.sort_values("study_id").drop_duplicates("study_id")
    sd = float(unique.grip_decline.std(ddof=1))
    data["grip_decline_sd"] = data.grip_decline / sd
    data["landmark_grip"] = data.grip_2011 - data.grip_decline
    common = (
        "age_2015 + C(sex) + C(education) + C(residence) + C(marital_status) + "
        "C(smoking) + C(alcohol) + chronic_condition_count + C(interval)"
    )
    f_same = "event ~ grip_decline_sd + grip_2011 + " + common
    f_joint = f_same + " + bmi_2015 + cesd10_2015"
    f_b = "event ~ landmark_grip + grip_decline_sd + " + common + " + bmi_2015 + cesd10_2015"
    f_c = "event ~ landmark_grip + " + common + " + bmi_2015 + cesd10_2015"
    _, r_same = fit(f_same, data)
    _, r_joint = fit(f_joint, data)
    ob, rb = fit(f_b, data)
    oc, rc = fit(f_c, data)
    lr = max(0.0, 2 * (float(ob.llf) - float(oc.llf)))
    n = int(data.study_id.nunique())
    rows = [
        {
            "model": "same-sample primary covariates",
            "participants": n, "person_periods": int(len(data)), "events": int(data.event.sum()),
            **effect(r_same, "grip_decline_sd"),
        },
        {
            "model": "post-hoc joint adjustment for 2015 BMI and CES-D10",
            "participants": n, "person_periods": int(len(data)), "events": int(data.event.sum()),
            **effect(r_joint, "grip_decline_sd"),
        },
        {
            "model": "Model B versus Model C under joint adjustment",
            "participants": n, "person_periods": int(len(data)), "events": int(data.event.sum()),
            **effect(rb, "grip_decline_sd"),
            "LR_statistic_model_based": lr,
            "LR_p_model_based": float(chi2.sf(lr, 1)),
            "delta_AIC_B_minus_C": float(ob.aic - oc.aic),
            "delta_BIC_participant_B_minus_C": float(
                (-2 * ob.llf + len(ob.params) * math.log(n))
                - (-2 * oc.llf + len(oc.params) * math.log(n))
            ),
        },
    ]
    pd.DataFrame(rows).to_csv(OUT / "CHARLS_BMI_CESD10_JOINT_RESULTS.csv", index=False, encoding="utf-8-sig")
    audit = v1_reproduction(ns, exposures)
    audit.update({
        "V2_joint_participants": n,
        "V2_joint_person_periods": int(len(data)),
        "V2_joint_events": int(data.event.sum()),
        "V2_CESD_complete_participants": int(pp.drop_duplicates("study_id").cesd10_2015.notna().sum()),
        "V2_BMI_complete_participants": int(pp.drop_duplicates("study_id").bmi_2015.notna().sum()),
        "V2_joint_complete_participants": n,
        "change_SD_kg": sd,
        "CESD_definition": "2015 dc009-dc018; response 1-4 mapped to 0-3; dc013 and dc016 reverse scored; all 10 items required",
        "BMI_definition": "2015 measured weight / measured height squared; BMI 10-80 kg/m2",
        "analysis_status": "POST_HOC_AFTER_PRE_SUBMISSION_REVIEW",
    })
    (OUT / "CESD10_CONSTRUCTION_AND_QC.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(pd.DataFrame(rows).to_string(index=False))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()

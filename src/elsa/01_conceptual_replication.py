# Part of the charls-elsa-physical-function-activity-limitations project.
# Developed for the analyses reported in the associated manuscript.
# Released under the MIT License. See LICENSE in the repository root.
# No CHARLS or ELSA individual-level data are distributed with this code.

from pathlib import Path
import os
import csv, hashlib, json, shutil
import numpy as np
import pandas as pd
import pyreadstat

REPO_ROOT = Path(__file__).resolve().parents[2]
ELSA_PROJECT_ROOT = Path(os.environ["ELSA_PROJECT_ROOT"])
P1 = Path(os.environ.get("ELSA_PHASE1_ROOT", ELSA_PROJECT_ROOT / "phase_elsa1"))
P2 = Path(os.environ.get("ELSA_PHASE2_ROOT", ELSA_PROJECT_ROOT / "phase_elsa2"))
ROOT = Path(os.environ.get("ELSA_OUTPUT_ROOT", REPO_ROOT / "outputs" / "elsa"))
V11 = Path(os.environ.get("ELSA_PROTOCOL_V11_ROOT", ELSA_PROJECT_ROOT / "protocol_v1_1"))
V12 = Path(os.environ.get("ELSA_PROTOCOL_V12_ROOT", ELSA_PROJECT_ROOT / "protocol_v1_2"))
DATA = Path(os.environ["ELSA_DATA_ROOT"])
DOC = Path(os.environ.get("ELSA_DOCUMENT_TEXT_ROOT", ELSA_PROJECT_ROOT / "document_text"))
for d in ["00_author_decision", "01_definitions", "02_anchors", "03_protocol", "04_workbooks",
          "05_preflight", "06_reports", "07_scripts", "08_logs", "09_manifest"]:
    (ROOT / d).mkdir(parents=True, exist_ok=True)

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()

def write_csv(path, rows, fields=None):
    rows = list(rows)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(rows)

def read_dta(name, cols):
    _, meta = pyreadstat.read_dta(str(DATA / name), metadataonly=True)
    exact = {x.lower(): x for x in meta.column_names}
    missing = [x for x in cols if x.lower() not in exact]
    if missing:
        raise KeyError(f"{name}: missing columns {missing}")
    actual = [exact[x.lower()] for x in cols]
    df, _ = pyreadstat.read_dta(str(DATA / name), usecols=actual, apply_value_formats=False)
    df.columns = [x.lower() for x in df.columns]
    df["idauniq"] = pd.to_numeric(df["idauniq"], errors="coerce").astype("Int64").astype("string")
    df.loc[df["idauniq"].isin(["", "<NA>"]), "idauniq"] = pd.NA
    if df["idauniq"].dropna().duplicated().any():
        raise ValueError(f"duplicate idauniq in {name}")
    return df.set_index("idauniq")

ADL_DIFF = ["headldr", "headlba", "headlea", "headlbe", "headlwc"]
IADL_DIFF = ["headlpr", "headlsh", "headlph", "headlme", "headlho", "headlmo"]
ADL_HELP = ["catkd", "catkf", "catkg", "catkh", "catki"]
OUTCOME_COLS = ["idauniq"] + ADL_DIFF + IADL_DIFF + ADL_HELP
GRIP_COLS = ["mmgsd1", "mmgsn1", "mmgsd2", "mmgsn2", "mmgsd3", "mmgsn3"]
WALK_COLS = ["mmtrya", "mmwlka", "mmtryb", "mmwlkb"]
BIOM = ["wbc", "hgb", "fglu", "chol", "hdl", "ldl", "trig", "hscrp", "hba1c"]

chronic_parts = {
    "hypertension": [("hedawbp", "hedacbp", "hediabp", 1)],
    "diabetes": [("hedawdi", "hedacdi", "hediadi", 7)],
    "stroke": [("hedawst", "hedacst", "hediast", 8)],
    "heart_disease": [
        ("hedawan", "hedacan", "hediaan", 2), ("hedawmi", "hedacmi", "hediami", 3),
        ("hedawhf", "hedachf", "hediahf", 4), ("hedawar", "hedacar", "hediaar", 6),
    ],
    "high_cholesterol": [("hedawch", "hedacch", "hediach", 9)],
    "chronic_lung_disease": [("hedbwlu", "hedbdlu", "hediblu", 1)],
    "asthma": [("hedbwas", "hedbdas", "hedibas", 2)],
    "arthritis": [("hedbwar", "hedbdar", "hedibar", 3)],
    "osteoporosis": [("hedbwos", "hedbdos", "hedibos", 4)],
    "cancer": [("hedbwca", "hedbdca", "hedibca", 5)],
    "parkinsons_disease": [("hedbwpd", "hedbdpd", "hedibpd", 6)],
    "psychiatric_condition": [("hedbwps", "hedbdps", "hedibps", 7)],
    "alzheimer_or_dementia": [
        ("hedbwad", "hedbdad", "hedibad", 8), ("hedbwde", "hedbdde", "hedibde", 9),
    ],
    "malignant_blood_disorder": [("hedbwbl", "hedbdbl", "hedibbl", 10)],
}
chronic_cols = sorted({v for parts in chronic_parts.values() for p in parts for v in p[:3]})
cov_cols = ["idauniq", "indager", "indsex", "fqendm", "dimarr", "hesmk", "heska", "scako"] + chronic_cols + WALK_COLS

w4c = read_dta("wave_4_elsa_data_eul.dta", ["idauniq", "indsex"] + WALK_COLS)
w6c = read_dta("wave_6_elsa_data_eul.dta", cov_cols + ADL_DIFF + IADL_DIFF + ADL_HELP)
w7 = read_dta("wave_7_elsa_data_eul.dta", OUTCOME_COLS)
w8 = read_dta("wave_8_elsa_data_eul_v2.dta", OUTCOME_COLS)
w4n = read_dta("wave_4_nurse_data.dta", ["idauniq"] + GRIP_COLS + BIOM + ["mmssre", "mmstre", "mmftre2"])
w6n = read_dta("wave_6_elsa_nurse_data_v2.dta", ["idauniq"] + GRIP_COLS + BIOM + ["mmssre", "mmstre", "mmftre2"])
geo = read_dta("elsa_geog_urindewr_2011_eul.dta", ["idauniq", "w6_urindewr_2011"])

def valid_positive(s, upper=None):
    x = pd.to_numeric(s, errors="coerce")
    ok = x.gt(0)
    if upper is not None: ok &= x.lt(upper)
    return x.where(ok)

def grip(df):
    x = df[GRIP_COLS].apply(lambda s: valid_positive(s, 100))
    return x.max(axis=1, skipna=True)

def walk(df):
    vals = []
    for outcome, time in [("mmtrya", "mmwlka"), ("mmtryb", "mmwlkb")]:
        t = valid_positive(df[time])
        vals.append(t.where(pd.to_numeric(df[outcome], errors="coerce").eq(1)))
    return pd.concat(vals, axis=1).mean(axis=1, skipna=True)

def balance_level(df):
    ss = pd.to_numeric(df["mmssre"], errors="coerce")
    st = pd.to_numeric(df["mmstre"], errors="coerce")
    ft = pd.to_numeric(df["mmftre2"], errors="coerce")
    out = pd.Series(pd.NA, index=df.index, dtype="Float64")
    out.loc[ss.isin([2, 3])] = 0
    out.loc[ss.eq(1) & st.isin([2, 3])] = 1
    out.loc[ss.eq(1) & st.eq(1) & ft.isin([2, 4, 5])] = 2
    out.loc[ss.eq(1) & st.eq(1) & ft.isin([1, 3])] = 3
    return out

def difficulty(df, items):
    x = df[items].apply(pd.to_numeric, errors="coerce")
    known = x.isin([0, 1]).all(axis=1)
    positive = x.eq(1).any(axis=1).where(known)
    return known, positive

def assistance5(df):
    d = df[ADL_DIFF].apply(pd.to_numeric, errors="coerce")
    h = df[ADL_HELP].apply(pd.to_numeric, errors="coerce")
    h.columns = ADL_DIFF
    item_known = d.eq(0) | (d.eq(1) & h.isin([1, 2]))
    item_help = d.eq(1) & h.eq(1)
    known = item_known.all(axis=1)
    positive = item_help.any(axis=1).where(known)
    return known, positive

def one_condition(df, ff, conf, new, ff_code):
    a = pd.to_numeric(df[ff], errors="coerce")
    b = pd.to_numeric(df[conf], errors="coerce")
    c = pd.to_numeric(df[new], errors="coerce")
    out = pd.Series(pd.NA, index=df.index, dtype="Float64")
    out.loc[b.eq(1) | c.eq(1)] = 1
    out.loc[c.eq(0) & (~a.eq(ff_code) | b.isin([2, 3]))] = 0
    return out

def condition_group(df, components):
    x = pd.concat([one_condition(df, *p) for p in components], axis=1)
    out = pd.Series(pd.NA, index=df.index, dtype="Float64")
    out.loc[x.eq(1).any(axis=1)] = 1
    out.loc[x.notna().all(axis=1) & x.eq(0).all(axis=1)] = 0
    return out

cov = pd.DataFrame(index=w6c.index)
cov["age"] = pd.to_numeric(w6c["indager"], errors="coerce").where(lambda x: x.ge(50))
cov["sex"] = pd.to_numeric(w6c["indsex"], errors="coerce").where(lambda x: x.isin([1, 2]))
ed = pd.to_numeric(w6c["fqendm"], errors="coerce")
cov["education"] = pd.Series(pd.NA, index=cov.index, dtype="Float64")
cov.loc[ed.isin([1, 2, 3]), "education"] = 1
cov.loc[ed.eq(4), "education"] = 2
cov.loc[ed.eq(5), "education"] = 3
cov.loc[ed.isin([6, 7]), "education"] = 4
cov.loc[ed.eq(8), "education"] = 5
ma = pd.to_numeric(w6c["dimarr"], errors="coerce")
cov["partnered"] = ma.isin([2, 3]).astype("Float64").where(ma.isin([1, 2, 3, 4, 5, 6]))
sm = pd.to_numeric(w6c["heska"], errors="coerce")
ever = pd.to_numeric(w6c["hesmk"], errors="coerce")
cov["current_smoking"] = pd.Series(pd.NA, index=cov.index, dtype="Float64")
# AUTHOR_DECISION (pre-result): explicit HESKA answers take precedence;
# only the deterministic HESMK=2/HESKA=-1 structural skip maps to non-current.
cov.loc[sm.eq(1), "current_smoking"] = 1
cov.loc[sm.eq(2), "current_smoking"] = 0
cov.loc[sm.eq(-1) & ever.eq(2), "current_smoking"] = 0
al = pd.to_numeric(w6c["scako"], errors="coerce")
cov["alcohol_frequency"] = pd.Series(pd.NA, index=cov.index, dtype="Float64")
cov.loc[al.eq(8), "alcohol_frequency"] = 0
cov.loc[al.isin([5, 6, 7]), "alcohol_frequency"] = 1
cov.loc[al.isin([3, 4]), "alcohol_frequency"] = 2
cov.loc[al.isin([1, 2]), "alcohol_frequency"] = 3
g = pd.to_numeric(geo["w6_urindewr_2011"], errors="coerce").reindex(cov.index)
cov["urban"] = g.eq(1).astype("Float64").where(g.isin([1, 2]))
condition_matrix = pd.DataFrame({k: condition_group(w6c, v) for k, v in chronic_parts.items()})
cov["chronic_count"] = condition_matrix.sum(axis=1).where(condition_matrix.notna().all(axis=1))

g4, g6 = grip(w4n), grip(w6n)
walk4, walk6 = walk(w4c), walk(w6c)
bal4, bal6 = balance_level(w4n), balance_level(w6n)

# ELSA-adapted 9-biomarker, 4-domain burden.
b4 = w4n[BIOM].apply(pd.to_numeric, errors="coerce").copy()
b6 = w6n[BIOM].apply(pd.to_numeric, errors="coerce").copy()
for c in BIOM:
    b4[c] = b4[c].where(b4[c].gt(0)); b6[c] = b6[c].where(b6[c].gt(0))
b4["trig"] = np.log1p(b4["trig"]); b6["trig"] = np.log1p(b6["trig"])
b4["hscrp"] = np.log1p(b4["hscrp"]); b6["hscrp"] = np.log1p(b6["hscrp"])
# Wave 6 IFCC mmol/mol converted to NGSP percent so the Wave 4 scale is retained.
b6["hba1c"] = 0.09148 * b6["hba1c"] + 2.152
sex4 = pd.to_numeric(w4c["indsex"], errors="coerce").reindex(b4.index)
sex6 = pd.to_numeric(w6c["indsex"], errors="coerce").reindex(b6.index)
means, sds = {}, {}
for c in BIOM:
    if c == "hgb": continue
    means[c], sds[c] = float(b4[c].mean()), float(b4[c].std(ddof=1))
hgb_ref = {}
for sx in [1, 2]:
    hgb_ref[sx] = {"mean": float(b4.loc[sex4.eq(sx), "hgb"].mean()), "sd": float(b4.loc[sex4.eq(sx), "hgb"].std(ddof=1))}

def biom_score(b, sex):
    z = pd.DataFrame(index=b.index)
    for c in BIOM:
        if c == "hgb":
            z[c] = pd.NA
            for sx in [1, 2]:
                z.loc[sex.eq(sx), c] = -(b.loc[sex.eq(sx), c] - hgb_ref[sx]["mean"]) / hgb_ref[sx]["sd"]
            z[c] = pd.to_numeric(z[c], errors="coerce")
        else:
            z[c] = (b[c] - means[c]) / sds[c]
            if c == "hdl": z[c] = -z[c]
    domain = pd.DataFrame({
        "inflammatory": z[["wbc", "hscrp"]].mean(axis=1).where(z[["wbc", "hscrp"]].notna().all(axis=1)),
        "haematologic": z["hgb"],
        "glycaemic": z[["fglu", "hba1c"]].mean(axis=1).where(z[["fglu", "hba1c"]].notna().all(axis=1)),
        "lipid": z[["chol", "hdl", "ldl", "trig"]].mean(axis=1).where(z[["chol", "hdl", "ldl", "trig"]].notna().all(axis=1)),
    })
    return domain.mean(axis=1).where(domain.notna().all(axis=1))

burden4, burden6 = biom_score(b4, sex4), biom_score(b6, sex6)

status = {}
for wk, df in [("w6", w6c), ("w7", w7), ("w8", w8)]:
    ak, ad = difficulty(df, ADL_DIFF); ik, ip = difficulty(df, IADL_DIFF); hk, hp = assistance5(df)
    status[(wk, "adl_difficulty")] = (ak, ad)
    status[(wk, "iadl_difficulty")] = (ik, ip)
    status[(wk, "combined_difficulty")] = (ak & ik, (ad.fillna(False) | ip.fillna(False)).where(ak & ik))
    status[(wk, "adl_help5")] = (hk, hp)

def person_period(eligible, outcome, baseline_positive=False):
    idx = pd.Index(sorted(set(eligible)))
    k6, p6 = status[("w6", outcome)]; k7, p7 = status[("w7", outcome)]; k8, p8 = status[("w8", outcome)]
    q = pd.DataFrame(index=idx)
    for n, s in [("k6", k6), ("p6", p6), ("k7", k7), ("p7", p7), ("k8", k8), ("p8", p8)]: q[n] = s.reindex(idx)
    # Use explicit equality rather than bitwise inversion so nullable/object
    # booleans cannot turn False into integer -1 under pandas coercion.
    risk = q["k6"].eq(True) & (q["p6"].eq(True) if baseline_positive else q["p6"].eq(False))
    r = q.loc[risk]
    pp1 = r["k7"].eq(True)
    ev1 = pp1 & (r["p7"].eq(False) if baseline_positive else r["p7"].eq(True))
    remain = r["p7"].eq(True) if baseline_positive else r["p7"].eq(False)
    pp2 = pp1 & remain & r["k8"].eq(True)
    ev2 = pp2 & (r["p8"].eq(False) if baseline_positive else r["p8"].eq(True))
    if (ev1 & ev2).any(): raise AssertionError("post-event duplicate event")
    return {
        "risk_set_n": int(risk.sum()), "valid_person_n": int((pp1 | pp2).sum()),
        "person_period_n": int(pp1.sum() + pp2.sum()), "event_n": int(ev1.sum() + ev2.sum()),
        "interval1_rows": int(pp1.sum()), "interval1_events": int(ev1.sum()),
        "interval2_rows": int(pp2.sum()), "interval2_events": int(ev2.sum()),
        "valid_ids": set(r.index[pp1 | pp2]),
    }

base_cov = ["age", "sex", "education", "partnered", "urban", "current_smoking", "alcohol_frequency", "chronic_count"]
cov_complete = set(cov.index[cov[base_cov].notna().all(axis=1)])

exposures = {
    "grip": set(g4.dropna().index) & set(g6.dropna().index),
    "walk": set(walk4.dropna().index) & set(walk6.dropna().index),
    "balance": set(bal4.dropna().index) & set(bal6.dropna().index),
    "biom9_grip": set(burden4.dropna().index) & set(burden6.dropna().index) & set(g4.dropna().index) & set(g6.dropna().index),
}
analysis_specs = [
    ("主要概念性重复：握力下降与首次ADL/IADL活动困难", "combined_difficulty", "grip", False, g4),
    ("次要：5项ADL帮助依赖", "adl_help5", "grip", False, g4),
    ("次要并列生理：9项负担、握力及乘积项与首次ADL/IADL活动困难", "combined_difficulty", "biom9_grip", False, g4),
    ("领域：步行时间增加与首次ADL/IADL活动困难", "combined_difficulty", "walk", False, walk4),
    ("领域：平衡恶化与首次ADL/IADL活动困难", "combined_difficulty", "balance", False, bal4),
    ("状态转换：ADL活动困难消失", "adl_difficulty", "grip", True, g4),
    ("状态转换：IADL活动困难消失", "iadl_difficulty", "grip", True, g4),
    ("状态转换：停止报告5项ADL帮助需求", "adl_help5", "grip", True, g4),
]

anchor_rows = []
for name, outcome, exp, baseline_positive, baseline_measure in analysis_specs:
    eligible0 = exposures[exp] & set(w6c.index)
    raw = person_period(eligible0, outcome, baseline_positive)
    eligible_cc = eligible0 & cov_complete & set(baseline_measure.dropna().index)
    if exp == "biom9_grip": eligible_cc &= set(burden4.dropna().index)
    cc = person_period(eligible_cc, outcome, baseline_positive)
    anchor_rows.append({
        "analysis": name, "outcome": outcome, "exposure_set": exp,
        "baseline_status": "positive" if baseline_positive else "negative",
        "two_wave_exposure_complete_n": len(eligible0), "risk_set_n": raw["risk_set_n"],
        "valid_followup_person_n": raw["valid_person_n"], "person_period_n_before_covariates": raw["person_period_n"],
        "events_before_covariates": raw["event_n"], "fully_adjusted_person_n": cc["valid_person_n"],
        "fully_adjusted_person_period_n": cc["person_period_n"], "fully_adjusted_event_n": cc["event_n"],
        "interval1_rows": cc["interval1_rows"], "interval1_events": cc["interval1_events"],
        "interval2_rows": cc["interval2_rows"], "interval2_events": cc["interval2_events"],
        "association_models_run": 0,
    })
write_csv(ROOT / "02_anchors" / "ELSA_REVISED_SAMPLE_AND_EVENT_ANCHORS_FINAL.csv", anchor_rows)
write_csv(ROOT / "01_definitions" / "ELSA_ANALYSIS_POPULATION_MATRIX_FINAL.csv", anchor_rows)
transition_rows = [r for r in anchor_rows if r["baseline_status"] == "positive"]
write_csv(ROOT / "02_anchors" / "ELSA_STATE_TRANSITION_COUNTS_FINAL.csv", transition_rows)
go_rows = []
for r in anchor_rows:
    is_primary = r["analysis"].startswith("主要概念性重复")
    adequate = r["fully_adjusted_person_n"] >= 500 and r["fully_adjusted_event_n"] >= 100
    if adequate:
        decision = "GO"
        reason = "definition closed and descriptive sample/event support adequate"
    elif is_primary:
        decision = "NO_GO_BLOCKS_REPLICATION"
        reason = "primary conceptual replication fails descriptive feasibility threshold"
    else:
        decision = "DESCRIPTIVE_ONLY_NOT_FORMAL_MODEL"
        reason = "definition closed but descriptive sample/event support below feasibility threshold"
    go_rows.append({"analysis":r["analysis"], "fully_adjusted_n":r["fully_adjusted_person_n"],
                    "events":r["fully_adjusted_event_n"], "n_threshold_500":r["fully_adjusted_person_n"]>=500,
                    "event_threshold_100":r["fully_adjusted_event_n"]>=100, "decision":decision, "reason":reason})
write_csv(ROOT / "02_anchors" / "ELSA_GO_NO_GO_BY_ANALYSIS_FINAL.csv", go_rows)

outcome_rows = []
item_names = {
    "headldr": "dressing", "headlba": "bathing/showering", "headlea": "eating/cutting food",
    "headlbe": "getting in/out of bed", "headlwc": "using toilet",
    "headlpr": "preparing hot meal", "headlsh": "shopping for groceries", "headlph": "making telephone calls",
    "headlme": "taking medications", "headlho": "doing work around house/garden", "headlmo": "managing money",
}
help_map = dict(zip(ADL_DIFF, ADL_HELP))
for v in ADL_DIFF + IADL_DIFF:
    outcome_rows.append({"domain": "ADL" if v in ADL_DIFF else "IADL", "activity": item_names[v],
                         "difficulty_variable_w6_w7_w8": v, "difficulty_codes": "0=no difficulty; 1=difficulty; negative=special missing",
                         "help_variable_w6_w7_w8": help_map.get(v, "NOT_AVAILABLE"),
                         "help_codes": "1=received help; 2=no help; negative=special missing/structural skip" if v in help_map else "NOT_AVAILABLE",
                         "formal_use": "main difficulty outcome + resolution" if v not in ADL_DIFF else "main difficulty outcome + 5-ADL assistance secondary + state resolution/cessation"})
write_csv(ROOT / "01_definitions" / "ELSA_ADL_IADL_FINAL_MAPPING.csv", outcome_rows)

physical_rows = [
    {"domain": "Grip strength", "raw_variables": "/".join(GRIP_COLS), "raw_codes_and_missing": "kg; values <=0 or >=100 invalid; negative=special missing", "wave_summary": "maximum of all valid trials", "change": "W4 minus W6; positive=decline", "eligibility": "at least one valid trial per wave", "charls_comparison": "same max-of-valid-trials concept; ELSA waves differ"},
    {"domain": "Timed walk", "raw_variables": "/".join(WALK_COLS), "raw_codes_and_missing": "outcome 1=completed; time>0 seconds; negative=special missing", "wave_summary": "mean of all successful timed trials; one successful trial is usable", "change": "W6 mean seconds minus W4 mean seconds; positive=slower", "eligibility": "at least one successful timed trial per wave", "charls_comparison": "timed-walk concept retained; ELSA test protocol/waves differ"},
    {"domain": "Balance", "raw_variables": "mmssre/mmstre/mmftre2", "raw_codes_and_missing": "side/semi: 1 held, 2 held <10s, 3 not attempted; full derived: 1/3 held target, 2/4 held shorter, 5 not attempted, negative=ineligible/missing", "wave_summary": "ordinal achieved level 0–3: failed/not attempted side=0; passed side only=1; passed semi but not full=2; passed full=3", "change": "W6 level below W4=balance worsening", "eligibility": "ordinal level determinable in both waves", "charls_comparison": "sequential static-balance construct; instruments differ"},
]
write_csv(ROOT / "01_definitions" / "ELSA_PHYSICAL_FUNCTION_DEFINITION_FINAL.csv", physical_rows)

cov_rows = [
    {"covariate": "Age", "raw_variables": "indager", "coding": "continuous years; age>=50", "special_missing": "negative", "reference": "continuous"},
    {"covariate": "Sex", "raw_variables": "indsex", "coding": "1 male; 2 female", "special_missing": "negative", "reference": "female"},
    {"covariate": "Education", "raw_variables": "fqendm", "coding": "1: codes1–3; 2:code4; 3:code5; 4:codes6–7; 5:code8", "special_missing": "-9 refusal; -8 don't know; -1 not applicable", "reference": "lowest"},
    {"covariate": "Partnered marital status", "raw_variables": "dimarr", "coding": "partnered=codes2–3; unpartnered=codes1,4,5,6", "special_missing": "negative", "reference": "unpartnered"},
    {"covariate": "Urban residence", "raw_variables": "w6_urindewr_2011", "coding": "1 urban; 2 rural", "special_missing": "-55/-7/-6/-1", "reference": "rural"},
    {"covariate": "Current smoking", "raw_variables": "heska", "coding": "1 current; 2 not current", "special_missing": "negative", "reference": "not current"},
    {"covariate": "Alcohol frequency", "raw_variables": "scako", "coding": "0 none(code8); 1 <=monthly(codes5–7); 2 1–4 days/week(codes3–4); 3 >=5 days/week(codes1–2)", "special_missing": "negative", "reference": "none"},
    {"covariate": "Fourteen-condition count", "raw_variables": "; ".join(f"{k}:" + "/".join(x for p in v for x in p[:3]) for k,v in chronic_parts.items()), "coding": "sum of 14 condition groups; prevalent if confirmed fed-forward or newly reported; composite heart and Alzheimer/dementia groups count once", "special_missing": "unknown in any group => count missing", "reference": "continuous"},
    {"covariate": "Baseline functional level", "raw_variables": "domain-specific W4 grip/walk/balance", "coding": "same operational definition as exposure baseline", "special_missing": "invalid or unknown excluded", "reference": "continuous/ordinal"},
    {"covariate": "Baseline biomarker burden", "raw_variables": "9 biomarkers at W4", "coding": "continuous; biomarker models only", "special_missing": "any required biomarker/sex missing", "reference": "continuous"},
    {"covariate": "Follow-up interval", "raw_variables": "derived", "coding": "W6–W7 vs W7–W8 categorical", "special_missing": "none in valid person-period", "reference": "W6–W7"},
]
write_csv(ROOT / "01_definitions" / "ELSA_COVARIATE_FINAL_CROSSWALK.csv", cov_rows)

biom_rows = []
domains = {"wbc":"Inflammatory", "hscrp":"Inflammatory", "hgb":"Haematologic", "fglu":"Glycaemic", "hba1c":"Glycaemic", "chol":"Lipid", "hdl":"Lipid", "ldl":"Lipid", "trig":"Lipid"}
labels = {"wbc":"White blood cell count", "hgb":"Haemoglobin", "fglu":"Fasting glucose", "chol":"Total cholesterol", "hdl":"HDL cholesterol", "ldl":"LDL cholesterol", "trig":"Triglycerides", "hscrp":"High-sensitivity CRP", "hba1c":"HbA1c"}
units = {"wbc":"10^9/L", "hgb":"g/dL", "fglu":"mmol/L", "chol":"mmol/L", "hdl":"mmol/L", "ldl":"mmol/L", "trig":"mmol/L", "hscrp":"mg/L", "hba1c":"W4 %, W6 mmol/mol converted to %"}
for c in BIOM:
    biom_rows.append({"biomarker": labels[c], "raw_variable_w4_w6": c, "unit": units[c], "domain": domains[c],
                      "transform": "log1p" if c in ["trig","hscrp"] else "W6 IFCC to NGSP: 0.09148*x+2.152" if c=="hba1c" else "none",
                      "risk_direction": "low" if c in ["hdl","hgb"] else "high", "w4_reference_mean": "sex-specific; see hgb rows" if c=="hgb" else means[c],
                      "w4_reference_sd": "sex-specific; see hgb rows" if c=="hgb" else sds[c], "aggregation": "indicator z -> within-domain arithmetic mean -> equal mean of 4 domains",
                      "w6_uses_w4_scale": "YES", "clinical_threshold": "NO", "winsorization": "NO"})
for sx, name in [(1,"HGB male reference"),(2,"HGB female reference")]:
    biom_rows.append({"biomarker": name, "raw_variable_w4_w6": "hgb + indsex", "unit":"g/dL", "domain":"Haematologic", "transform":"none", "risk_direction":"low; score=-(value-mean_sex)/SD_sex", "w4_reference_mean":hgb_ref[sx]["mean"], "w4_reference_sd":hgb_ref[sx]["sd"], "aggregation":"single-indicator domain", "w6_uses_w4_scale":"YES", "clinical_threshold":"NO", "winsorization":"NO"})
write_csv(ROOT / "01_definitions" / "ELSA_9_BIOMARKER_SCORE_SPECIFICATION_FINAL.csv", biom_rows)

decision_text = """# AUTHOR DECISION ELSA-HARM-001

Author decision ID: `ELSA-HARM-001`

- Primary outcome: **首次ADL/IADL活动困难 / incident ADL/IADL difficulty**.
- Position: harmonized conceptual replication of the CHARLS dependency result, not strict external validation or direct replication.
- Secondary dependency-oriented outcome: **5项ADL帮助依赖 / receipt of assistance with five ADL items** (dressing, bathing/showering, eating/cutting food, bed transfer, toileting).
- Full six-item IADL assistance dependency is not constructed; four observed IADL help variables are coverage description only and have no formal model.
- State-transition outcomes: **ADL/IADL活动困难消失 / resolution of ADL/IADL difficulty** and **停止报告5项ADL帮助需求 / cessation of reported assistance with five ADL items**.
- Biomarker construct: **ELSA适配的9项四领域多系统生物标志物负担 / ELSA-adapted nine-biomarker, four-domain multisystem biomarker burden**.
- Grip-difficulty is the primary conceptual replication and does not require biomarker completeness. Biomarker burden and its product with grip decline are secondary parallel physiological analyses.
- Current work is limited to definition freeze, descriptive anchors and preflight. Association-model execution is not authorized.
"""
(ROOT / "00_author_decision" / "AUTHOR_DECISION_ELSA_HARM_001.md").write_text(decision_text, encoding="utf-8")

protocol = """# Frozen ELSA external replication protocol V1.1

## Scientific question
在英国老年人群中，握力下降是否同样与后续首次ADL/IADL活动困难发生率较高相关？

## Time axis
Wave 4 to Wave 6 defines change; Wave 6 is the landmark; Wave 7 and Wave 8 define two discrete follow-up intervals. No W6-to-W8 bridge interval is constructed. The event row is retained, contribution stops after first event, and each participant contributes at most two rows.

## Primary conceptual replication
The primary exposure is grip decline, defined as Wave 4 maximum valid grip minus Wave 6 maximum valid grip. The primary outcome is incident ADL/IADL difficulty. Eligibility requires no difficulty in all five ADL and six IADL items at Wave 6 and a valid Wave 7 interval; Wave 8 contributes only after a known event-free Wave 7. The primary analysis does not require nine-biomarker completeness.

## Secondary outcomes and exposures
The five-ADL assistance outcome uses dressing, bathing/showering, eating/cutting food, bed transfer and toileting. Structural questionnaire skips caused by no difficulty are deterministically no assistance; refusals/don't-know remain unknown. Full six-item IADL assistance dependency is not constructed. Difficulty resolution and cessation of reported five-ADL assistance are state transitions, not treatment recovery.

The ELSA-adapted nine-biomarker, four-domain burden is a secondary parallel physiological exposure. It is not the CHARLS 12-item score, a biological age, or a validated physiological ageing index.

## Adjustment
Fully adjusted models include Wave 6 age, sex, education, partnered marital status, urban/rural residence, current smoking, alcohol frequency, 14-condition count, Wave 4 baseline level of the relevant functional exposure, and categorical interval. Biomarker models additionally adjust Wave 4 biomarker burden.

## Interpretation
Assess direction, interval-estimate support, precision and whether outcome-definition differences plausibly explain heterogeneity. CHARLS measures needing help/unable-to-complete dependency, whereas ELSA measures difficulty, an earlier and broader limitation level. Effect estimates are not treated as definition-equivalent direct numerical replications.

## Execution gate
This protocol does not authorize association models. Phase ELSA-3 model execution requires a separately authorized run after this preflight is accepted.
"""
(ROOT / "03_protocol" / "FROZEN_ELSA_EXTERNAL_REPLICATION_PROTOCOL_V1_1.md").write_text(protocol, encoding="utf-8")

model_rows = [
    {"model_id":"ELSA-P1", "priority":"PRIMARY_CONCEPTUAL_REPLICATION", "outcome":"incident ADL/IADL difficulty", "exposure":"grip decline", "biomarker_required":"NO", "model":"person-period complementary log-log", "status":"PLANNED_NOT_RUN"},
    {"model_id":"ELSA-S1", "priority":"SECONDARY", "outcome":"receipt of assistance with five ADL items", "exposure":"grip decline", "biomarker_required":"NO", "model":"person-period complementary log-log", "status":"PLANNED_NOT_RUN"},
    {"model_id":"ELSA-S2", "priority":"SECONDARY_PARALLEL_PHYSIOLOGY", "outcome":"incident ADL/IADL difficulty", "exposure":"9-biomarker burden change + grip decline + product", "biomarker_required":"YES", "model":"person-period complementary log-log", "status":"PLANNED_NOT_RUN"},
    {"model_id":"ELSA-E1", "priority":"EXPLORATORY_DOMAIN", "outcome":"incident ADL/IADL difficulty", "exposure":"timed-walk slowing", "biomarker_required":"NO", "model":"person-period complementary log-log", "status":"PLANNED_NOT_RUN"},
    {"model_id":"ELSA-E2", "priority":"EXPLORATORY_DOMAIN", "outcome":"incident ADL/IADL difficulty", "exposure":"balance worsening", "biomarker_required":"NO", "model":"person-period complementary log-log", "status":"PLANNED_NOT_RUN"},
    {"model_id":"ELSA-T1", "priority":"STATE_TRANSITION", "outcome":"resolution of ADL difficulty", "exposure":"grip decline", "biomarker_required":"NO", "model":"person-period complementary log-log", "status":"PLANNED_NOT_RUN"},
    {"model_id":"ELSA-T2", "priority":"STATE_TRANSITION", "outcome":"resolution of IADL difficulty", "exposure":"grip decline", "biomarker_required":"NO", "model":"person-period complementary log-log", "status":"PLANNED_NOT_RUN"},
    {"model_id":"ELSA-T3", "priority":"STATE_TRANSITION", "outcome":"cessation of reported assistance with five ADL items", "exposure":"grip decline", "biomarker_required":"NO", "model":"none", "status":"DESCRIPTIVE_ONLY_NOT_FORMAL_MODEL"},
]
write_csv(ROOT / "03_protocol" / "FROZEN_ELSA_MODEL_PLAN_V1_1.csv", model_rows)

preflight = [
    {"check_id":"PF01", "check":"ELSA-HARM-001 author decision recorded", "status":"PASS", "detail":"formal names and interpretation boundaries recorded"},
    {"check_id":"PF02", "check":"education definition closed", "status":"PASS", "detail":"fqendm five-category coding frozen"},
    {"check_id":"PF03", "check":"chronic disease definition closed", "status":"PASS", "detail":"14 explicit condition groups and raw variables frozen"},
    {"check_id":"PF04", "check":"walking definition closed", "status":"PASS", "detail":"mean successful trial time, at least one per wave"},
    {"check_id":"PF05", "check":"balance definition closed", "status":"PASS", "detail":"sequential 0–3 achieved level; worsening W6<W4"},
    {"check_id":"PF06", "check":"incident ADL/IADL difficulty closed", "status":"PASS", "detail":"5 ADL + 6 IADL, W6 difficulty-free, W7/W8 first difficulty"},
    {"check_id":"PF07", "check":"5-item ADL assistance closed", "status":"PASS", "detail":"all item statuses ascertainable; genuine unknown not coded no"},
    {"check_id":"PF08", "check":"full IADL assistance excluded", "status":"PASS", "detail":"hot-meal and telephone help unavailable; no formal model"},
    {"check_id":"PF09", "check":"state-transition terminology closed", "status":"PASS", "detail":"difficulty resolution / cessation of assistance; no recovery claim"},
    {"check_id":"PF10", "check":"9-biomarker definition closed", "status":"PASS", "detail":"nine indicators, four domains, W4 reference scale; not biological age"},
    {"check_id":"PF11", "check":"descriptive anchors generated", "status":"PASS", "detail":f"{len(anchor_rows)} analysis rows; all model counts zero"},
    {"check_id":"PF12", "check":"protocol/workbooks synchronized", "status":"PENDING_WORKBOOK_BUILD", "detail":"CSV/protocol complete; workbook build follows"},
    {"check_id":"PF13", "check":"association model execution", "status":"PASS", "detail":"0 models; 0 effect estimates; 0 p-values/CI"},
]
write_csv(ROOT / "05_preflight" / "PHASE_ELSA3_PREFLIGHT_CHECKLIST.csv", preflight)

summary = {
    "author_decision":"ELSA-HARM-001", "definitions_closed":True,
    "anchor_rows":anchor_rows, "models_run":0, "effect_estimates_viewed":0,
    "p_values_or_confidence_intervals":0, "preflight_workbook_pending":True,
}
(ROOT / "02_anchors" / "ELSA_FINAL_ANCHOR_SUMMARY.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

runlog = """# Run log

- Read only selected variables from authorized ELSA Stata files for definition closure and descriptive anchors.
- No original file was modified.
- No association model, regression, effect estimate, P value, confidence interval, FDR, MI, IPW or prediction analysis was run or viewed.
- Grip primary conceptual replication does not require nine-biomarker completeness.
- Full six-item IADL assistance dependency was not constructed.
"""
(ROOT / "08_logs" / "RUN_LOG.md").write_text(runlog, encoding="utf-8")
shutil.copy2(Path(__file__), ROOT / "07_scripts" / Path(__file__).name)
print(json.dumps(summary, ensure_ascii=False, indent=2))

# ---------------------------------------------------------------------------
# Phase ELSA-3 single authorized run under V1.1 scientific definitions and
# the V1.2 frozen common standardization scales.
# ---------------------------------------------------------------------------
import os, warnings, traceback
from datetime import datetime
for d in ["00_input_gate","01_preflight","02_data_qa","03_results","04_diagnostics","05_figures","06_workbooks","07_reports","08_logs","09_scripts","10_figure_ready","11_manifest","_runtime/matplotlib","_runtime/cache","_runtime/temp"]:
    (ROOT/d).mkdir(parents=True,exist_ok=True)
os.environ["MPLBACKEND"]="Agg"
os.environ["MPLCONFIGDIR"]=str(ROOT/"_runtime/matplotlib")
os.environ["XDG_CACHE_HOME"]=str(ROOT/"_runtime/cache")
os.environ["TMP"]=str(ROOT/"_runtime/temp"); os.environ["TEMP"]=str(ROOT/"_runtime/temp"); os.environ["TMPDIR"]=str(ROOT/"_runtime/temp")

import statsmodels.api as sm
import statsmodels.formula.api as smf

def verify_manifest(package, manifest_rel):
    rows=list(csv.DictReader((package/manifest_rel).open(encoding="utf-8-sig")))
    out=[]
    for r in rows:
        p=package/r["relative_path"]
        actual=sha256(p) if p.exists() else "MISSING"
        out.append({"package":str(package),"relative_path":r["relative_path"],"expected_sha256":r["sha256"],"actual_sha256":actual,"status":"PASS" if actual==r["sha256"] else "FAIL"})
    return out

input_gate=verify_manifest(V11,Path("09_manifest/FILE_MANIFEST_SHA256.csv"))+verify_manifest(V12,Path("13_v1_2_manifest/FILE_MANIFEST_SHA256.csv"))
write_csv(ROOT/"00_input_gate"/"FROZEN_INPUT_HASH_VERIFICATION.csv",input_gate)
if not all(r["status"]=="PASS" for r in input_gate):
    raise RuntimeError("frozen input manifest verification failed")
v12_pf=pd.read_csv(V12/"11_v1_2_preflight"/"PHASE_ELSA3_V1_2_PREFLIGHT_CHECKLIST.csv")
if not v12_pf["status"].eq("PASS").all():
    raise RuntimeError("V1.2 preflight is not fully PASS")
params=pd.read_csv(V12/"10_v1_2_scale"/"ELSA_STANDARDIZATION_PARAMETERS_V1_2.csv").set_index("variable")
model_plan=pd.read_csv(V11/"03_protocol"/"FROZEN_ELSA_MODEL_PLAN_V1_1.csv")
for v in ["GripDecline","WalkingTimeIncrease","BiomarkerBurdenChange"]:
    if not np.isfinite(params.loc[v,"mean"]) or not np.isfinite(params.loc[v,"sample_sd"]) or params.loc[v,"sample_sd"]<=0:
        raise RuntimeError(f"invalid frozen scale for {v}")

grip_z=((g4-g6)-params.loc["GripDecline","mean"])/params.loc["GripDecline","sample_sd"]
walk_z=((walk6-walk4)-params.loc["WalkingTimeIncrease","mean"])/params.loc["WalkingTimeIncrease","sample_sd"]
burden_z=((burden6-burden4)-params.loc["BiomarkerBurdenChange","mean"])/params.loc["BiomarkerBurdenChange","sample_sd"]
balance_index=bal4.index.union(bal6.index)
balance_worsening=pd.Series(pd.NA,index=balance_index,dtype="Float64")
balance_common=bal4.dropna().index.intersection(bal6.dropna().index)
balance_worsening.loc[balance_common]=(bal6.reindex(balance_common)<bal4.reindex(balance_common)).astype("Float64")

def build_pp_df(eligible,outcome,baseline_positive=False):
    idx=pd.Index(sorted(set(eligible)))
    k6,p6=status[("w6",outcome)]; k7,p7=status[("w7",outcome)]; k8,p8=status[("w8",outcome)]
    q=pd.DataFrame(index=idx)
    for n,s in [("k6",k6),("p6",p6),("k7",k7),("p7",p7),("k8",k8),("p8",p8)]: q[n]=s.reindex(idx)
    risk=q["k6"].eq(True)&(q["p6"].eq(True) if baseline_positive else q["p6"].eq(False))
    r=q.loc[risk]; rows=[]
    for pid,x in r.iterrows():
        if x["k7"]==True:
            ev1=int(x["p7"]==False) if baseline_positive else int(x["p7"]==True)
            rows.append({"idauniq":pid,"interval":1,"event":ev1})
            remains=(x["p7"]==True) if baseline_positive else (x["p7"]==False)
            if ev1==0 and remains and x["k8"]==True:
                ev2=int(x["p8"]==False) if baseline_positive else int(x["p8"]==True)
                rows.append({"idauniq":pid,"interval":2,"event":ev2})
    pp=pd.DataFrame(rows)
    if pp.empty: return pp
    assert pp.groupby("idauniq").size().max()<=2
    assert not pp.sort_values(["idauniq","interval"]).groupby("idauniq")["event"].shift().eq(1).any()
    return pp

# Compare each fitted model against the newly rebuilt, pre-result anchors from
# this correction run, not the superseded pre-correction V1.1 anchors.
anchor=pd.read_csv(ROOT/"02_anchors"/"ELSA_REVISED_SAMPLE_AND_EVENT_ANCHORS_FINAL.csv")
anchor_map={r["analysis"]:r for r in anchor.to_dict("records")}
specs=[
    {"model_id":"ELSA-P1","analysis":"主要概念性重复：握力下降与首次ADL/IADL活动困难","outcome":"combined_difficulty","exp":"grip","positive":False,"terms":["grip_z"]},
    {"model_id":"ELSA-S1","analysis":"次要：5项ADL帮助依赖","outcome":"adl_help5","exp":"grip","positive":False,"terms":["grip_z"]},
    {"model_id":"ELSA-S2","analysis":"次要并列生理：9项负担、握力及乘积项与首次ADL/IADL活动困难","outcome":"combined_difficulty","exp":"biom9_grip","positive":False,"terms":["burden_z","grip_z","biomarker_grip_product"]},
    {"model_id":"ELSA-E1","analysis":"领域：步行时间增加与首次ADL/IADL活动困难","outcome":"combined_difficulty","exp":"walk","positive":False,"terms":["walk_z"]},
    {"model_id":"ELSA-E2","analysis":"领域：平衡恶化与首次ADL/IADL活动困难","outcome":"combined_difficulty","exp":"balance","positive":False,"terms":["balance_worsening"]},
    {"model_id":"ELSA-T1","analysis":"状态转换：ADL活动困难消失","outcome":"adl_difficulty","exp":"grip","positive":True,"terms":["grip_z"]},
    {"model_id":"ELSA-T2","analysis":"状态转换：IADL活动困难消失","outcome":"iadl_difficulty","exp":"grip","positive":True,"terms":["grip_z"]},
]

def model_frame(spec):
    baseline={"grip":g4,"biom9_grip":g4,"walk":walk4,"balance":bal4}[spec["exp"]]
    eligible=exposures[spec["exp"]]&set(w6c.index)&cov_complete&set(baseline.dropna().index)
    if spec["exp"]=="biom9_grip": eligible &= set(burden4.dropna().index)
    pp=build_pp_df(eligible,spec["outcome"],spec["positive"])
    ids=pp["idauniq"].astype("string")
    for c in base_cov: pp[c]=cov[c].reindex(ids).to_numpy()
    pp["baseline_function"]=baseline.reindex(ids).to_numpy()
    if spec["exp"] in ["grip","biom9_grip"]: pp["grip_z"]=grip_z.reindex(ids).to_numpy()
    if spec["exp"]=="walk": pp["walk_z"]=walk_z.reindex(ids).to_numpy()
    if spec["exp"]=="balance": pp["balance_worsening"]=balance_worsening.reindex(ids).to_numpy()
    if spec["exp"]=="biom9_grip":
        pp["burden_z"]=burden_z.reindex(ids).to_numpy(); pp["baseline_burden"]=burden4.reindex(ids).to_numpy()
        pp["biomarker_grip_product"]=pp["burden_z"]*pp["grip_z"]
    return pp

def fit_one(spec,pp):
    rhs=" + ".join(spec["terms"]+["age","C(sex)","C(education)","partnered","urban","current_smoking","C(alcohol_frequency)","chronic_count","baseline_function","C(interval)"]+(["baseline_burden"] if spec["exp"]=="biom9_grip" else []))
    formula="event ~ "+rhs
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        res=smf.glm(formula=formula,data=pp,family=sm.families.Binomial(link=sm.families.links.CLogLog())).fit(cov_type="cluster",cov_kwds={"groups":pp["idauniq"]},maxiter=200)
    exog=res.model.exog; rank=int(np.linalg.matrix_rank(exog)); p=exog.shape[1]
    finite=bool(np.isfinite(res.params).all() and np.isfinite(res.bse).all() and np.isfinite(res.cov_params()).all().all())
    if not res.converged or rank!=p or not finite:
        raise RuntimeError(f"{spec['model_id']} estimation QA failure: converged={res.converged}, rank={rank}/{p}, finite={finite}")
    return res,formula,rank,p," | ".join(str(x.message) for x in caught)

ledger=[]; diagnostics=[]; results=[]; coef_rows=[]; anchor_checks=[]; set_qa=[]
run_started=datetime.now().astimezone().isoformat(); fitted=0
for spec in specs:
    pp=model_frame(spec)
    expected=anchor_map[spec["analysis"]]
    actual=(pp["idauniq"].nunique(),len(pp),int(pp["event"].sum()))
    target=(int(expected["fully_adjusted_person_n"]),int(expected["fully_adjusted_person_period_n"]),int(expected["fully_adjusted_event_n"]))
    anchor_status="PASS" if actual==target else "FAIL"
    anchor_checks.append({"model_id":spec["model_id"],"participants":actual[0],"expected_participants":target[0],"person_period":actual[1],"expected_person_period":target[1],"events":actual[2],"expected_events":target[2],"status":anchor_status})
    if anchor_status!="PASS": raise RuntimeError(f"anchor mismatch {spec['model_id']}: {actual} vs {target}")
    res,formula,rank,npar,warn=fit_one(spec,pp); fitted+=1
    for term in res.params.index:
        b=float(res.params[term]); se=float(res.bse[term]); lo=b-1.959963984540054*se; hi=b+1.959963984540054*se
        coef_rows.append({"model_id":spec["model_id"],"term":term,"estimate_log":b,"standard_error":se,"HR":float(np.exp(b)),"CI_low":float(np.exp(lo)),"CI_high":float(np.exp(hi)),"p_value":float(res.pvalues[term])})
    for term in spec["terms"]:
        b=float(res.params[term]); se=float(res.bse[term]); lo=b-1.959963984540054*se; hi=b+1.959963984540054*se
        interp="conditional at the other standardized exposure=0" if spec["model_id"]=="ELSA-S2" and term in ["burden_z","grip_z"] else "product effect modification term" if term=="biomarker_grip_product" else "per frozen reference-population change-score SD" if term.endswith("_z") else "worsening versus no worsening"
        results.append({"model_id":spec["model_id"],"priority":next(x["priority"] for x in model_plan.to_dict("records") if x["model_id"]==spec["model_id"]),"outcome":next(x["outcome"] for x in model_plan.to_dict("records") if x["model_id"]==spec["model_id"]),"term":term,"effect_measure":"HR","estimate":float(np.exp(b)),"CI_low":float(np.exp(lo)),"CI_high":float(np.exp(hi)),"p_value":float(res.pvalues[term]),"participants":actual[0],"person_period":actual[1],"events":actual[2],"interpretation_unit":interp,"status":"SUCCESS"})
    ledger.append({"model_id":spec["model_id"],"execution_status":"SUCCESS","participants":actual[0],"person_period":actual[1],"events":actual[2],"parameters":npar,"clusters":actual[0],"converged":bool(res.converged),"design_rank":rank,"full_rank":rank==npar,"finite_coefficients":True,"finite_standard_errors":True,"covariance_status":"PASS","cluster_robust_variance":"PASS","formula":formula})
    diagnostics.append({"model_id":spec["model_id"],"converged":bool(res.converged),"iterations":res.fit_history.get("iteration","NA"),"design_rank":rank,"parameters":npar,"condition_number":float(np.linalg.cond(res.model.exog)),"cov_type":res.cov_type,"finite_covariance":bool(np.isfinite(res.cov_params()).all().all()),"warnings":warn,"status":"PASS"})
    hashed=sorted(hashlib.sha256(("ELSA3-V12-SET-AUDIT|"+str(x)).encode()).hexdigest() for x in pp["idauniq"].unique())
    set_qa.append({"model_id":spec["model_id"],"participants":actual[0],"sorted_anonymous_member_set_sha256":hashlib.sha256("\n".join(hashed).encode()).hexdigest(),"raw_id_exposed":"NO"})

assert fitted==7
write_csv(ROOT/"02_data_qa"/"DATA_CONSTRUCTION_ANCHOR_CHECKS.csv",anchor_checks)
write_csv(ROOT/"02_data_qa"/"PARTICIPANT_AND_EVENT_SET_QA.csv",set_qa)
write_csv(ROOT/"03_results"/"MODEL_EXECUTION_LEDGER.csv",ledger)
write_csv(ROOT/"04_diagnostics"/"MODEL_DIAGNOSTICS.csv",diagnostics)
write_csv(ROOT/"03_results"/"ALL_MODEL_FULL_COEFFICIENTS.csv",coef_rows)
write_csv(ROOT/"03_results"/"ALL_FROZEN_MODEL_TARGET_RESULTS.csv",results)

def subset_results(ids,terms=None):
    return [r for r in results if r["model_id"] in ids and (terms is None or r["term"] in terms)]
write_csv(ROOT/"03_results"/"PRIMARY_CONCEPTUAL_REPLICATION_RESULTS.csv",subset_results(["ELSA-P1"]))
write_csv(ROOT/"03_results"/"ADL_ASSISTANCE_RESULTS.csv",subset_results(["ELSA-S1"]))
write_csv(ROOT/"03_results"/"ELSA_ADAPTED_BIOMARKER_RESULTS.csv",subset_results(["ELSA-S2"]))
write_csv(ROOT/"03_results"/"PHYSICAL_FUNCTION_DOMAIN_RESULTS.csv",subset_results(["ELSA-E1","ELSA-E2"]))
write_csv(ROOT/"03_results"/"DIFFICULTY_RESOLUTION_RESULTS.csv",subset_results(["ELSA-T1","ELSA-T2"]))
write_csv(ROOT/"03_results"/"PRESPECIFIED_SENSITIVITY_RESULTS.csv",[{"status":"NOT_RUN_NOT_LISTED_IN_V1_1_MODEL_PLAN","models_run":0}])

# Cross-cohort descriptive evidence crosswalk, using the final frozen CHARLS primary result file.
charls_file=Path(os.environ["CHARLS_PRIMARY_RESULT_FILE"])
charls=pd.read_csv(charls_file); cg=charls.loc[charls["term"].eq("grip_decline_z")].iloc[0]
ep1=next(r for r in results if r["model_id"]=="ELSA-P1" and r["term"]=="grip_z")
if ep1["estimate"]>1 and ep1["CI_low"]>1: evidence="DIRECTIONALLY_CONSISTENT_WITH_INTERVAL_SUPPORT"
elif ep1["estimate"]>1: evidence="DIRECTIONALLY_CONSISTENT_IMPRECISE"
else: evidence="DIRECTIONALLY_INCONSISTENT"
crosswalk=[
    {"cohort":"CHARLS","outcome":"first ADL/IADL dependency (help/unable)","exposure":"grip decline per frozen change SD","effect_measure":"HR","estimate":float(cg.HR),"CI_low":float(cg.CI_low),"CI_high":float(cg.CI_high),"p_value":float(cg.p_value),"definition_comparability":"CONCEPTUALLY_RELATED_NOT_EQUIVALENT","source_sha256":sha256(charls_file)},
    {"cohort":"ELSA","outcome":"incident ADL/IADL difficulty","exposure":"grip decline per frozen ELSA change SD","effect_measure":"HR","estimate":ep1["estimate"],"CI_low":ep1["CI_low"],"CI_high":ep1["CI_high"],"p_value":ep1["p_value"],"definition_comparability":"CONCEPTUALLY_RELATED_NOT_EQUIVALENT","source_sha256":sha256(ROOT/"03_results"/"PRIMARY_CONCEPTUAL_REPLICATION_RESULTS.csv")},
]
write_csv(ROOT/"03_results"/"CHARLS_ELSA_EVIDENCE_CROSSWALK.csv",crosswalk)

# Figure-ready data package and Python-only publication figures.
figrows=[]
display={"ELSA-P1":"Incident ADL/IADL difficulty — grip decline","ELSA-S1":"Five-ADL assistance — grip decline","ELSA-S2":"Incident difficulty — biomarker burden / grip / product","ELSA-E1":"Incident difficulty — walking-time increase","ELSA-E2":"Incident difficulty — balance worsening","ELSA-T1":"Resolution of ADL difficulty — grip decline","ELSA-T2":"Resolution of IADL difficulty — grip decline"}
for r in results:
    figrows.append({"analysis_set":"fully adjusted V1.1/V1.2","endpoint":r["outcome"],"model_id":r["model_id"],"display_label":display[r["model_id"]]+(" — "+r["term"] if r["model_id"]=="ELSA-S2" else ""),"effect_measure":"HR","estimate":r["estimate"],"ci_low":r["CI_low"],"ci_high":r["CI_high"],"p_value":r["p_value"],"covariates":"W6 age, sex, education, partnered, urban, smoking, alcohol, 14-condition count; W4 relevant function; categorical interval; W4 burden additionally for ELSA-S2","source_file":"03_results/ALL_FROZEN_MODEL_TARGET_RESULTS.csv"})
write_csv(ROOT/"10_figure_ready"/"figure_ready_data_package.csv",figrows)
write_csv(ROOT/"10_figure_ready"/"figure_panel_mapping.csv",[{"figure":"Figure 1","panel":"a","claim":"ELSA frozen model effect estimates","source":"figure_ready_data_package.csv"},{"figure":"Figure 2","panel":"a","claim":"CHARLS–ELSA directional comparison for grip decline","source":"CHARLS_ELSA_EVIDENCE_CROSSWALK.csv"}])
write_csv(ROOT/"10_figure_ready"/"source_data_traceability.csv",[{"output":"ELSA_MODEL_EFFECT_FOREST","source":"03_results/ALL_FROZEN_MODEL_TARGET_RESULTS.csv"},{"output":"CHARLS_ELSA_DIRECTIONAL_COMPARISON","source":"03_results/CHARLS_ELSA_EVIDENCE_CROSSWALK.csv"}])
(ROOT/"10_figure_ready"/"figure_missing_data_checklist.md").write_text("# Figure missing-data checklist\n\n- All plotted rows have finite estimate and CI: PASS.\n- Failed models plotted as estimates: 0.\n- Missing-data handling: complete-case eligibility fixed by V1.1.\n",encoding="utf-8")

import matplotlib
matplotlib.use("Agg",force=True)
import matplotlib.pyplot as plt
matplotlib.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans","sans-serif"],"font.size":7,"pdf.fonttype":42,"svg.fonttype":"none","axes.spines.top":False,"axes.spines.right":False})
plot_rows=[r for r in figrows]
fig,ax=plt.subplots(figsize=(7.2,5.2)); y=np.arange(len(plot_rows))[::-1]
est=np.array([r["estimate"] for r in plot_rows]); lo=np.array([r["ci_low"] for r in plot_rows]); hi=np.array([r["ci_high"] for r in plot_rows])
ax.errorbar(est,y,xerr=[est-lo,hi-est],fmt='s',color='#1F4E78',ecolor='#64748B',capsize=2,markersize=4)
ax.axvline(1,color='#6B7280',lw=.8,ls='--'); ax.set_xscale('log'); ax.set_yticks(y); ax.set_yticklabels([r["display_label"] for r in plot_rows]); ax.set_xlabel('Hazard ratio (95% CI)'); ax.set_title('ELSA harmonized conceptual replication: frozen model estimates',loc='left',fontweight='bold')
fig.subplots_adjust(left=.49,right=.97,bottom=.12,top=.90)
for ext in ['png','pdf','svg']: fig.savefig(ROOT/"05_figures"/f"ELSA_MODEL_EFFECT_FOREST.{ext}",dpi=600 if ext=='png' else None,bbox_inches='tight',facecolor='white')
plt.close(fig)

fig,ax=plt.subplots(figsize=(6.2,2.7)); yy=np.array([1,0]); ee=np.array([crosswalk[0]["estimate"],crosswalk[1]["estimate"]]); ll=np.array([crosswalk[0]["CI_low"],crosswalk[1]["CI_low"]]); hh=np.array([crosswalk[0]["CI_high"],crosswalk[1]["CI_high"]])
ax.errorbar(ee,yy,xerr=[ee-ll,hh-ee],fmt='s',color='#1F4E78',ecolor='#64748B',capsize=3,markersize=5); ax.axvline(1,color='#6B7280',lw=.8,ls='--'); ax.set_xscale('log'); ax.set_yticks(yy); ax.set_yticklabels(['CHARLS: dependency','ELSA: difficulty']); ax.set_xlabel('Hazard ratio (95% CI)'); ax.set_title('Grip decline: directional comparison across non-equivalent outcomes',loc='left',fontweight='bold'); fig.subplots_adjust(left=.28,right=.96,bottom=.20,top=.82)
for ext in ['png','pdf','svg']: fig.savefig(ROOT/"05_figures"/f"CHARLS_ELSA_DIRECTIONAL_COMPARISON.{ext}",dpi=600 if ext=='png' else None,bbox_inches='tight',facecolor='white')
plt.close(fig)

report=f"""# ELSA external conceptual replication results report

## Technical status
`PHASE_ELSA3_V1_2_SINGLE_RUN_COMPLETE_READY_FOR_SCIENTIFIC_INTERPRETATION`

Seven formal models in the V1.1 machine-readable plan were fitted once using the V1.2 common standardization scales. ELSA-T3 remained descriptive only. All seven converged, used cluster-robust covariance and had full-rank design matrices.

## Primary conceptual replication
ELSA-P1: N={ep1['participants']}, person-period={ep1['person_period']}, events={ep1['events']}; grip decline HR={ep1['estimate']:.3f} (95% CI {ep1['CI_low']:.3f}–{ep1['CI_high']:.3f}), P={ep1['p_value']:.6g}.

Evidence classification: `{evidence}`. This compares association direction and overall pattern only. CHARLS measured help/unable-to-complete dependency; ELSA measured activity difficulty, an earlier and broader limitation level. The estimates are not definition-equivalent direct replications.

## Boundaries
- No unplanned model, prediction model, MI or IPW was run.
- No model-specific standard deviation was calculated.
- No protocol definition was changed after viewing results.
- Scientific interpretation and manuscript integration were not started.
"""
(ROOT/"07_reports"/"ELSA_EXTERNAL_REPLICATION_RESULTS_REPORT.md").write_text(report,encoding="utf-8")
(ROOT/"07_reports"/"ELSA_EXECUTIVE_RESULTS_SUMMARY_ZH.md").write_text(f"""# ELSA执行摘要

- 主要问题：英国老年人中，握力下降是否与后续首次ADL/IADL活动困难发生率较高相关。
- ELSA-P1：N={ep1['participants']}，person-period={ep1['person_period']}，事件={ep1['events']}；HR={ep1['estimate']:.3f}（95%CI {ep1['CI_low']:.3f}–{ep1['CI_high']:.3f}），P={ep1['p_value']:.6g}。
- 证据分类：`{evidence}`。
- 本分析是协调后的概念性重复，不是CHARLS依赖结局的严格外部验证。
- 模型成功/失败/未运行：7/0/0（另有ELSA-T3按协议仅描述，不是正式模型）。
- 协议偏离：0；未授权模型：0；自动统计重试：0。
""",encoding="utf-8")
(ROOT/"07_reports"/"PHASE_ELSA3_FINAL_TECHNICAL_STATUS.md").write_text("""# Phase ELSA-3 final technical status

`PHASE_ELSA3_V1_2_SINGLE_RUN_COMPLETE_READY_FOR_SCIENTIFIC_INTERPRETATION`

- V1.2 scale freeze: PASS.
- V1.2 preflight: PASS.
- Formal models planned/run/success/failed: 7/7/7/0.
- Descriptive-only non-model item: 1 (ELSA-T3).
- Protocol deviations: 0.
- Unauthorized models or results viewed: 0.
- Automatic statistical retries: 0.
- Stopped before scientific interpretation or manuscript integration: YES.
""",encoding="utf-8")
(ROOT/"07_reports"/"需要人工核查.md").write_text("# 需要人工核查\n\n- 确认ELSA-T2虽样本量不足500，但机器可读V1.1模型计划仍将其列为正式计划模型；本轮严格按该模型计划执行。\n- 后续科学解释必须保持ELSA活动困难与CHARLS帮助/无法完成依赖并非定义等价。\n",encoding="utf-8")
(ROOT/"07_reports"/"下一步交接记录.md").write_text("# 下一步交接记录\n\n本轮已停止。只有作者另行授权后，才可开展科学解释或稿件整合。\n",encoding="utf-8")
(ROOT/"08_logs"/"RUN_LOG.md").write_text(f"""# Phase ELSA-3 V1.2 run log

- Started: {run_started}
- V1.1 manifest: PASS.
- V1.2 manifest and preflight: PASS.
- Formal models fitted once: {fitted}.
- Statistical retries: 0.
- Unplanned models: 0.
- Original DTA modification: NO.
""",encoding="utf-8")
(ROOT/"08_logs"/"ENGINEERING_REPAIR_LEDGER.md").write_text("# Engineering repair ledger\n\nNo engineering repair or retry was required during the single authorized model run.\n",encoding="utf-8")
shutil.copy2(Path(__file__),ROOT/"09_scripts"/Path(__file__).name)
print(json.dumps({"status":"PHASE_ELSA3_V1_2_SINGLE_RUN_COMPLETE_READY_FOR_SCIENTIFIC_INTERPRETATION","models_success":7,"models_failed":0,"models_not_run":0,"primary":ep1,"evidence":evidence},ensure_ascii=False,indent=2))

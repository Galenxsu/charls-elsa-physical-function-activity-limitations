from pathlib import Path
import os

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
OUT = Path(
    os.environ.get(
        "SENSITIVITY_VERIFY_ROOT",
        Path(os.environ.get("SENSITIVITY_OUTPUT_ROOT", REPO_ROOT / "outputs" / "sensitivity")),
    )
)
r = pd.read_csv(OUT / "SENSITIVITY_MODEL_RESULTS.csv")
s = pd.read_csv(OUT / "CHARLS_2015_SKIP_PATH_AUDIT.csv")

checks = []


def check(name, ok, detail):
    checks.append({"check": name, "status": "PASS" if ok else "FAIL", "detail": detail})


check("expected result row count", len(r) == 15, len(r))
check("all estimates finite", r[["HR", "CI_low", "CI_high", "p_value"]].notna().all().all(), "finite fields")
check("all interval ordering valid", ((r.CI_low < r.HR) & (r.HR < r.CI_high)).all(), "CI_low < HR < CI_high")
check("CHARLS anchor HR", abs(r.query("cohort == 'CHARLS' and analysis == 'primary reproduction'").HR.iloc[0] - 1.2137716269311625) < 1e-12, "exact reproduction")
check("ELSA anchor HR", abs(r.query("cohort == 'ELSA' and analysis == 'primary reproduction'").HR.iloc[0] - 1.1152778685466338) < 1e-12, "exact reproduction")
check("skip categories reconcile", int(s.N.iloc[1:].sum()) == int(s.N.iloc[0]) == 3307, f"{int(s.N.iloc[1:].sum())}/3307")
check("only one skip category eligible", int((s.eligible_for_recode == 'YES').sum()) == 1, int((s.eligible_for_recode == 'YES').sum()))

for cohort in ["CHARLS", "ELSA"]:
    for analysis in ["earlier ADL/IADL status", "BMI adjustment"]:
        z = r[(r.cohort == cohort) & (r.analysis == analysis)]
        same = z.participants.nunique() == z.person_periods.nunique() == z.events.nunique() == 1
        check(f"{cohort} {analysis} same-sample comparison", same, f"N={z.participants.tolist()}, rows={z.person_periods.tolist()}, events={z.events.tolist()}")

check("all sensitivity directions positive", (r.HR > 1).all(), f"min HR={r.HR.min():.6f}")

pd.DataFrame(checks).to_csv(OUT / "INDEPENDENT_OUTPUT_VERIFICATION.csv", index=False, encoding="utf-8-sig")
if any(x["status"] == "FAIL" for x in checks):
    raise SystemExit("Sensitivity-output verification failed")
print(pd.DataFrame(checks).to_string(index=False))

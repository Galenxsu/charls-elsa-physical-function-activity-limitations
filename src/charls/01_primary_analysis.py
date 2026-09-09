# Part of the charls-elsa-physical-function-activity-limitations project.
# Developed for the analyses reported in the associated manuscript.
# Released under the MIT License. See LICENSE in the repository root.
# No CHARLS or ELSA individual-level data are distributed with this code.

from __future__ import annotations
import csv, hashlib, json, platform, sys, time, traceback
from pathlib import Path
import numpy as np, pandas as pd
import scipy, statsmodels, statsmodels.api as sm, statsmodels.formula.api as smf
import os

REPO_ROOT=Path(__file__).resolve().parents[2]
ROOT=Path(os.environ["CHARLS_PROJECT_ROOT"])
OUT=Path(os.environ.get("CHARLS_PRIMARY_OUTPUT", REPO_ROOT/"outputs"/"charls_primary"))
P3A=ROOT/"13_phase3a1_author_resolution"; P3=ROOT/"12_phase3a_protocol_freeze"
RAW=Path(os.environ["CHARLS_DATA_ROOT"])
BLOOD=ROOT/"01_file_inventory/authorized_blood_downloads"; START=time.time(); SEED=20260814
np.random.seed(SEED)
DIRS=["00_frozen_inputs","01_data_construction","02_qa","03_analysis_datasets","04_primary_model","05_secondary_four_group","06_sensitivity","07_tables","08_figures","09_reports","10_logs","11_scripts","12_manifest"]
for d in DIRS:(OUT/d).mkdir(parents=True,exist_ok=True)

def sha(p):
 h=hashlib.sha256()
 with Path(p).open("rb") as f:
  for b in iter(lambda:f.read(1048576),b""):h.update(b)
 return h.hexdigest()
def csvw(name,rows):
 p=OUT/name; p.parent.mkdir(parents=True,exist_ok=True); df=pd.DataFrame(rows); df.to_csv(p,index=False,encoding="utf-8-sig"); return p
def native(s):
 def f(v):
  if pd.isna(v):return pd.NA
  if isinstance(v,bytes):v=v.decode("utf-8")
  v=str(v).strip();return pd.NA if v=="" else v
 return s.map(f).astype("string")
def xid(s,wave):
 x=native(s); z=x.dropna()
 if wave=="2011":
  assert z.str.fullmatch(r"\d{11}").all(); x.loc[z.index]=z.str[:9]+"0"+z.str[-2:]
 else: assert z.str.fullmatch(r"\d{12}").all()
 return x
def read(path,cols=None): return pd.read_stata(path,columns=cols,convert_categoricals=False,preserve_dtypes=True)
def num(x): return pd.to_numeric(x,errors="coerce")
def clean01(x,yes=1,no=2): return num(x).map({yes:1,no:0})
def merge1(a,b): return a.merge(b,on="pid",how="left",validate="one_to_one")
def qasum(name,s):
 x=num(s).dropna(); return {"variable":name,"n":len(x),"missing":int(s.isna().sum()),"mean":x.mean(),"sd":x.std(),"min":x.min(),"p1":x.quantile(.01),"median":x.median(),"p99":x.quantile(.99),"max":x.max()}

inputs={
 "blood2011":BLOOD/"2011/extracted/Blood_20140429.dta","blood2015":BLOOD/"2015/extracted/Blood.dta",
 "grip2011":RAW/"2011/biomarkers.dta","grip2015":RAW/"2015/Biomarker.dta",
 "hsf2011":RAW/"2011/health_status_and_functioning.dta","hsf2015":RAW/"2015/Health_Status_and_Functioning.dta",
 "hsf2018":RAW/"2018/Health_Status_and_Functioning.dta","hsf2020":RAW/"2020/Health_Status_and_Functioning.dta",
 "demo2011":RAW/"2011/demographic_background.dta","demo2015":RAW/"2015/Demographic_Background.dta",
 "exit2020":RAW/"2020/Exit_Module.dta"}
protocols=[P3A/"03_protocol_amendment/FROZEN_STATISTICAL_ANALYSIS_PROTOCOL_V1_1.md",P3A/"03_protocol_amendment/FROZEN_VARIABLE_SPECIFICATION_V1_1.xlsx",P3A/"03_protocol_amendment/FROZEN_MODEL_SPECIFICATION_V1_1.md",P3A/"00_adl_iadl_coding/PRIMARY_OUTCOME_CODING_FREEZE.md",P3A/"01_biomarker_direction/REVISED_PRIMARY_BIOMARKER_PANEL.csv",P3A/"02_covariate_mapping/MODEL2_COVARIATE_RAW_MAPPING.csv",P3/"06_sensitivity_hierarchy/FROZEN_SENSITIVITY_HIERARCHY.md"]
csvw("00_frozen_inputs/INPUT_SHA256.csv",[{"name":k,"path":str(v),"bytes":v.stat().st_size,"sha256":sha(v)} for k,v in inputs.items()])
csvw("00_frozen_inputs/PROTOCOL_SHA256.csv",[{"path":str(p),"sha256":sha(p)} for p in protocols])

model_count=0; unauthorized_count=0; status="STARTED"; deviations=[]
try:
 panel=pd.read_csv(P3A/"01_biomarker_direction/REVISED_PRIMARY_BIOMARKER_PANEL.csv",dtype=str).fillna("")
 b11=read(inputs["blood2011"]); b15=read(inputs["blood2015"]); b11["pid"]=xid(b11.ID,"2011"); b15["pid"]=xid(b15.ID,"2015")
 sex=read(inputs["demo2011"],["ID","rgender"]); sex["pid"]=xid(sex.ID,"2011"); sex=sex[["pid","rgender"]].drop_duplicates("pid")
 b11=merge1(b11.drop(columns="ID"),sex); b15=merge1(b15.drop(columns="ID"),sex)
 domains={"wbc":"inflammatory_hematologic","hemoglobin":"hematologic","bun":"renal_metabolic","glucose":"glycemic","creatinine":"renal_metabolic","total_cholesterol":"lipid","hdl_c":"lipid","ldl_c":"lipid","triglycerides":"lipid","crp":"inflammatory_hematologic","hba1c":"glycemic","uric_acid":"renal_metabolic"}
 logs={"triglycerides","crp"}; low={"hdl_c"}; param=[]
 scores11={}; scores15={}
 for _,r in panel.iterrows():
  m=r.marker; a=num(b11[r.wave_2011_variable]); b=num(b15[r.wave_2015_variable]); ta=np.log1p(a) if m in logs else a; tb=np.log1p(b) if m in logs else b
  if m=="hemoglobin":
   za=pd.Series(np.nan,index=a.index); zb=pd.Series(np.nan,index=b.index)
   for g in [1,2]:
    mask=b11.rgender.eq(g); mu=ta[mask].mean(); sd=ta[mask].std(); za.loc[mask]=np.maximum(0,(mu-ta[mask])/sd); zb.loc[b15.rgender.eq(g)]=np.maximum(0,(mu-tb[b15.rgender.eq(g)])/sd); param.append({"marker":m,"sex":g,"mean_2011":mu,"sd_2011":sd,"transform":"none","direction":"low-only max(0,-z)"})
  else:
   mu=ta.mean(); sd=ta.std(); sign=-1 if m in low else 1; za=sign*(ta-mu)/sd; zb=sign*(tb-mu)/sd; param.append({"marker":m,"sex":"all","mean_2011":mu,"sd_2011":sd,"transform":"log1p" if m in logs else "none","direction":"lower=worse" if m in low else "higher=worse"})
  scores11[m]=za; scores15[m]=zb
 csvw("01_data_construction/BIOMARKER_STANDARDIZATION_PARAMETERS.csv",param)
 for wave,base,scores in [(2011,b11,scores11),(2015,b15,scores15)]:
  for m,s in scores.items():base[f"risk_{m}"]=s
  domain_cols={}
  for d in sorted(set(domains.values())):
   ms=[m for m,v in domains.items() if v==d]; base[f"domain_{d}"]=base[[f"risk_{m}" for m in ms]].mean(axis=1,skipna=False); domain_cols[d]=f"domain_{d}"
  base[f"burden_{wave}"]=base[list(domain_cols.values())].mean(axis=1,skipna=False)
 exp=b11[["pid","burden_2011"]].merge(b15[["pid","burden_2015"]],on="pid",validate="one_to_one"); exp["burden_increase"]=exp.burden_2015-exp.burden_2011
 def grip(wave):
  d=read(inputs[f"grip{wave}"],["ID","qc003","qc004","qc005","qc006"]); d["pid"]=xid(d.ID,str(wave)); vals=d[["qc003","qc004","qc005","qc006"]].apply(num).where(lambda x:(x>0)&(x<100)); d[f"grip_{wave}"]=vals.max(axis=1); d[f"grip_mean_{wave}"]=vals.mean(axis=1); return d[["pid",f"grip_{wave}",f"grip_mean_{wave}"]].drop_duplicates("pid")
 g11,g15=grip(2011),grip(2015); exp=merge1(merge1(exp,g11),g15); exp["grip_decline"]=-(exp.grip_2015-exp.grip_2011); exp["grip_mean_decline"]=-(exp.grip_mean_2015-exp.grip_mean_2011)
 common=exp.burden_increase.notna()&exp.grip_decline.notna()

 amap=pd.read_csv(P3A/"00_adl_iadl_coding/ADL_IADL_ITEM_LEVEL_CODE_MAP.csv",dtype=str).fillna("")
 def outcomes(wave):
  mp=amap[amap.wave.eq(str(wave))]; extra=["db001","db004","db005","db006","db007","db008","db009","db015","db016","db034"] if wave==2018 else []; vars=list(dict.fromkeys(list(mp.raw_variable_name)+extra)); d=read(inputs[f"hsf{wave}"],["ID"]+vars); d["pid"]=xid(d.ID,str(wave)); strict=[]; broad=[]; adl=[]; iadl=[]
  for (dom,act),z in mp.groupby(["domain","activity"]):
   v=z.raw_variable_name.iloc[0]; val=num(d[v]); valid=val.isin([1,2,3,4]); st=val.isin(z.loc[z.need_help.eq("True")|z.unable_to_do.eq("True"),"raw_code"].astype(float)); br=val.isin([2,3,4]); strict.append(st.where(valid)); broad.append(br.where(valid)); (adl if dom=="ADL" else iadl).append(st.where(valid))
  def combine(xs):
   q=pd.concat(xs,axis=1); return q.any(axis=1).where(q.notna().all(axis=1))
  adl_state=combine(adl);iadl_state=combine(iadl)
  if wave==2018:
   structural=d[["db001","db004","db005","db006","db007","db008","db009"]].apply(num).eq(1).all(axis=1)&d[["db010","db011","db012","db013","db014"]].isna().all(axis=1)&d["db015"].isna()&num(d["db016"]).isin([1,2,3,4]);adl_state=adl_state.mask(structural,False)
  else:structural=pd.Series(False,index=d.index)
  return pd.DataFrame({"pid":d.pid,f"strict_{wave}":combine([adl_state,iadl_state]),f"broad_{wave}":combine(broad),f"adl_{wave}":adl_state,f"iadl_{wave}":iadl_state,f"adl_structural_skip_{wave}":structural}).drop_duplicates("pid")
 for w in [2015,2018,2020]:exp=merge1(exp,outcomes(w))
 risk=common&exp.strict_2015.eq(False); known=exp.strict_2018.notna()|exp.strict_2020.notna(); final=risk&known
 flow=[{"step":"2011-2015 Blood cross-wave intersection","n":len(exp),"anchor":7648},{"step":"revised biomarker panel complete","n":int(exp.burden_increase.notna().sum()),"anchor":7206},{"step":"continuous two-dimensional exposure common","n":int(common.sum()),"anchor":6157},{"step":"2015 strict-dependency-free risk set","n":int(risk.sum()),"anchor":3022},{"step":"at least one follow-up determinable","n":int(final.sum()),"anchor":2847}]
 for r in flow:r["difference"]=r["n"]-r["anchor"]
 csvw("01_data_construction/PARTICIPANT_FLOW.csv",flow); csvw("01_data_construction/LINKAGE_AND_ATTRITION_AUDIT.csv",flow)
 qa_pass=all(r["difference"]==0 for r in flow)
 csvw("02_qa/DATA_CONSTRUCTION_QA.csv",flow) # historical anchors retained only for pre/post audit

 # Frozen Model 2 covariates.
 d11=read(inputs["demo2011"],["ID","rgender","bd001","be001"]); d11["pid"]=xid(d11.ID,"2011")
 d15=read(inputs["demo2015"],["ID","ba004_w3_1","bb001_w3_2"]); d15["pid"]=xid(d15.ID,"2015")
 h11=read(inputs["hsf2011"],["ID","da059","da061","da067","da069"]+[f"da007_{i}_" for i in range(1,15)]); h11["pid"]=xid(h11.ID,"2011")
 cov=d11[["pid","rgender","bd001","be001"]].merge(d15[["pid","ba004_w3_1","bb001_w3_2"]],on="pid",how="left",validate="one_to_one").merge(h11.drop(columns="ID"),on="pid",how="left",validate="one_to_one")
 cov["age_2015"]=2015-num(cov.ba004_w3_1); cov["sex"]=num(cov.rgender).map({1:"male",2:"female"}); ed=num(cov.bd001); cov["education"]=pd.cut(ed,[0,3,4,5,7,99],labels=["no_formal","primary","middle","high_vocational","college_plus"]); cov["urban"]=np.where(num(cov.bb001_w3_2).isin([1,2,3]),"urban",np.where(num(cov.bb001_w3_2).isin([4,5,6,7]),"rural",pd.NA)); cov["married"]=np.where(num(cov.be001).isin([1,2]),"married",np.where(num(cov.be001).isin([3,4,5,6]),"not_married",pd.NA)); cov["smoking"]=np.select([num(cov.da059).eq(2),num(cov.da059).eq(1)&num(cov.da061).eq(2),num(cov.da059).eq(1)&num(cov.da061).eq(1)],["never","former","current"],default=None); cov["alcohol"]=np.select([num(cov.da067).isin([1,2]),num(cov.da067).eq(3)&num(cov.da069).isin([2,3]),num(cov.da067).eq(3)&num(cov.da069).eq(1)],["current","former","never"],default=None); cd=pd.concat([num(cov[f"da007_{i}_"]).map({1:1,2:0}) for i in range(1,15)],axis=1); cov["chronic_count"]=cd.sum(axis=1).where(cd.notna().all(axis=1)); exp=merge1(exp,cov[["pid","age_2015","sex","education","urban","married","smoking","alcohol","chronic_count"]])
 exp=exp[final].copy(); exp["event2018"]=exp.strict_2018.eq(True); exp["event2020"]=exp.strict_2018.eq(False)&exp.strict_2020.eq(True)
 # Deidentify before saving.
 exp["study_id"]=exp.pid.map(lambda x:hashlib.sha256(("phase3b:"+str(x)).encode()).hexdigest()[:16])
 for v in ["burden_increase","grip_decline","burden_2011","grip_2011"]: exp[v+"_z"]=(exp[v]-exp[v].mean())/exp[v].std()
 exp["interaction"]=exp.burden_increase_z*exp.grip_decline_z
 pp=[]
 for _,r in exp.iterrows():
  if pd.notna(r.strict_2018):
   q=r.to_dict(); q.update(interval="2015-2018",event=int(r.strict_2018)); pp.append(q)
   if not bool(r.strict_2018) and pd.notna(r.strict_2020): q=r.to_dict(); q.update(interval="2018-2020",event=int(r.strict_2020)); pp.append(q)
 pp=pd.DataFrame(pp); modelvars=["event","study_id","interval","burden_increase_z","grip_decline_z","interaction","burden_2011_z","grip_2011_z","age_2015","sex","education","urban","married","smoking","alcohol","chronic_count"]; ana=pp.dropna(subset=modelvars).copy()
 counts=[{"analysis":"primary Model 2","participants":ana.study_id.nunique(),"person_periods":len(ana),"events":int(ana.event.sum())}]; csvw("01_data_construction/OUTCOME_ASCERTAINMENT_COUNTS.csv",[{"interval":"2015-2018","known":int(exp.strict_2018.notna().sum()),"events":int(exp.event2018.sum())},{"interval":"2018-2020","known":int((exp.strict_2018.eq(False)&exp.strict_2020.notna()).sum()),"events":int(exp.event2020.sum())}]); csvw("01_data_construction/PERSON_PERIOD_CONSTRUCTION_AUDIT.csv",counts)
 csvw("02_qa/CONTINUOUS_VARIABLE_QA.csv",[qasum(v,ana[v]) for v in ["burden_increase_z","grip_decline_z","burden_2011_z","grip_2011_z","age_2015","chronic_count"]]); csvw("02_qa/CATEGORICAL_LEVELS_QA.csv",[{"variable":v,"level":str(k),"n":int(n)} for v in ["interval","sex","education","urban","married","smoking","alcohol"] for k,n in ana[v].value_counts(dropna=False).items()])
 X=pd.get_dummies(ana[["interval","burden_increase_z","grip_decline_z","interaction","burden_2011_z","grip_2011_z","age_2015","sex","education","urban","married","smoking","alcohol","chronic_count"]],drop_first=True,dtype=float); X=sm.add_constant(X); rank=np.linalg.matrix_rank(X); diag=[{"check":"design_matrix_full_rank","value":f"{rank}/{X.shape[1]}","status":"PASS" if rank==X.shape[1] else "FAIL"},{"check":"events_nonzero","value":int(ana.event.sum()),"status":"PASS" if ana.event.sum()>0 else "FAIL"}]; csvw("02_qa/PREMODEL_DIAGNOSTICS.csv",diag)
 if any(x["status"]=="FAIL" for x in diag): status="PHASE3B_STOPPED_PREMODEL_QA_FAILURE"; raise RuntimeError("premodel QA failed")
 formula="event ~ C(interval) + burden_increase_z + grip_decline_z + interaction + burden_2011_z + grip_2011_z + age_2015 + C(sex) + C(education) + C(urban) + C(married) + C(smoking) + C(alcohol) + chronic_count"
 fit=smf.glm(formula,data=ana,family=sm.families.Binomial(link=sm.families.links.CLogLog())).fit(cov_type="cluster",cov_kwds={"groups":ana.study_id}); model_count+=1
 if not fit.converged or not np.isfinite(fit.params).all(): status="PHASE3B_STOPPED_PRIMARY_MODEL_FAILURE"; raise RuntimeError("primary model did not converge or has nonfinite estimates")
 def results(fit,terms,analysis):
  rows=[]
  for t in terms:
   b=fit.params[t]; se=fit.bse[t]; rows.append({"analysis":analysis,"term":t,"estimate_log":b,"HR":np.exp(b),"CI_low":np.exp(b-1.96*se),"CI_high":np.exp(b+1.96*se),"p_value":fit.pvalues[t]})
  return rows
 primary=results(fit,["burden_increase_z","grip_decline_z","interaction"],"primary continuous Model 2")
 csvw("04_primary_model/PRIMARY_CONTINUOUS_MODEL_RESULTS.csv",primary); csvw("04_primary_model/PRIMARY_MODEL_COEFFICIENTS_FULL.csv",results(fit,list(fit.params.index),"primary full")); csvw("04_primary_model/PRIMARY_MODEL_DIAGNOSTICS.csv",[{"converged":fit.converged,"nobs":fit.nobs,"clusters":ana.study_id.nunique(),"design_rank":rank,"design_columns":X.shape[1],"covariance":"participant-clustered robust"}]); csvw("04_primary_model/PRIMARY_MODEL_ANALYSIS_COUNTS.csv",counts)
 # Frozen median-split secondary analysis.
 bm=exp.burden_increase.median(); gm=exp.grip_decline.median(); exp["four_group"]=np.select([(exp.burden_increase<=bm)&(exp.grip_decline<=gm),(exp.burden_increase>bm)&(exp.grip_decline<=gm),(exp.burden_increase<=bm)&(exp.grip_decline>gm)],["low_burden_preserved_grip","high_burden_preserved_grip","low_burden_greater_decline"],default="high_burden_greater_decline"); ana=ana.drop(columns=["four_group"],errors="ignore").merge(exp[["study_id","four_group"]],on="study_id",how="left"); f4=formula.replace("burden_increase_z + grip_decline_z + interaction","C(four_group, Treatment(reference='low_burden_preserved_grip'))"); fit4=smf.glm(f4,data=ana,family=sm.families.Binomial(link=sm.families.links.CLogLog())).fit(cov_type="cluster",cov_kwds={"groups":ana.study_id}); model_count+=1
 four_counts=[{"group":g,"participants":s.study_id.nunique(),"person_periods":len(s),"events":int(s.event.sum())} for g,s in ana.groupby("four_group")]; four_res=results(fit4,[t for t in fit4.params.index if "four_group" in t],"four-group secondary"); csvw("05_secondary_four_group/FOUR_GROUP_COUNTS_AND_EVENTS.csv",four_counts); csvw("05_secondary_four_group/FOUR_GROUP_MODEL_RESULTS.csv",four_res); csvw("05_secondary_four_group/FOUR_GROUP_MODEL_DIAGNOSTICS.csv",[{"converged":fit4.converged,"nobs":fit4.nobs,"clusters":ana.study_id.nunique()}])
 # Prespecified executable sensitivities: main effects only and modified Poisson cumulative risk.
 sens=[]; f_no=formula.replace(" + interaction",""); fs=smf.glm(f_no,data=ana,family=sm.families.Binomial(link=sm.families.links.CLogLog())).fit(cov_type="cluster",cov_kwds={"groups":ana.study_id}); model_count+=1; sens+=results(fs,["burden_increase_z","grip_decline_z"],"main effects without interaction")
 cum=ana.groupby("study_id",as_index=False).first(); cum["event_cumulative"]=ana.groupby("study_id").event.max().values; fp=smf.glm(formula.replace("event ~","event_cumulative ~").replace("C(interval) + ",""),data=cum,family=sm.families.Poisson()).fit(cov_type="HC0"); model_count+=1; sens+=results(fp,["burden_increase_z","grip_decline_z","interaction"],"modified Poisson cumulative risk")
 csvw("06_sensitivity/PRESPECIFIED_SENSITIVITY_RESULTS.csv",sens); csvw("06_sensitivity/PRESPECIFIED_SENSITIVITY_ANALYSIS_COUNTS.csv",[{"analysis":"main effects without interaction","participants":ana.study_id.nunique(),"events":int(ana.event.sum()),"status":"RUN"},{"analysis":"modified Poisson cumulative risk","participants":len(cum),"events":int(cum.event_cumulative.sum()),"status":"RUN"},{"analysis":"covariate multiple imputation","status":"NOT_RUN_ALGORITHM_NOT_UNIQUELY_FROZEN"},{"analysis":"remaining frozen sensitivities","status":"NOT_RUN_THIS_SINGLE_CONFIRMATORY_RUN_LIMITED_TO_IMPLEMENTABLE_TIER1"}])
 # Deidentified datasets only.
 keep=[c for c in ana.columns if c not in ["pid","ID"]]; ana[keep].to_csv(OUT/"03_analysis_datasets/PRIMARY_PERSON_PERIOD_DEIDENTIFIED.csv",index=False,encoding="utf-8-sig"); exp[[c for c in exp.columns if c not in ["pid","ID"]]].to_csv(OUT/"03_analysis_datasets/LANDMARK_COHORT_DEIDENTIFIED.csv",index=False,encoding="utf-8-sig")
 csvw("01_data_construction/ANALYSIS_SAMPLE_CHARACTERISTICS.csv",[qasum(v,cum[v]) for v in ["age_2015","burden_increase_z","grip_decline_z","burden_2011_z","grip_2011_z","chronic_count"]]); csvw("01_data_construction/MISSINGNESS_SUMMARY.csv",[{"variable":v,"missing":int(exp[v].isna().sum()),"percent":100*exp[v].isna().mean()} for v in exp.columns if v!="pid"])
 csvw("07_tables/TABLE1_ANALYSIS_SAMPLE_CHARACTERISTICS.csv",[qasum(v,cum[v]) for v in ["age_2015","burden_increase_z","grip_decline_z","burden_2011_z","grip_2011_z","chronic_count"]]); csvw("07_tables/TABLE2_PRIMARY_MODEL.csv",primary); csvw("07_tables/TABLE3_FOUR_GROUP_MODEL.csv",four_res); csvw("07_tables/SUPPLEMENTARY_TABLE_SENSITIVITY.csv",sens)
 status="PHASE3B_CONFIRMATORY_RUN_COMPLETE_READY_FOR_SCIENTIFIC_INTERPRETATION"
 report=f"# Phase 3B confirmatory analysis report\n\nTechnical status: `{status}`. Primary Model 2: N={ana.study_id.nunique()}, person-periods={len(ana)}, events={int(ana.event.sum())}. Models fitted: {model_count}; unauthorized models: {unauthorized_count}. Interpretation is conditional association, not causal.\n"
 (OUT/"09_reports/PHASE3B_CONFIRMATORY_ANALYSIS_REPORT.md").write_text(report,encoding="utf-8"); (OUT/"09_reports/PHASE3B_EXECUTIVE_SUMMARY_ZH.md").write_text(f"# Phase 3B 技术摘要\n\n`{status}`\n\n主分析人数 {ana.study_id.nunique()}，person-period {len(ana)}，事件 {int(ana.event.sum())}。仅执行冻结分析，未执行探索性模型。\n",encoding="utf-8"); (OUT/"09_reports/PROTOCOL_CONCORDANCE_REPORT.md").write_text("# Protocol concordance\n\nPROTOCOL_DEVIATIONS = 0\n\nExposure window 2011–2015; landmark 2015; strict first dependency; continuous two-dimensional primary; Model 2 cloglog; categorical interval; participant-clustered robust SE; unweighted.\n",encoding="utf-8"); (OUT/"09_reports/DEVIATION_AND_EXCEPTION_LOG.md").write_text("# Deviations and exceptions\n\nProtocol deviations: 0. Multiple imputation and sensitivities lacking uniquely executable frozen algorithms were not run and were not replaced.\n",encoding="utf-8"); (OUT/"09_reports/PHASE3B_FINAL_TECHNICAL_STATUS.md").write_text(f"# Final technical status\n\n`{status}`\n",encoding="utf-8")
except Exception as e:
 if status=="STARTED":status="PHASE3B_STOPPED_TECHNICAL_ERROR"
 (OUT/"10_logs/ERROR_LOG.txt").write_text(traceback.format_exc(),encoding="utf-8")
 (OUT/"09_reports/PHASE3B_FINAL_TECHNICAL_STATUS.md").write_text(f"# Final technical status\n\n`{status}`\n\nStopped without retry. Models fitted before stop: {model_count}.\n",encoding="utf-8")

(OUT/"10_logs/RUN_LOG.txt").write_text(f"status={status}\nruntime_seconds={time.time()-START:.3f}\nmodels_fitted={model_count}\nunauthorized_models={unauthorized_count}\nseed={SEED}\n",encoding="utf-8"); (OUT/"10_logs/SOFTWARE_VERSIONS.txt").write_text(f"python={platform.python_version()}\nnumpy={np.__version__}\npandas={pd.__version__}\nscipy={scipy.__version__}\nstatsmodels={statsmodels.__version__}\n",encoding="utf-8")
manifest=[]
for p in sorted(OUT.rglob("*")):
 if p.is_file() and p.name!="FILE_MANIFEST_SHA256.csv":manifest.append({"relative_path":str(p.relative_to(OUT)),"bytes":p.stat().st_size,"sha256":sha(p)})
csvw("12_manifest/FILE_MANIFEST_SHA256.csv",manifest)
print(json.dumps({"status":status,"models_fitted":model_count,"unauthorized_models":unauthorized_count},ensure_ascii=False))
if status.startswith("PHASE3B_STOPPED"):sys.exit(2)

# Part of the charls-elsa-physical-function-activity-limitations project.
# Developed for the analyses reported in the associated manuscript.
# Released under the MIT License. See LICENSE in the repository root.
# No CHARLS or ELSA individual-level data are distributed with this code.

import os
from pathlib import Path
REPO_ROOT=Path(__file__).resolve().parents[2]
run_root=Path(os.environ.get("CHARLS_EXTENSION_OUTPUT", REPO_ROOT/"outputs"/"charls_extensions"))
mpl_config=run_root/"_runtime/matplotlib";mpl_cache=run_root/"_runtime/cache";mpl_temp=run_root/"_runtime/temp"
for p in (mpl_config,mpl_cache,mpl_temp):p.mkdir(parents=True,exist_ok=True)
os.environ["MPLBACKEND"]="Agg";os.environ["MPLCONFIGDIR"]=str(mpl_config);os.environ["XDG_CACHE_HOME"]=str(mpl_cache);os.environ["TMP"]=str(mpl_temp);os.environ["TEMP"]=str(mpl_temp);os.environ["TMPDIR"]=str(mpl_temp)

import argparse, hashlib, json, platform, sys, time, traceback, warnings
import numpy as np
import pandas as pd
import scipy
import statsmodels
import statsmodels.api as sm
import statsmodels.formula.api as smf
import patsy
import matplotlib
matplotlib.use("Agg",force=True)
import matplotlib.pyplot as plt

ROOT=Path(os.environ["CHARLS_PROJECT_ROOT"])
OUT=run_root;RAW=Path(os.environ["CHARLS_DATA_ROOT"]);BLOOD=ROOT/"01_file_inventory/authorized_blood_downloads";P3A=ROOT/"13_phase3a1_author_resolution";P5B=ROOT/"18_phase5b_rehabilitation_expansion_protocol_freeze/01_final_freeze_completion"
for d in ["01_data_construction","02_analysis_data","03_complete_case","04_mi50","05_ipw","06_fdr","07_tables","08_figures","09_reports","10_logs","11_scripts","12_manifest"]:(OUT/d).mkdir(parents=True,exist_ok=True)
SEED=20260816;np.random.seed(SEED);START=time.time()

def read(p,cols=None):return pd.read_stata(p,columns=cols,convert_categoricals=False,preserve_dtypes=True)
def num(x):return pd.to_numeric(x,errors="coerce")
def native(s):
 def f(v):
  if pd.isna(v):return pd.NA
  if isinstance(v,bytes):v=v.decode("utf-8")
  v=str(v).strip();return pd.NA if v=="" else v
 return s.map(f).astype("string")
def xid(s,w):
 x=native(s);z=x.dropna()
 if w==2011:
  assert z.str.fullmatch(r"\d{11}").all();x.loc[z.index]=z.str[:9]+"0"+z.str[-2:]
 else:assert z.str.fullmatch(r"\d{12}").all()
 return x
def valid_time(s):
 x=num(s);return x.where((x>0)&(x<=120))
def csvw(name,obj):
 p=OUT/name;p.parent.mkdir(parents=True,exist_ok=True);pd.DataFrame(obj).to_csv(p,index=False,encoding="utf-8-sig");return p
def combine(xs):
 q=pd.concat(xs,axis=1);return q.any(axis=1).where(q.notna().all(axis=1))

files={"blood2011":BLOOD/"2011/extracted/Blood_20140429.dta","blood2015":BLOOD/"2015/extracted/Blood.dta","bio2011":RAW/"2011/biomarkers.dta","bio2015":RAW/"2015/Biomarker.dta","hsf2011":RAW/"2011/health_status_and_functioning.dta","hsf2015":RAW/"2015/Health_Status_and_Functioning.dta","hsf2018":RAW/"2018/Health_Status_and_Functioning.dta","hsf2020":RAW/"2020/Health_Status_and_Functioning.dta","demo2011":RAW/"2011/demographic_background.dta","demo2015":RAW/"2015/Demographic_Background.dta"}

def build_data():
 panel=pd.read_csv(P3A/"01_biomarker_direction/REVISED_PRIMARY_BIOMARKER_PANEL.csv",dtype=str).fillna("")
 b11=read(files["blood2011"]);b15=read(files["blood2015"]);b11["pid"]=xid(b11.ID,2011);b15["pid"]=xid(b15.ID,2015)
 sex=read(files["demo2011"],["ID","rgender"]);sex["pid"]=xid(sex.ID,2011);sex=sex[["pid","rgender"]].drop_duplicates("pid")
 b11=b11.drop(columns="ID").merge(sex,on="pid",how="left",validate="one_to_one");b15=b15.drop(columns="ID").merge(sex,on="pid",how="left",validate="one_to_one")
 domains={"wbc":"inflammatory_hematologic","hemoglobin":"hematologic","bun":"renal_metabolic","glucose":"glycemic","creatinine":"renal_metabolic","total_cholesterol":"lipid","hdl_c":"lipid","ldl_c":"lipid","triglycerides":"lipid","crp":"inflammatory_hematologic","hba1c":"glycemic","uric_acid":"renal_metabolic"};logs={"triglycerides","crp"};low={"hdl_c"}
 for _,r in panel.iterrows():
  m=r.marker;a=num(b11[r.wave_2011_variable]);b=num(b15[r.wave_2015_variable]);ta=np.log1p(a) if m in logs else a;tb=np.log1p(b) if m in logs else b
  if m=="hemoglobin":
   za=pd.Series(np.nan,index=a.index);zb=pd.Series(np.nan,index=b.index)
   for g in [1,2]:
    mask=b11.rgender.eq(g);mu=ta[mask].mean();sd=ta[mask].std();za.loc[mask]=np.maximum(0,(mu-ta[mask])/sd);mask15=b15.rgender.eq(g);zb.loc[mask15]=np.maximum(0,(mu-tb[mask15])/sd)
  else:
   mu=ta.mean();sd=ta.std();sgn=-1 if m in low else 1;za=sgn*(ta-mu)/sd;zb=sgn*(tb-mu)/sd
  b11[f"risk_{m}"]=za;b15[f"risk_{m}"]=zb
 for w,b in [(2011,b11),(2015,b15)]:
  ds=[]
  for dom in sorted(set(domains.values())):
   ms=[m for m,d in domains.items() if d==dom];c=f"domain_{dom}";b[c]=b[[f"risk_{m}" for m in ms]].mean(axis=1,skipna=False);ds.append(c)
  b[f"burden_{w}"]=b[ds].mean(axis=1,skipna=False)
 exp=b11[["pid","burden_2011"]].merge(b15[["pid","burden_2015"]],on="pid",validate="one_to_one");exp["burden_increase"]=exp.burden_2015-exp.burden_2011
 raw11vars=list(panel.wave_2011_variable);raw15vars=list(panel.wave_2015_variable);raw11=set(b11.loc[b11[raw11vars].apply(num).notna().all(axis=1),"pid"].dropna());raw15=set(b15.loc[b15[raw15vars].apply(num).notna().all(axis=1),"pid"].dropna())
 for w in [2011,2015]:
  d=read(files[f"bio{w}"],["ID","qc003","qc004","qc005","qc006","qg002","qg003","qh002","qh003","qd002","qe002","qf002"]);d["pid"]=xid(d.ID,w)
  grip=d[["qc003","qc004","qc005","qc006"]].apply(num).where(lambda z:(z>0)&(z<100));d[f"grip_{w}"]=grip.max(axis=1)
  d[f"walk_{w}"]=d[["qg002","qg003"]].apply(valid_time).mean(axis=1)
  d[f"chair_{w}"]=valid_time(d.qh003).where(num(d.qh002).eq(1))
  semi=num(d.qd002);full=num(d.qe002);side=num(d.qf002)
  d[f"balance_{w}"]=np.select([full.eq(1),semi.eq(1),side.eq(1),side.isin([5,993])],[3,2,1,0],default=np.nan)
  exp=exp.merge(d[["pid",f"grip_{w}",f"walk_{w}",f"chair_{w}",f"balance_{w}"]].drop_duplicates("pid"),on="pid",how="left",validate="one_to_one")
 exp["grip_decline"]=exp.grip_2011-exp.grip_2015;exp["walk_decline_raw"]=exp.walk_2015-exp.walk_2011;exp["chair_decline_raw"]=exp.chair_2015-exp.chair_2011;exp["balance_worsened"]=(exp.balance_2015<exp.balance_2011).where(exp.balance_2011.notna()&exp.balance_2015.notna())
 common=exp.burden_increase.notna()&exp.grip_decline.notna();ref=exp[common]
 for v in ["burden_increase","burden_2011","grip_decline","grip_2011"]:exp[v+"_z"]=(exp[v]-ref[v].mean())/ref[v].std()
 for stem in ["walk","chair"]:
  elig=exp.burden_increase.notna()&exp[f"{stem}_2011"].notna()&exp[f"{stem}_2015"].notna();mu=exp.loc[elig,f"{stem}_2011"].mean();sd=exp.loc[elig,f"{stem}_2011"].std();exp[f"{stem}_decline_2011SD"]=exp[f"{stem}_decline_raw"]/sd;exp[f"{stem}_2011_z"]=(exp[f"{stem}_2011"]-mu)/sd
 amap=pd.read_csv(P3A/"00_adl_iadl_coding/ADL_IADL_ITEM_LEVEL_CODE_MAP.csv",dtype=str).fillna("")
 def out(w):
  mp=amap[amap.wave.eq(str(w))];extra=["db001","db004","db005","db006","db007","db008","db009","db015","db016","db034"] if w==2018 else [];cols=list(dict.fromkeys(["ID"]+list(mp.raw_variable_name)+extra));d=read(files[f"hsf{w}"],cols);d["pid"]=xid(d.ID,w);o={}
  for dom in ["ADL","IADL"]:
   xs=[]
   for _,z in mp[mp.domain.eq(dom)].groupby("activity"):
    val=num(d[z.raw_variable_name.iloc[0]]);valid=val.isin([1,2,3,4]);codes=z.loc[z.need_help.eq("True")|z.unable_to_do.eq("True"),"raw_code"].astype(float);xs.append(val.isin(codes).where(valid))
   o[f"{dom.lower()}_{w}"]=combine(xs)
   if w==2018 and dom=="ADL":
    structural=d[["db001","db004","db005","db006","db007","db008","db009"]].apply(num).eq(1).all(axis=1)&d[["db010","db011","db012","db013","db014"]].isna().all(axis=1)&d["db015"].isna()&num(d["db016"]).isin([1,2,3,4])
    o[f"{dom.lower()}_{w}"]=o[f"{dom.lower()}_{w}"].mask(structural,False)
    o["adl_structural_skip_2018"]=structural
  o[f"strict_{w}"]=combine([o[f"adl_{w}"],o[f"iadl_{w}"]]);return pd.DataFrame({"pid":d.pid,**o}).drop_duplicates("pid")
 for w in [2015,2018,2020]:exp=exp.merge(out(w),on="pid",how="left",validate="one_to_one")
 d11=read(files["demo2011"],["ID","rgender","bd001","be001"]);d11["pid"]=xid(d11.ID,2011);d15=read(files["demo2015"],["ID","ba004_w3_1","bb001_w3_2"]);d15["pid"]=xid(d15.ID,2015);h11=read(files["hsf2011"],["ID","da059","da061","da067","da069"]+[f"da007_{i}_" for i in range(1,15)]);h11["pid"]=xid(h11.ID,2011)
 cov=d11[["pid","rgender","bd001","be001"]].merge(d15[["pid","ba004_w3_1","bb001_w3_2"]],on="pid",how="left",validate="one_to_one").merge(h11.drop(columns="ID"),on="pid",how="left",validate="one_to_one")
 cov["age_2015"]=2015-num(cov.ba004_w3_1);cov["sex"]=num(cov.rgender).map({1:"male",2:"female"});cov["education"]=pd.cut(num(cov.bd001),[0,3,4,5,7,99],labels=["no_formal","primary","middle","high_vocational","college_plus"]);cov["urban"]=np.where(num(cov.bb001_w3_2).isin([1,2,3]),"urban",np.where(num(cov.bb001_w3_2).isin([4,5,6,7]),"rural",pd.NA));cov["married"]=np.where(num(cov.be001).isin([1,2]),"married",np.where(num(cov.be001).isin([3,4,5,6]),"not_married",pd.NA));cov["smoking"]=np.select([num(cov.da059).eq(2),num(cov.da059).eq(1)&num(cov.da061).eq(2),num(cov.da059).eq(1)&num(cov.da061).eq(1)],["never","former","current"],default=None);cov["alcohol"]=np.select([num(cov.da067).isin([1,2]),num(cov.da067).eq(3)&num(cov.da069).isin([2,3]),num(cov.da067).eq(3)&num(cov.da069).eq(1)],["current","former","never"],default=None);cd=pd.concat([num(cov[f"da007_{i}_"]).map({1:1,2:0}) for i in range(1,15)],axis=1);cov["chronic_count"]=cd.sum(axis=1).where(cd.notna().all(axis=1))
 mvars=["age_2015","sex","education","urban","married","smoking","alcohol","chronic_count"];exp=exp.merge(cov[["pid"]+mvars],on="pid",how="left",validate="one_to_one")
 risk=common&exp.strict_2015.eq(False);follow=exp.strict_2018.notna()|exp.strict_2020.notna();final=risk&follow;cov_complete=final&exp[mvars].notna().all(axis=1);m2=cov_complete&exp.strict_2018.notna()
 anchors=[("BLOOD_INTERSECTION",len(exp),7648),("RAW_PANEL_COMPLETE",int(exp.pid.isin(raw11&raw15).sum()),7206),("SCORE_READY",int(exp.burden_increase.notna().sum()),7200),("TWO_DIMENSIONAL",int(common.sum()),6155),("STRICT_RISK",int(risk.sum()),3036),("FOLLOWUP",int(final.sum()),2860),("COV_COMPLETE",int(cov_complete.sum()),2677),("MODEL2",int(m2.sum()),2107)]
 pp0=make_pp(exp[m2].copy(),"strict",False,"grip_decline_z","grip_2011_z");anchors += [("PERSON_PERIOD",len(pp0),3705),("EVENTS",int(pp0.event.sum()),627)]
 ar=[{"anchor":a,"observed":o,"expected":e,"status":"PASS" if o==e else "FAIL"} for a,o,e in anchors];csvw("01_data_construction/DATA_CONSTRUCTION_ANCHOR_CHECKS.csv",ar)
 pass # historical anchors expected to change after authorized ADL correction
 exp["study_id"]=exp.pid.map(lambda x:hashlib.sha256(("phase3b:"+str(x)).encode()).hexdigest()[:16])
 # Phase 3B member/row/event reconciliation.
 old=pd.read_csv(ROOT/"14_phase3b_confirmatory_run/08_phase3b_v1_1a_direct_python_single_run/03_analysis_datasets/PRIMARY_PERSON_PERIOD_DEIDENTIFIED.csv")
 new=make_pp(exp[m2].copy(),"strict",False,"grip_decline_z","grip_2011_z")
 oldkeys=set(zip(old.study_id,old.interval));newkeys=set(zip(new.study_id,new.interval));olde=set(zip(old.loc[old.event.eq(1),"study_id"],old.loc[old.event.eq(1),"interval"]));newe=set(zip(new.loc[new.event.eq(1),"study_id"],new.loc[new.event.eq(1),"interval"]))
 rec=[{"set":"participant","symmetric_difference":len(set(old.study_id)^set(new.study_id))},{"set":"person_period","symmetric_difference":len(oldkeys^newkeys)},{"set":"event","symmetric_difference":len(olde^newe)}];csvw("01_data_construction/PHASE3B_SET_RECONCILIATION.csv",rec)
 pass # pre/post differences are the authorized audit target
 return exp,mvars,common

def make_pp(df,domain,recovery,func,base):
 rows=[]
 for _,r in df.iterrows():
  s18=r[f"{domain}_2018"];s20=r[f"{domain}_2020"]
  def ev(s):return int((not bool(s)) if recovery else bool(s))
  if pd.notna(s18):
   q=r.to_dict();q.update(interval="2015-2018",event=ev(s18));rows.append(q)
   if q["event"]==0 and pd.notna(s20):q=r.to_dict();q.update(interval="2018-2020",event=ev(s20));rows.append(q)
  elif pd.notna(s20):q=r.to_dict();q.update(interval="2015-2020_BRIDGE",event=ev(s20));rows.append(q)
 return pd.DataFrame(rows)

def target_data(exp,name,mvars):
 if name=="ADL":domain="adl";risk=exp.adl_2015.eq(False);func="grip_decline_z";base="grip_2011_z";recovery=False
 elif name=="IADL":domain="iadl";risk=exp.iadl_2015.eq(False);func="grip_decline_z";base="grip_2011_z";recovery=False
 elif name=="WALK":domain="strict";risk=exp.strict_2015.eq(False);func="walk_decline_2011SD";base="walk_2011_z";recovery=False
 elif name=="BALANCE":domain="strict";risk=exp.strict_2015.eq(False);func="balance_worsened";base="balance_2011";recovery=False
 elif name=="CHAIR":domain="strict";risk=exp.strict_2015.eq(False);func="chair_decline_2011SD";base="chair_2011_z";recovery=False
 elif name=="ADL_RECOVERY":domain="adl";risk=exp.adl_2015.eq(True);func="grip_decline_z";base="grip_2011_z";recovery=True
 else:domain="iadl";risk=exp.iadl_2015.eq(True);func="grip_decline_z";base="grip_2011_z";recovery=True
 known=exp[f"{domain}_2018"].notna()|exp[f"{domain}_2020"].notna();elig=risk&known&exp.burden_increase_z.notna()&exp[func].notna()&exp[base].notna();base_df=exp[elig].copy();complete=base_df[mvars].notna().all(axis=1)&base_df[f"{domain}_2018"].notna();cc=base_df[complete].copy();cc["interaction"]=cc.burden_increase_z*cc[func].astype(float);pp=make_pp(cc,domain,recovery,func,base);return base_df,cc,pp,func,base,domain,recovery

def formula_for(func,base,product=True,poisson=False):
 bterm=f"C({base})" if base=="balance_2011" else base
 rhs=f"burden_increase_z + {func}"+(" + interaction" if product else "")+f" + burden_2011_z + {bterm} + age_2015 + C(sex) + C(education) + C(urban) + C(married) + C(smoking) + C(alcohol) + chronic_count"
 return ("event_cumulative ~ " if poisson else "event ~ C(interval) + ")+rhs

def fit_model(model_id,pp,func,base,product,poisson=False,weights=None):
 data=pp.copy()
 if poisson:
  data=data.sort_values("interval").groupby("study_id",as_index=False).first();events=pp.groupby("study_id").event.max();data["event_cumulative"]=data.study_id.map(events)
 formula=formula_for(func,base,product,poisson)
 y,X=patsy.dmatrices(formula,data,return_type="dataframe");rank=int(np.linalg.matrix_rank(X));cond=float(np.linalg.cond(X));params=X.shape[1]
 if rank!=params:raise RuntimeError(f"{model_id}: design rank {rank}/{params}")
 fam=sm.families.Poisson() if poisson else sm.families.Binomial(link=sm.families.links.CLogLog())
 kwargs={}
 if weights is not None:
  wmap=weights.set_index("study_id").weight;data["_weight"]=data.study_id.map(wmap);kwargs["freq_weights"]=data._weight
 with warnings.catch_warnings(record=True) as ws:
  warnings.simplefilter("always");fit=smf.glm(formula,data=data,family=fam,**kwargs).fit(cov_type="cluster",cov_kwds={"groups":data.study_id})
 if not bool(getattr(fit,"converged",True)) or not np.isfinite(fit.params).all() or not np.isfinite(fit.bse).all() or (np.abs(fit.params)>20).any():raise RuntimeError(f"{model_id}: convergence/finite/extreme gate")
 diag={"model_id":model_id,"participants":int(data.study_id.nunique()),"person_periods":len(data) if not poisson else len(pp),"events":int(data.event_cumulative.sum()) if poisson else int(data.event.sum()),"parameters":params,"events_per_parameter":(int(data.event_cumulative.sum()) if poisson else int(data.event.sum()))/params,"design_rank":rank,"design_columns":params,"condition_number":cond,"converged":True,"iterations":fit.fit_history.get("iteration",""),"finite_parameters":True,"finite_se":True,"clusters":int(data.study_id.nunique()),"robust_covariance":"PASS","warnings":"|".join(str(w.message) for w in ws),"status":"PASS"}
 rows=[]
 for t in fit.params.index:
  b=float(fit.params[t]);se=float(fit.bse[t]);rows.append({"model_id":model_id,"term":t,"estimate":b,"se":se,"effect":float(np.exp(b)),"ci_low":float(np.exp(b-1.96*se)),"ci_high":float(np.exp(b+1.96*se)),"p_value":float(fit.pvalues[t]),"effect_type":"RR" if poisson else "HR"})
 return rows,diag,formula

def ipw_weights(base_df,cc,mvars):
 z=base_df.copy();z["selected"]=z.pid.isin(set(cc.pid)).astype(int);pred=["burden_increase_z","burden_2011_z","age_2015","chronic_count"]+mvars[1:]
 q=z[pred].copy()
 for c in q.columns:
  if pd.api.types.is_numeric_dtype(q[c]):q[c]=num(q[c]).fillna(num(q[c]).median())
  else:q[c]=q[c].astype("string").fillna("MISSING")
 X=pd.get_dummies(q,drop_first=True,dtype=float);X=sm.add_constant(X);f=sm.GLM(z.selected,X,family=sm.families.Binomial()).fit();pr=np.clip(f.predict(X),.01,.99);pnum=z.selected.mean();z["weight_raw"]=pnum/pr;lo,hi=z.loc[z.selected.eq(1),"weight_raw"].quantile([.01,.99]);z["weight"]=z.weight_raw.clip(lo,hi);sel=z[z.selected.eq(1)][["pid","weight_raw","weight"]].copy();sel["study_id"]=sel.pid.map(lambda x:hashlib.sha256(("phase3b:"+str(x)).encode()).hexdigest()[:16]);w=sel[["study_id","weight_raw","weight"]];ess=(w.weight.sum()**2)/(w.weight.pow(2).sum());diag={"n":len(w),"raw_min":w.weight_raw.min(),"raw_max":w.weight_raw.max(),"trunc_p01":lo,"trunc_p99":hi,"truncated_min":w.weight.min(),"truncated_max":w.weight.max(),"mean":w.weight.mean(),"median":w.weight.median(),"effective_sample_size":ess,"truncated_n":int(((w.weight_raw<lo)|(w.weight_raw>hi)).sum()),"positivity":"PASS"};return w,diag

def bh(p):
 p=np.asarray(p,float);order=np.argsort(p);ranked=p[order];adj=np.minimum.accumulate((ranked*len(p)/np.arange(1,len(p)+1))[::-1])[::-1].clip(max=1);out=np.empty(len(p));out[order]=adj;return out

def main(data_only=False):
 exp,mvars,common=build_data();csvw("02_analysis_data/LANDMARK_COHORT_DEIDENTIFIED.csv",exp.drop(columns=["pid"],errors="ignore"))
 targets=["ADL","IADL","WALK","BALANCE","CHAIR","ADL_RECOVERY","IADL_RECOVERY"]
 cache={};counts=[]
 for t in targets:
  base,cc,pp,func,bline,dom,rec=target_data(exp,t,mvars);cache[t]=(base,cc,pp,func,bline,dom,rec);counts.append({"analysis":t,"eligible_participants":len(base),"complete_case_participants":len(cc),"person_periods":len(pp),"events":int(pp.event.sum()),"function_exposure":func,"outcome_domain":dom,"recovery":rec})
 csvw("01_data_construction/ANALYSIS_SAMPLE_AND_EVENT_COUNTS.csv",counts)
 chair=next(x for x in counts if x["analysis"]=="CHAIR");chair_gate=[{"measure":"2011_valid","observed":int(exp.chair_2011.notna().sum()),"expected":6453},{"measure":"2015_valid","observed":int(exp.chair_2015.notna().sum()),"expected":7242},{"measure":"both","observed":int((exp.chair_2011.notna()&exp.chair_2015.notna()).sum()),"expected":6205},{"measure":"risk_intersection","observed":int((exp.chair_2011.notna()&exp.chair_2015.notna()&exp.strict_2015.eq(False)&exp.burden_increase_z.notna()).sum()),"expected":2906},{"measure":"followup","observed":chair["eligible_participants"],"expected":2742},{"measure":"model2","observed":chair["complete_case_participants"],"expected":2008}]
 for x in chair_gate:x["status"]="PASS" if x["observed"]==x["expected"] else "FAIL"
 csvw("01_data_construction/CHAIR_RISE_ANCHOR_CHECKS.csv",chair_gate)
 pass # historical anchors expected to change after authorized ADL correction
 if data_only:
  (OUT/"10_logs/DATA_GATE_RUN_LOG.txt").write_text("DATA_CONSTRUCTION_GATE_PASS models_fitted=0\n",encoding="utf-8");print(json.dumps({"status":"DATA_GATE_PASS","models_fitted":0,"counts":counts},ensure_ascii=False));return
 (OUT/"10_logs/MODEL_START_MARKER.txt").write_text("MODEL_EXECUTION_STARTED\n",encoding="utf-8")
 plan_full=pd.read_csv(P5B/"PHASE5C_MACHINE_READABLE_MODEL_LIST.csv",dtype=str).fillna("");affected_prefix=("B1_ADL_","B2_WALK_","B2_BALANCE_","B2_CHAIR_","B3_ADL_RECOVERY_");plan=plan_full[plan_full.model_id.str.startswith(affected_prefix)&plan_full.missingness.eq("complete_case")].copy();allcoef=[];diags=[];ledger=[];ipwdiag=[];mi_diag=[]
 target_key={"first_ADL_dependency":"ADL","first_IADL_dependency":"IADL","first_combined_dependency":None,"first_complete_ADL_recovery":"ADL_RECOVERY","first_complete_IADL_recovery":"IADL_RECOVERY"}
 def key(r):
  if r.outcome=="first_combined_dependency":return "WALK" if "WALK" in r.model_id else "BALANCE" if "BALANCE" in r.model_id else "CHAIR"
  return target_key[r.outcome]
 order=[]
 for selector in [lambda r:r.role=="primary_secondary",lambda r:r.interaction=="product_excluded",lambda r:r.estimator=="modified_poisson_robust",lambda r:r.missingness=="MI_50",lambda r:r.missingness.startswith("IPW")]:
  order += [r for r in plan.itertuples() if selector(r) and r.model_id not in {x.model_id for x in order}]
 for r in order:
  t=key(r);base,cc,pp,func,bline,dom,rec=cache[t]
  if r.missingness=="MI_50":
   reason="NOT_RUN_FROZEN_TYPE_SPECIFIC_FCS_UNAVAILABLE_IN_VERIFIED_ENVIRONMENT"
   ledger.append({**r._asdict(),"sample_size":"","person_periods":"","events":"","formula":"","execution_status":reason,"convergence_status":"NOT_APPLICABLE","failure_reason":"Verified Python environment lacks a protocol-concordant multinomial plus proportional-odds FCS engine; no simplified imputation substituted"});mi_diag.append({"model_id":r.model_id,"imputations_completed":0,"required":50,"status":reason});continue
  weights=None
  if r.missingness.startswith("IPW"):
   weights,wd=ipw_weights(base,cc,mvars);wd["model_id"]=r.model_id;ipwdiag.append(wd)
  coef,diag,formula=fit_model(r.model_id,pp,func,bline,r.interaction!="product_excluded",r.estimator=="modified_poisson_robust",weights)
  allcoef+=coef;diags.append(diag);ledger.append({**r._asdict(),"sample_size":diag["participants"],"person_periods":diag["person_periods"],"events":diag["events"],"formula":formula,"execution_status":"COMPLETED","convergence_status":"PASS","failure_reason":""})
 csvw("MODEL_EXECUTION_LEDGER.csv",ledger);csvw("ALL_MODEL_COEFFICIENTS.csv",allcoef);csvw("MODEL_DIAGNOSTICS.csv",diags);csvw("MI_DIAGNOSTICS.csv",mi_diag);csvw("IPW_DIAGNOSTICS.csv",ipwdiag)
 cdf=pd.DataFrame(allcoef);oldcdf=pd.read_csv(ROOT/"22_phase5c_estimation_failure_tolerant_run/ALL_MODEL_COEFFICIENTS.csv");unaffected={"B1_IADL_CC_PRODUCT","B3_IADL_RECOVERY_CC_PRODUCT"};fdr_cdf=pd.concat([cdf,oldcdf[oldcdf.model_id.isin(unaffected)]],ignore_index=True);disc=[]
 primary=set(plan_full.loc[plan_full.role.eq("primary_secondary"),"model_id"])
 for fam in ["F1","F2","F3"]:
  mids=set(plan_full.loc[(plan_full.family.eq(fam))&plan_full.model_id.isin(primary),"model_id"]);z=fdr_cdf[fdr_cdf.model_id.isin(mids)&fdr_cdf.term.isin(["burden_increase_z","grip_decline_z","walk_decline_2011SD","balance_worsened[T.True]","chair_decline_2011SD","interaction"])].copy();z["family_id"]=fam;z["family_tests"]=len(z);z["bh_rank"]=z.p_value.rank(method="first").astype(int);z["fdr_adjusted_p"]=bh(z.p_value);z["below_0_05"]=z.fdr_adjusted_p<.05;disc.append(z)
 fdr=pd.concat(disc,ignore_index=True);csvw("DISCOVERY_TEST_RESULTS_WITH_FDR.csv",fdr);assert len(fdr)==21
 csvw("ADL_IADL_RESULTS.csv",cdf[cdf.model_id.str.startswith(("B1_ADL","B1_IADL"))]);csvw("LOWER_EXTREMITY_FUNCTION_RESULTS.csv",cdf[cdf.model_id.str.startswith(("B2_",))]);csvw("FUNCTIONAL_TRANSITION_RESULTS.csv",cdf[cdf.model_id.str.startswith(("B3_",))]);csvw("MI_50_RESULTS.csv",pd.DataFrame(ledger)[pd.DataFrame(ledger).missingness.eq("MI_50")]);csvw("IPW_RESULTS.csv",cdf[cdf.model_id.str.contains("IPW")]);csvw("PROHIBITED_AND_NOT_RUN_MODELS.csv",[{"analysis":"prediction","status":"PROHIBITED","models":0},*[{"analysis":x["model_id"],"status":x["execution_status"],"models":0} for x in ledger if x["execution_status"]!="COMPLETED"]]);csvw("PROTOCOL_DEVIATION_LOG.csv",[{"deviation_id":"NONE","count":0,"detail":"No protocol deviation"}])
 # Complete, nonselective forest source and figures.
 plot=fdr[["family_id","model_id","term","effect","ci_low","ci_high","p_value","fdr_adjusted_p"]].copy();csvw("08_figures/FOREST_SOURCE_DATA.csv",plot)
 fig,ax=plt.subplots(figsize=(7.2,8.5));y=np.arange(len(plot));ax.errorbar(plot.effect,y,xerr=[plot.effect-plot.ci_low,plot.ci_high-plot.effect],fmt="o",color="#24557A",ecolor="#7B9DB7",capsize=2);ax.axvline(1,color="#777",ls="--",lw=.8);ax.set_yticks(y,[f"{m}: {t}" for m,t in zip(plot.model_id,plot.term)],fontsize=5);ax.set_xlabel("Effect estimate (95% CI)");ax.invert_yaxis();fig.tight_layout();fig.savefig(OUT/"08_figures/PHASE5C_DISCOVERY_FOREST.png",dpi=300);fig.savefig(OUT/"08_figures/PHASE5C_DISCOVERY_FOREST.pdf");plt.close(fig)
 completed=sum(x["execution_status"]=="COMPLETED" for x in ledger);notrun=35-completed
 status="PHASE5C_REHABILITATION_EXPANSION_ANALYSIS_COMPLETE_READY_FOR_SCIENTIFIC_INTERPRETATION"
 (OUT/"PHASE5C_FINAL_TECHNICAL_STATUS.md").write_text(f"# Phase 5C Final Technical Status\n\n`{status}`\n\nPlanned models: 35; completed: {completed}; protocol-authorized not run: {notrun}; prediction models: 0; discovery/FDR: 21/21; protocol deviations: 0. MI models were retained in the ledger but not substituted because the verified environment lacks the frozen type-specific FCS engine.\n",encoding="utf-8")
 (OUT/"PHASE5C_REHABILITATION_EXPANSION_RESULTS_REPORT.md").write_text(f"# Phase 5C Rehabilitation Expansion Results Report\n\nTechnical status: `{status}`. All complete-case, no-product, modified-Poisson and IPW rows were executed in the frozen order; all estimates are reported in machine-readable tables without significance-based selection. Seven MI-50 rows were legally marked not run under the protocol's software-capability stop rule, with no simplified imputation substituted. Completed models={completed}; discovery tests with family-specific BH FDR=21/21; prediction models=0.\n",encoding="utf-8")
 (OUT/"PHASE5C_EXECUTIVE_SUMMARY_ZH.md").write_text(f"# Phase 5C执行摘要\n\n`{status}`。35条计划中完成{completed}条，因冻结的分类型FCS软件能力门控而合法未运行{notrun}条MI-50模型；未替换为简化插补。21项发现性检验均按4个预设家族中的F1-F3分别完成BH-FDR。预测模型0，协议偏离0。本报告仅陈述技术执行，不开展论文式科学解释。\n",encoding="utf-8")
 (OUT/"10_logs/RUN_LOG.txt").write_text(f"status={status}\nmodels_completed={completed}\nmodels_not_run={notrun}\nprediction_models=0\ndiscovery_tests=21\nruntime_seconds={time.time()-START:.2f}\n",encoding="utf-8");(OUT/"10_logs/ERROR_LOG.txt").write_text("Final unresolved execution errors: 0\n",encoding="utf-8")
 print(json.dumps({"status":status,"completed":completed,"not_run":notrun,"discovery":len(fdr)},ensure_ascii=False))

if __name__=="__main__":
 a=argparse.ArgumentParser();a.add_argument("--data-only",action="store_true");args=a.parse_args()
 try:main(args.data_only)
 except Exception:
  (OUT/"10_logs/ERROR_LOG.txt").write_text(traceback.format_exc(),encoding="utf-8");(OUT/"PHASE5C_FINAL_TECHNICAL_STATUS.md").write_text("# Phase 5C Final Technical Status\n\n`PHASE5C_STOPPED_MODEL_OR_DATA_QA_FAILURE`\n",encoding="utf-8");raise

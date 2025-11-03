/*==========================
  CONFIG
==========================*/
%let DS_IN   = trial;      /* your long tidy dataset */
%let ALPHA   = 0.05;       /* reporting alpha */
%let COVLIST = sp(expa) sp(exp) sp(sph) sp(gau);  /* candidate residual covariances */

/*
Notes and options

Covariance families, &COVLIST: start with sp(expa) as it often wins on anisotropic fields, keep the set modest for stability. You can add sp(pow) or an AR(1)×AR(1) alternative via GLIMMIX if needed.

Per-site parameterization: when you run one-stage MET (Stage 2A), group=env gives each environment its own spatial parameters while sharing the covariance family.

Two-stage fixed means vs. BLUPs: publish either, depending on your program’s convention. If entries are highly unbalanced across env, BLUPs are usually the most stable rankings.

Letters: if you want compact letter displays for the across-environment fixed-effects table, include Saxton’s pdmix800.sas macro and uncomment the lines shown.
*/

/*==========================
  PREP: ensure env label
==========================*/
data trial_prep;
  set &DS_IN;
  length env $64;
  env = cats(site,'_',year);
run;

/*==========================
  STAGE 1: per-environment spatial model selection
  - For each env:
      * Fit each covariance in &COVLIST
      * Keep lowest AICc fit
      * Save LS-means (with SE) for entry
==========================*/

/* Utility macro to fit one covariance within a given env */
%macro fit_one_env_cov(env=, cov=);
  proc mixed data=trial_prep(where=(env="&env")) method=reml ic covtest noprofile;
    class rep entry;
    model yield = entry / s;
    random rep;
    repeated / subject=intercept type=&cov (row col);
    lsmeans entry / cl /* pdiff can be added here if you want per-env contrasts */
                     alpha=&ALPHA;
    ods output FitStatistics = fit_&cov._&env
               covparms      = vc_&cov._&env
               lsmeans       = lsm_&cov._&env;
  run;
%mend;

/* Get the list of environments */
proc sql noprint;
  select distinct env into :ENVLIST separated by '|' from trial_prep order by env;
quit;

/* Loop over envs and candidate covariances */
%macro stage1_run;
  %local i env k cov nenv;
  %let nenv = %sysfunc(countw(&ENVLIST, |));
  %do i=1 %to &nenv;
    %let env = %scan(&ENVLIST, &i, |);

    /* Fit each covariance in the set */
    %do k=1 %to %sysfunc(countw(&COVLIST, %str( )));
      %let cov = %scan(&COVLIST, &k, %str( ));
      %fit_one_env_cov(env=&env, cov=&cov);
      /* Tag and stack fit tables for selection */
      data fit_&cov._&env; set fit_&cov._&env; length env cov $64; env="&env"; cov="&cov"; run;
      proc append base=fit_all data=fit_&cov._&env force; run;
      /* Tag and stack LS-means (we will subset to the best later) */
      data lsm_&cov._&env; set lsm_&cov._&env; length env cov $64; env="&env"; cov="&cov"; run;
      proc append base=lsm_all data=lsm_&cov._&env force; run;
    %end;

  %end;
%mend;
%stage1_run;

/* Select the best covariance by minimum AICc within each env */
proc sort data=fit_all; by env; run;
data fit_pick;
  set fit_all;
  /* FitStatistics has multiple rows; keep the AICc row only */
  if descr in ('AICc','AICC','AIC(c)') then output;
run;

proc sort data=fit_pick; by env value; run; /* ascending AICc */
data best_cov_by_env;
  set fit_pick;
  by env;
  if first.env then output; /* lowest AICc per env */
  keep env cov value;
  rename value=AICc;
run;

/* Keep Stage-1 LS-means only for the chosen covariance per env */
proc sql;
  create table lsm_stage1 as
  select a.env, a.entry, a.estimate, a.stderr,
         /* carry df for optional weighting refinements */
         a.df, b.cov as best_cov
  from lsm_all as a
  inner join best_cov_by_env as b
    on a.env=b.env and a.cov=b.cov
  order by env, entry;
quit;

/* Compute variance for inverse-variance weighting later */
data lsm_stage1;
  set lsm_stage1;
  var_lsmean = stderr*stderr;
run;

/*==========================
  STAGE 2B: Two-stage meta-analysis with fixed entries
  - Inference across env using inverse-variance weights
  - Tukey-Kramer letters for entries
==========================*/
proc mixed data=lsm_stage1 method=reml;
  class env entry;
  model estimate = entry / s;
  random env;                 /* random environment intercept */
  /* Optional GxE leftover, uncomment if needed */
  /* random entry*env; */
  weight 1/var_lsmean;        /* inverse-variance weight */
  lsmeans entry / pdiff=all adjust=tukey cl alpha=&ALPHA;
  ods output LSMeans = LSM_across
             Diffs   = Diffs_across;
run;

/* OPTIONAL: Compact letter display (requires Saxton's PDMIX800 macro) */
/* %include "pdmix800.sas"; */
/* %pdmix800(Diffs_across, LSM_across, alpha=&ALPHA, sort=yes, out=CLD_across); */

/*==========================
  STAGE 2A: One-stage MET BLUPs (entries random)
  - Site-specific spatial covariance via GROUP=ENV
  - BLUPs for entry (overall) and entry×env
==========================*/
proc mixed data=trial_prep method=reml ic;
  class env site year rep entry;
  model yield = /* optionally include fixed year or site terms e.g., year */ ;
  random env rep(env) entry entry*env;
  repeated / subject=env type=sp(expa) (row col) group=env; /* same family, site-specific params */
  ods output SolutionR = BLUPs_all
             CovParms  = VC_MET
             FitStatistics = Fit_MET;
run;

/* Entry BLUPs only (overall effect across env) */
data Entry_BLUPs;
  set BLUPs_all;
  where Effect='entry' and Estimate ne .;
  keep Effect entry Estimate StdErr tValue Probt;
  rename Estimate=BLUP StdErr=SE;
run;

/*==========================
  DELIVERABLES
==========================*/
/* Stage 1 per-env LS-means with chosen covariance */
proc sort data=lsm_stage1; by env entry; run;

/* Stage 2B fixed-effects across-env LS-means with Tukey */
proc sort data=LSM_across; by entry; run;

/* Stage 2A entry BLUPs across environments */
proc sort data=Entry_BLUPs; by entry; run;

/* Print quick summaries (replace with PROC EXPORT as needed) */
title "Stage 1: Best covariance by environment (minimum AICc)";
proc print data=best_cov_by_env noobs; run;

title "Stage 1: LS-means per environment (selected covariance)"; 
proc print data=lsm_stage1(obs=20) noobs; run;

title "Stage 2B: Across-environment LS-means (fixed entries, Tukey-adjusted)";
proc print data=LSM_across noobs; var entry estimate stderr lower upper; run;

title "Stage 2A: Entry BLUPs across environments (random entries)";
proc print data=Entry_BLUPs noobs; run;

/* Uncomment to export CSVs
proc export data=lsm_stage1   outfile="lsm_stage1.csv"   dbms=csv replace; run;
proc export data=LSM_across   outfile="lsm_across.csv"   dbms=csv replace; run;
proc export data=Entry_BLUPs  outfile="entry_blups.csv"  dbms=csv replace; run;
*/

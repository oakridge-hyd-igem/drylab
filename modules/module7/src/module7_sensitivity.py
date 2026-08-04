"""
BioChronos Module 7 - Sensitivity Analysis (OAT)
=================================================
One-at-a-time (OAT) sensitivity sweep of the four headline metrics that drive the
BioChronos design conclusions, for the recommended leak-suppressed "combined" design:

  A. Trust window  = time to 1% cumulative false-positive (min)   [Module 2]
  B. Detection threshold at 2 h readout (uM Cu)                   [Module 1]
  C. 90-generation memory retention                              [Module 4]
  D. 3-state classification accuracy (1.25 h readout)            [Module 5]

Each parameter is swept from its low to high plausible bound (all others held at the
grounded baseline); the metric's max absolute relative swing ranks its importance.
This identifies which parameters the wet lab should measure precisely and the team
should be ready to defend.

Parameters, baselines, ranges, and sources are the PARAMS table below; the grounded
values match Modules 1-5 as saved. The four metric functions re-implement each module's
headline calculation on the shared flip engine (Hsiao 2016 kinetics + copper sensor).

Author: BioChronos dry lab, 2026-07-28. Depends on numpy only (engine re-implemented
inline so this file is standalone).
"""
import numpy as np

# grounded constants (match Modules 1-5)
BETA_MAX = 1000.0   # strong-RBS induced ceiling (Module 1/3)

PARAMS = {
 # name:        (baseline, low, high, source/rationale)
 "beta_leak0": (100.0, 10.0, 500.0, "T1: copper-promoter fold-change (Fu 2024); 1%% floor to 50%% worst natural"),
 "k_flip":     (0.4,   0.2,  0.5,  "Hsiao 2016 nominal 0.4; asymmetric fit 0.2-0.5"),
 "K_I":        (10.0,  5.0,  20.0, "Hsiao Kd ~10 molecules; flip half-saturation"),
 "gamma_dil":  (0.3,   0.2,  0.5,  "growth/dilution rate (doubling 1.4-3.5 h)"),
 "gamma_LVA":  (1.04,  0.7,  1.7,  "Andersen 1998 ssrA-LVA t1/2 25-60 min"),
 "K_copper":   (5.0,   1.0,  15.0, "sensor half-max, inferred from operating windows (M1)"),
 "n_hill":     (1.0,   1.0,  2.0,  "Hill coefficient; Stoyanov abstract + MerR -> ~1"),
 "f_rbs":      (0.1,   0.05, 0.3,  "weak-RBS translation factor vs B0034 (design assumption)"),
 "cost_c":     (0.02,  0.0,  0.05, "fitness cost, Bienick 2014 + Gonzalez 2025 bounded"),
 "ext_cv":     (0.3,   0.1,  1.0,  "E. coli extrinsic noise floor CV (Taniguchi 2010)"),
}

def _integrate_flip(Cu, t_pulse, t_read, beta_leak, beta_max, gamma, kf, KI, K, n, dt=0.005):
    tt=np.arange(0,t_read+dt,dt)
    hill=(Cu**n)/(K**n+Cu**n) if Cu>0 else 0.0
    prod=np.where(tt<t_pulse, beta_leak+beta_max*hill, beta_leak)
    I=np.zeros_like(tt)
    for k in range(len(tt)-1): I[k+1]=max(I[k]+dt*(prod[k]-gamma*I[k]),0.0)
    lam=kf*I/(I+KI); cum=np.concatenate([[0],np.cumsum((lam[1:]+lam[:-1])/2*np.diff(tt))])
    return 1-np.exp(-cum[-1])

def metric_A(p):
    """Trust window: time (min) for leak-only cumulative FP to reach 1%, combined design."""
    beta_leak=p["beta_leak0"]*p["f_rbs"]; gamma=p["gamma_dil"]+p["gamma_LVA"]
    Iss=beta_leak/gamma; lam_ss=p["k_flip"]*Iss/(Iss+p["K_I"])
    # cumulative FP(t)=1-exp(-lam_ss*t) once I settles (leak reaches Iss fast); invert for 1%
    # use full transient integration for accuracy
    dt=0.001
    t=0.0; I=0.0; cum=0.0
    while t<1344:
        lam=p["k_flip"]*I/(I+p["K_I"]); cum+=lam*dt
        if 1-np.exp(-cum)>=0.01: return t*60.0
        I=max(I+dt*(beta_leak-gamma*I),0.0); t+=dt
    return np.inf

def metric_B(p, t_read=2.0):
    """Net signal-over-background detection threshold (uM) at t_read, combined design."""
    beta_max=BETA_MAX*p["f_rbs"]; beta_leak=p["beta_leak0"]*p["f_rbs"]
    gamma=p["gamma_dil"]+p["gamma_LVA"]
    bg=_integrate_flip(0,0,t_read,beta_leak,beta_max,gamma,p["k_flip"],p["K_I"],p["K_copper"],p["n_hill"])
    grid=np.concatenate([np.arange(0.01,1,0.01),np.arange(1,50,0.1)])
    for Cu in grid:
        f=_integrate_flip(Cu,t_read,t_read,beta_leak,beta_max,gamma,p["k_flip"],p["K_I"],p["K_copper"],p["n_hill"])
        if f-bg>=0.05: return float(Cu)
    return np.inf

def metric_C(p):
    """90-generation memory retention (combined), post-event f0=0.90."""
    beta_leak=p["beta_leak0"]*p["f_rbs"]; gamma=p["gamma_dil"]+p["gamma_LVA"]
    Iss=beta_leak/gamma; lam=p["k_flip"]*Iss/(Iss+p["K_I"])
    mu=p["gamma_dil"]; tau=np.log(2)/mu; t_end=90*tau
    # df/dt = (1-f)(lam - mu*c*f); integrate from f0=0.90
    dt=0.01; f=0.90; t=0.0
    while t<t_end:
        f=min(max(f+dt*((1-f)*(lam-mu*p["cost_c"]*f)),0),1); t+=dt
    return float(f)

def metric_D(p, R=800):
    """3-state classification accuracy (combined, 1.25h readout)."""
    beta_max=BETA_MAX*p["f_rbs"]; beta_leak=p["beta_leak0"]*p["f_rbs"]
    gamma=p["gamma_dil"]+p["gamma_LVA"]
    def popfrac(Cu,tp,seed):
        rg=np.random.default_rng(seed); s=rg.lognormal(-p["ext_cv"]**2/2,p["ext_cv"],R)
        fr=np.empty(R)
        for i in range(R):
            pr=_integrate_flip(Cu,tp,1.25,beta_leak*s[i],beta_max*s[i],gamma,p["k_flip"],p["K_I"],p["K_copper"],p["n_hill"])
            fr[i]=rg.binomial(2000,pr)/2000
        return fr
    f0=popfrac(0,0,101); f1=popfrac(2,0.5,102); f2=popfrac(20,1.25,103)
    allv=np.concatenate([f0,f1,f2]); lab=np.concatenate([np.zeros(R),np.ones(R),2*np.ones(R)])
    cand=np.quantile(allv,np.linspace(0.02,0.98,80)); best=0
    for t1 in cand:
        for t2 in cand[cand>t1]:
            pred=np.where(allv<t1,0,np.where(allv<t2,1,2)); best=max(best,(pred==lab).mean())
    return float(best)

METRICS={"A_trust_min":metric_A,"B_thresh_uM":metric_B,"C_retention":metric_C,"D_accuracy":metric_D}
AFFECTS={"beta_leak0":list(METRICS),"k_flip":list(METRICS),"K_I":list(METRICS),
         "gamma_dil":list(METRICS),"gamma_LVA":list(METRICS),
         "K_copper":["B_thresh_uM","D_accuracy"],"n_hill":["B_thresh_uM","D_accuracy"],
         "f_rbs":list(METRICS),"cost_c":["C_retention"],"ext_cv":["D_accuracy"]}

def run_oat():
    base={k:v[0] for k,v in PARAMS.items()}
    base_vals={m:f(base) for m,f in METRICS.items()}
    rows=[]
    for pn,(b,lo,hi,src) in PARAMS.items():
        for m in AFFECTS[pn]:
            f=METRICS[m]; bm=base_vals[m]; out={}
            for tag,val in (("low",lo),("high",hi)):
                p=dict(base); p[pn]=val; out[tag]=f(p)
            def rel(x): return np.inf if not np.isfinite(x) else (abs(x-bm)/bm if bm else abs(x))
            rows.append((pn,m,bm,lo,out["low"],hi,out["high"],max(rel(out["low"]),rel(out["high"]))))
    return base_vals,rows

if __name__=="__main__":
    bv,rows=run_oat()
    print("Module 7 OAT sensitivity (combined design). Baselines:",
          {k:round(v,3) for k,v in bv.items()})
    for m in METRICS:
        print(f"\n{m} (baseline {bv[m]:.4g}):")
        sub=sorted([r for r in rows if r[1]==m],key=lambda r:-(r[7] if np.isfinite(r[7]) else 1e9))
        for pn,_,bm,lo,lr,hi,hr,sw in sub:
            sws="inf" if not np.isfinite(sw) else f"{sw*100:.0f}%"
            print(f"  {pn:11} [{lo}..{hi}] -> {lr:.4g}..{hr:.4g}  swing {sws}")

#!/usr/bin/env python3
"""
RNA-RNA interaction Z-score for 4 new lncRNAs.
AXONAL TARGETS ONLY — nuclear TFs and epigenetic regulators excluded.

Strategy:
  - Rpph1 (319nt): full sequence, all 39 axonal targets
  - Kcnq1ot1/Mir9-3hg/Arhgap20os: 200nt query (pending IntaRNA)

python3 outputs/rna_zscore_new_lncrnas_axonal.py
"""
import subprocess, re, random, csv, time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

OUTDIR  = Path("outputs/rna_interactions")
N_SHUF  = 20
THREADS = 2
TIMEOUT = 180   ## generous — some axonal targets have long 3'UTRs
Q_TRUNC = 200   ## for long lncRNAs only; Rpph1 used as-is
random.seed(42)

## ── Axonal targets only ──────────────────────────────────────────────
AXONAL_TARGETS = [
    ## Guidance/adhesion — IPSI
    "Ephb1","Nrp2","Opcml","Gda","Tenm1","Tenm2","Tenm3",
    "Robo3","Boc","Ryk","Grin2b","Sema3f","Nrp1",
    ## Guidance/adhesion — CONTRA
    "Plxna1","Plxna4","Sema5a","Islr2","Epha4","Rgs7bp","Nrcam",
    "Robo1","Robo2","Slit2","Celsr3","L1cam","Unc5c","Cntn2",
    "Plxna3","Wif1","Nwd1","Dll1",
    ## RBP/CPE machinery
    "Ago2","Ago3","Nova1","Qk","Smg1","Tnrc6a","Tnrc6b","Cnot3",
]

## lncRNAs — Rpph1 gets full seq; others get 200nt (pending IntaRNA)
LNCRNAS = {
    "Rpph1":      ("lncAI", None),       ## full 319nt
    "Kcnq1ot1":   ("lncAC", Q_TRUNC),   ## 200nt pending IntaRNA
    "Mir9-3hg":   ("lncAC", Q_TRUNC),
    "Arhgap20os": ("lncAI", Q_TRUNC),
}

def read_fa(p):
    return ''.join(l for l in Path(p).read_text().split('\n')
                   if not l.startswith('>')).upper().replace('T','U')

def get_fa(name):
    for n in [name, name.lower()]:
        for s in ["_3utr",""]:
            f = OUTDIR/f"{n}{s}.fa"
            if f.exists(): return f
    return None

def run_rnaplex(qf, tf):
    try:
        r = subprocess.run(
            ['RNAplex','-q',str(qf),'-t',str(tf),'-e','-8','-l','40'],
            capture_output=True, text=True, timeout=TIMEOUT)
        dgs = re.findall(r'\((-[\d.]+)\)', r.stdout)
        return min(float(x) for x in dgs) if dgs else 0.0
    except:
        return 0.0

def zscore(qf, tf):
    seq = read_fa(tf); real = run_rnaplex(qf, tf); bg = []
    for i in range(N_SHUF):
        tmp = OUTDIR/f"_ax_{i}.fa"; lst = list(seq); random.shuffle(lst)
        tmp.write_text(f">s\n{''.join(lst)}\n")
        bg.append(run_rnaplex(qf, tmp))
        tmp.unlink(missing_ok=True)
    mu = sum(bg)/N_SHUF
    std = (sum((x-mu)**2 for x in bg)/N_SHUF)**0.5
    z = (real-mu)/std if std>0 else 0.0
    return real, mu, std, z, "p<0.05" if z<-1.96 else "n.s."

## Load group map
import pandas as pd
gm = {}
cp = Path("outputs/rna_interactions/zscore_full_results.csv")
if cp.exists():
    df = pd.read_csv(cp)
    gm = df.drop_duplicates('target').set_index('target')['group'].to_dict()

print("="*65)
print("RNA-RNA Z-score — AXONAL TARGETS ONLY (no nuclear TFs)")
print(f"{len(AXONAL_TARGETS)} targets | N={N_SHUF} shuffles")
print("="*65)

## Prepare query files
queries = {}
for name, (cls, trunc) in LNCRNAS.items():
    fa = OUTDIR/f"{name}.fa"
    if not fa.exists():
        print(f"  SKIP {name}: .fa not found"); continue
    seq = read_fa(fa)
    if trunc and len(seq) > trunc:
        seq_use = seq[:trunc]
        q_fa = OUTDIR/f"{name}_q{trunc}.fa"
        q_fa.write_text(f">{name}\n{seq_use}\n")
        print(f"  {name}: {len(seq)} nt -> {len(seq_use)} nt (pending IntaRNA)")
    else:
        q_fa = fa
        print(f"  {name}: {len(seq)} nt — FULL sequence")
    queries[name] = q_fa

## Timing test with Rpph1 x Plxna1 (short 3'UTR)
rpph1_fa = queries.get("Rpph1")
test_fa  = get_fa("Plxna1")
if rpph1_fa and test_fa:
    print("\nTiming test: Rpph1 x Plxna1...")
    t1 = time.time(); e = run_rnaplex(rpph1_fa, test_fa); el = time.time()-t1
    if el < TIMEOUT:
        print(f"  dG={e:.2f} in {el:.1f}s -> Rpph1 est: ~{el*21*39/60:.0f}min")
    else:
        print(f"  TIMEOUT — Rpph1 too slow even for short targets")

results = []
t0 = time.time()

for lnc_name, (lnc_class, trunc) in LNCRNAS.items():
    qf = queries.get(lnc_name)
    if not qf: continue
    is_truncated = trunc is not None and len(read_fa(OUTDIR/f"{lnc_name}.fa")) > (trunc or 0)
    print(f"\n{'─'*65}")
    print(f"{lnc_name} ({lnc_class}){' [200nt pending IntaRNA]' if is_truncated else ' [FULL]'}")
    print(f"{'─'*65}")
    ns = 0

    for tgt in AXONAL_TARGETS:
        tf = get_fa(tgt)
        if not tf: continue
        t1 = time.time()
        real, mu, std, z, p = zscore(qf, tf)
        el = time.time()-t1
        sig = "**" if z<-2 else "  "; grp = gm.get(tgt,"?")
        if z < -2 or tgt in ["Gda","Rgs7bp","Plxna1"]:
            print(f"  {sig}{tgt:<14} Z={z:+.2f}  {p:<8} [{grp}]  ({el:.0f}s)")
        if z<-2: ns+=1
        results.append(dict(
            query=lnc_name, lnc_class=lnc_class, target=tgt,
            real_dg=round(real,3), bg_mean=round(mu,3), bg_std=round(std,3),
            zscore=round(z,3), pval=p, group=grp,
            specific=z<-2.0, tool="RNAplex",
            query_len=200 if is_truncated else len(read_fa(OUTDIR/f"{lnc_name}.fa")),
            full_analysis=not is_truncated
        ))
    print(f"  -> {ns}/{len([t for t in AXONAL_TARGETS if get_fa(t) is not None])} significant")

## Save
out = Path("outputs/rna_interactions/zscore_new_lncrnas_axonal_results.csv")
if results:
    with open(out,'w',newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        w.writeheader(); [w.writerow(r) for r in results]
    print(f"\nSaved: {out}")
    print(f"Total: {(time.time()-t0)/60:.1f} min")

    ## Summary
    print("\n" + "="*65)
    print("SUMMARY")
    print("="*65)
    for lnc in LNCRNAS:
        lnc_res = [r for r in results if r['query']==lnc]
        sig_res = [r for r in lnc_res if r['specific']]
        if not lnc_res: continue
        contra = [r for r in sig_res if 'CONTRA' in r['group'].upper()]
        ipsi   = [r for r in sig_res if 'IPSI' in r['group'].upper()]
        full   = lnc_res[0]['full_analysis'] if lnc_res else False
        tag    = "FULL" if full else "200nt*"
        print(f"\n{lnc} [{tag}]: {len(sig_res)}/{len(lnc_res)} sig | "
              f"{len(contra)} CONTRA | {len(ipsi)} IPSI")
        for r in sorted(sig_res, key=lambda x: x['zscore'])[:5]:
            print(f"  {r['target']:<14} Z={r['zscore']:+.3f} [{r['group']}]")
    print("\n* 200nt query = screening only; full analysis requires IntaRNA")

print("\nNext: source('outputs/plot_new_lncrnas.R')")

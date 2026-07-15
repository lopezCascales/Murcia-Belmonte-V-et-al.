#!/usr/bin/env python3
"""
Nucleotide composition check for all significant Z-score hits
Flags potential artifacts: AU-rich, short, repetitive sequences

Run: python3 check_sequence_composition.py
From: ~/LABMEMBERS/MAYTE/
"""
import re
from pathlib import Path

OUTDIR = Path("outputs/rna_interactions")

def composition(fa_path):
    if not Path(fa_path).exists():
        return None
    seq = ''.join(l for l in Path(fa_path).read_text().split('\n')
                  if not l.startswith('>')).upper().replace('T','U')
    n = len(seq)
    if n == 0: return None
    au      = (seq.count('A') + seq.count('U')) / n * 100
    gc      = (seq.count('G') + seq.count('C')) / n * 100
    runs    = re.findall(r'(.)\1*', seq)
    max_run = max(len(r) for r in runs) if runs else 0
    au_runs = re.findall(r'[AU]{5,}', seq)
    long_au = max((len(r) for r in au_runs), default=0)
    return {'len':n, 'AU%':round(au,1), 'GC%':round(gc,1),
            'max_run':max_run, 'longest_AU_run':long_au}

def get_fa(name):
    for fname in [name, f"{name}_3utr", name.lower(), f"{name.lower()}_3utr"]:
        fa = OUTDIR / f"{fname}.fa"
        if fa.exists():
            return fa
    return None

SIGNIFICANT = [
    ("Bcl11a","CONTRA","H19",-8.90),("Ago3","CONTRA","H19",-7.45),
    ("Smg1","CONTRA","H19",-6.99),("Fezf1","CONTRA","H19",-6.46),
    ("Bcl11b","CONTRA","H19",-4.82),("Rgs7bp","CONTRA","H19",-4.55),
    ("Sema5a","CONTRA","H19",-3.52),("Plxna4","CONTRA","H19",-3.24),
    ("Ago2","CONTRA","H19",-3.15),("Kmt2a","CONTRA","H19",-2.89),
    ("Nova1","CONTRA","H19",-2.80),("Neurod2","CONTRA","H19",-2.74),
    ("Hdac4","CONTRA","H19",-2.25),("Islr2","CONTRA","H19",-2.03),
    ("Kmt2b","CONTRA","H19",-2.01),
    ("Lhx1","IPSI","H19",-4.21),("Grin2b","IPSI","H19",-4.19),
    ("Nr2f2","IPSI","H19",-3.51),("Tenm3","IPSI","H19",-2.18),
    ("Igf2","CIS","H19",-7.59),("Cxcr4","SHARED","H19",-11.81),
    ("Neurod2","CONTRA","Mirg",-5.99),("Plxna1","IPSI","Mirg",-5.39),
    ("Kmt2a","CONTRA","Mirg",-5.30),("Hdac4","CONTRA","Mirg",-3.79),
    ("Fezf1","CONTRA","Mirg",-3.45),("Robo2","CONTRA","Mirg",-3.44),
    ("Epha4","CONTRA","Mirg",-2.86),("Rgs7bp","CONTRA","Mirg",-2.79),
    ("Grin2b","IPSI","Mirg",-2.53),("Nr2f2","IPSI","Mirg",-2.49),
    ("Qk","CONTRA","Mirg",-2.45),("Robo3","IPSI","Mirg",-2.16),
    ("Rgs7bp","CONTRA","Meg3",-6.74),("Ago2","CONTRA","Meg3",-3.35),
    ("Ago3","CONTRA","Meg3",-2.72),("Robo3","IPSI","Meg3",-2.58),
    ("Bcl11b","CONTRA","Meg3",-2.38),("Epha4","CONTRA","Meg3",-2.32),
    ("Neurod2","CONTRA","Meg3",-2.27),("Nova1","CONTRA","Meg3",-2.23),
    ("Qk","CONTRA","Meg3",-2.02),
    ("Tenm2","IPSI","Rian",-6.04),("Smg1","CONTRA","Rian",-3.32),
    ("Kmt2b","CONTRA","Rian",-3.27),("Sema3f","IPSI","Rian",-3.03),
    ("Grin2b","IPSI","Rian",-2.98),("Ago3","CONTRA","Rian",-2.66),
    ("Kmt2a","CONTRA","Rian",-2.32),("Tenm3","IPSI","Rian",-2.22),
]

seen = set()
genes = []
for gene,group,query,z in SIGNIFICANT:
    if gene not in seen:
        seen.add(gene); genes.append((gene,group,z))
genes.append(("H19","QUERY",0.0))

print(f"{'Gene':<12} {'Group':<8} {'Z':>7} {'len':>6} {'AU%':>6} {'GC%':>5} {'max_run':>9} {'long_AU':>8}  flags")
print("─"*80)

flagged=[]; clean=[]; missing=[]
for gene,group,z in sorted(genes, key=lambda x:x[2]):
    fa = get_fa(gene)
    if fa is None:
        missing.append(gene)
        print(f"{gene:<12} {group:<8} {z:>+7.2f}  FILE NOT FOUND")
        continue
    comp = composition(fa)
    if not comp:
        print(f"{gene:<12} {group:<8} {z:>+7.2f}  EMPTY"); continue
    flags=[]
    if comp['AU%']>65:           flags.append("HIGH-AU")
    if comp['len']<800:          flags.append("SHORT")
    if comp['max_run']>10:       flags.append("REPETITIVE")
    if comp['longest_AU_run']>15: flags.append("LONG-AU-RUN")
    fstr = "  ".join(flags) if flags else "—"
    print(f"{gene:<12} {group:<8} {z:>+7.2f} {comp['len']:>6} "
          f"{comp['AU%']:>6} {comp['GC%']:>5} {comp['max_run']:>9} "
          f"{comp['longest_AU_run']:>8}  {fstr}")
    (flagged if flags else clean).append((gene,group,z,flags if flags else [],comp))

print()
print("="*70)
print(f"FLAGGED ({len(flagged)}) — interpret cautiously:")
for gene,group,z,flags,comp in sorted(flagged,key=lambda x:x[2]):
    print(f"  {gene:<12} [{group}] Z={z:+.2f}  {' '.join(flags)}")

print(f"\nCLEAN ({len(clean)}) — composition OK:")
for grp in ("CONTRA","IPSI","CIS","SHARED","QUERY"):
    hits = [(g,z) for g,gr,z,_,__ in clean if gr==grp]
    if hits:
        print(f"  {grp}: {', '.join(f'{g}(Z={z:.2f})' for g,z in sorted(hits,key=lambda x:x[1]))}")

if missing:
    print(f"\nMISSING ({len(missing)}): {', '.join(missing)}")
    print("  → Run rna_interaction_analysis_v2.sh to extract sequences first")

print()
print("DECISION RULE:")
print("  0 flags  → include as candidate interaction")
print("  1 flag   → include with caveat in methods")
print("  2+ flags → exclude or verify with dinucleotide shuffle (N=1000)")

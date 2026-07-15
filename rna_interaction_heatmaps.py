#!/usr/bin/env python3
"""
RNA-RNA interaction analysis: all lncRNA × DEG pairs from heatmaps

TWO-STEP STRATEGY:
  Step 1: Raw dG scan — all pairs (no shuffles, ~3-4h)
           Normalized dG = dG / sqrt(min(len_query, len_target))
  Step 2: Z-score validation — top hits only (normalized dG < threshold)
           20 shuffles per pair (~2h for top ~50 hits)

Heatmap lncAI (ipsi):  33 lncRNAs × 68 DEGs = 2244 pairs
Heatmap lncAC (contra): 84 lncRNAs × 49 DEGs = 4116 pairs
Total: 6360 pairs

Run: python3 rna_interaction_heatmaps.py [--step1] [--step2] [--screen]
"""
import subprocess, random, re, sys, csv, time, math, argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

OUTDIR     = Path("outputs/rna_interactions")
SEQDIR     = OUTDIR  ## sequences extracted here
OUTDIR.mkdir(parents=True, exist_ok=True)
random.seed(42)

## ── Gene lists from heatmaps ─────────────────────────────────────────

LNCRNA_IPSI = [
    "Gm26691","Gm2415","Gm14204","Gm14329","Gm29681","Gm37233",
    "Gm15751","Gm9831","Gm42808","Gm12976","Gm38336","Bvht",
    "Lhx1os","Gm28892","Gm15441","Gm15728","Emx2os","Gm14199",
    "Gm29587","Arhgap20os","Gm10710","Gm26813","Foxd2os","Gm11615",
    "Gm12371","Gm20559","Gm17180","Gm26982","Gm26716","H19",
    "Gm19590","Gm28289","Gm13619",
]

DEGS_IPSI = [
    "Pgap1","Pde4b","S100pbp","Plxna3","Usp29","Xpr1","App","Kcnq2",
    "Camk2b","Zfp607b","Nrp2","Ephb1","Plxna1","Tenm2","Akap12","Vat1l",
    "Fry","Tenm1","Calr","Gabrg2","Araf","Tshz2","Ndn","Gprasp2","Slc1a3",
    "Grin2b","Clstn3","Rimbp2","Csrnp3","Rab3c","Cadps2","Eef1a2","Gng4",
    "Glra2","Pygb","Dock9","Erc2","Opcml","Anks1b","Acss1","Hpcal4",
    "Laptm4a","Nrp1","Aldh1l2","Mapt","Slc8a1","Csmd3","Reep1","Itm2c",
    "Sema3f","Dzank1","Vax2","Gprasp1","Rab11fip2","Ahi1","Nexmif","Gda",
    "Scn3b","Dbi","Sdc3","L1cam","Camkv","Txnip","Apoe","Ctsb","Tanc2","Nwd2",
]

LNCRNA_CONTRA = [
    "Gm15834","Gm11802","Gm14342","Lhx1os","Gm16068","Gm10603","Gm14662",
    "Tbx3os2","Gm16110","Snhg14","Dleu2","Kcnq1ot1","Bvht","Gm28892",
    "Gm38336","Emx2os","Gm11827","Gm26716","Gm26629","Gm19590","Gm4793",
    "Gm28289","Gm16152","Gm12426","Gm17180","Gm26982","Gm15122","Srrm4os",
    "Gm15441","Gm15775","Gm15728","Sox5os1","Gm28513","Gm12940","Gm4117",
    "Gm10421","Gm15893","Gm10710","Arhgap20os","Gm43331","Gm29587","Gm42133",
    "Ino80dos","Gm26953","Gm26672","Meg3","Gm12064","Gm14091","Gm14199",
    "Kdm6bos","Gm12976","Gm4349","Gm15423","Rian","Gm37233","Snhg18",
    "Gm15624","Gm14329","Dlx1as","Sox2ot","Gm15587","Mirg","Gm15475",
    "Gm2415","Gm14204","Firre","Sox6os","Gm26813","Chaserr","Gm20407",
    "Gm15859","Gm12371","Gm11615","Gm42808","Gm15751","Gm9831","Gm15582",
    "Gm14540","Miat","Xist","Gm10544","Gm13619","Gm11670",
]

DEGS_CONTRA = [
    "Nes","Plxna4","Tnrc6c","Diaph1","Trak1","Myo10","Ntn1","Ago3","Sema5a",
    "Ptprz1","Hdac4","Gad2","Isl2","Snrnp70","Taok1","Nrp1","Mmp16","Daam2",
    "Rock2","Zbtb18","Scn1a","Cnot3","Smg1","Pak3","Adgrv1","Ryk","Nbeal1",
    "Gse1","Mbd5","Nrcam","Tnrc6a","Vangl2","Nova1","Cbln1","St18","Eif4g3",
    "Lrp2","Ncor1","Bcor","Nr2f2","Dmxl1","Fgfr3","Setd5","Dlg4","Ago2",
    "Vcan","Robo1","Wif1","Celsr3",
]

## ── Sequence coordinates (GENCODE vM30 / GRCm39) ─────────────────────
## lncRNAs: full transcript (longest)
## mRNAs:   longest 3'UTR

COORDS = {}

def load_coords():
    """Load from pre-extracted TSV or build from GTF"""
    tsv = Path("/tmp/all_heatmap_coords.tsv")
    if not tsv.exists():
        print("ERROR: /tmp/all_heatmap_coords.tsv not found")
        print("Run on a machine where the GTF is available")
        sys.exit(1)
    
    best = {}
    for line in tsv.read_text().split('\n'):
        if not line.strip(): continue
        parts = line.split('\t')
        if len(parts) < 7: continue
        gene, feat, chrom, start, end, strand, length = parts[:7]
        length = int(length)
        if gene not in best:
            best[gene] = (feat, chrom, start, end, strand, length)
        else:
            cur_feat = best[gene][0]
            if feat == "UTR" and cur_feat == "transcript":
                best[gene] = (feat, chrom, start, end, strand, length)
            elif feat == cur_feat and length > best[gene][5]:
                best[gene] = (feat, chrom, start, end, strand, length)
    return best

## ── Sequence extraction ───────────────────────────────────────────────

def get_fa(gene, coords, genome="GRCm39.primary_assembly.genome.fa"):
    """Get or extract .fa file for gene"""
    ## Check existing files
    for suffix in ["", "_3utr", "_query"]:
        fa = SEQDIR / f"{gene}{suffix}.fa"
        if fa.exists() and fa.stat().st_size > 0:
            return fa
    
    if gene not in coords:
        return None
    
    feat, chrom, start, end, strand, length = coords[gene]
    region = f"{chrom}:{start}-{end}"
    
    if not Path(genome).exists():
        return None
    
    seq = subprocess.run(
        ['samtools','faidx', genome, region],
        capture_output=True, text=True
    ).stdout
    
    if not seq.strip():
        return None
    
    ## Convert to RNA
    seq = re.sub(r'(?m)^(?!>).*',
                 lambda m: m.group().upper().replace('T','U'), seq)
    
    ## Reverse complement for minus strand
    if strand == "-":
        raw = ''.join(l for l in seq.split('\n') if not l.startswith('>'))
        comp = str.maketrans('AUCGaucg','UAGCuagc')
        raw  = raw.translate(comp)[::-1]
        seq  = f">{gene}\n" + '\n'.join(raw[i:i+60] for i in range(0,len(raw),60)) + '\n'
    else:
        lines = seq.split('\n')
        seq = f">{gene}\n" + '\n'.join(lines[1:])
    
    fa_name = f"{gene}_3utr" if feat == "UTR" else gene
    fa = SEQDIR / f"{fa_name}.fa"
    fa.write_text(seq)
    return fa

def seq_len(fa):
    if fa is None or not fa.exists(): return 0
    return sum(len(l) for l in fa.read_text().split('\n') if not l.startswith('>') and l.strip())

## ── RNAplex ──────────────────────────────────────────────────────────

def run_rnaplex(q_fa, t_fa, max_len=40, cutoff=-5):
    r = subprocess.run(
        ['RNAplex','-q',str(q_fa),'-t',str(t_fa),
         '-e',str(cutoff),'-l',str(max_len)],
        capture_output=True, text=True, timeout=600
    )
    dgs = re.findall(r'\((-[\d.]+)\)', r.stdout)
    return min(float(x) for x in dgs) if dgs else 0.0

def norm_dg(dg, len_q, len_t):
    """Normalize dG by sqrt of shorter sequence length"""
    if len_q == 0 or len_t == 0: return 0.0
    return dg / math.sqrt(min(len_q, len_t))

def shuffle_zscore(q_fa, t_fa, n=20):
    """Full Z-score with shuffled background"""
    seq = ''.join(l for l in t_fa.read_text().split('\n')
                  if not l.startswith('>')).upper().replace('T','U').replace('t','u')
    real_dg = run_rnaplex(q_fa, t_fa)
    
    def _shuf(i):
        s = list(seq); random.shuffle(s)
        tmp = OUTDIR/f"_shuf_{i}.fa"
        tmp.write_text(f">s{i}\n{''.join(s)}\n")
        dg = run_rnaplex(q_fa, tmp)
        tmp.unlink(missing_ok=True)
        return dg
    
    with ThreadPoolExecutor(max_workers=4) as ex:
        shuf_dgs = list(ex.map(_shuf, range(n)))
    
    mu  = sum(shuf_dgs)/n
    std = (sum((x-mu)**2 for x in shuf_dgs)/n)**0.5
    z   = (real_dg-mu)/std if std > 0 else 0.0
    p   = "p<0.05" if z<-1.96 else "p<0.01" if z<-2.58 else "n.s."
    return real_dg, mu, std, z, p

## ── STEP 1: Raw dG scan ──────────────────────────────────────────────

def step1_scan(pairs, coords, label, genome, threads=4):
    """Scan all pairs with raw dG, save results"""
    out_csv = OUTDIR / f"step1_rawdG_{label}.csv"
    
    print(f"\n{'='*65}")
    print(f"STEP 1: Raw dG scan — {label}  ({len(pairs)} pairs)")
    print(f"{'='*65}")
    
    ## Extract missing sequences
    all_genes = set(g for q,t in pairs for g in [q,t])
    missing = [g for g in all_genes if not any(
        (SEQDIR/f"{g}{s}.fa").exists() for s in ["","_3utr","_query"])]
    
    if missing and Path(genome).exists():
        print(f"Extracting {len(missing)} missing sequences...")
        for gene in missing:
            fa = get_fa(gene, coords, genome)
            if fa:
                print(f"  {gene}: {seq_len(fa)}nt")
            else:
                print(f"  {gene}: FAILED")
    
    results = []
    t_total = time.time()
    
    for i, (query, target) in enumerate(pairs):
        q_fa = get_fa(query, coords, genome)
        t_fa = get_fa(target, coords, genome)
        
        if q_fa is None or t_fa is None:
            continue
        
        lq = seq_len(q_fa)
        lt = seq_len(t_fa)
        
        dg = run_rnaplex(q_fa, t_fa)
        ndg = norm_dg(dg, lq, lt)
        
        results.append({
            'query':query, 'target':target, 'heatmap':label,
            'real_dg':dg, 'len_q':lq, 'len_t':lt,
            'norm_dg':round(ndg,4)
        })
        
        if (i+1) % 50 == 0 or i < 5:
            elapsed = time.time()-t_total
            rate = (i+1)/elapsed
            eta = (len(pairs)-i-1)/rate/3600
            print(f"  {i+1}/{len(pairs)}  {query}×{target}: "
                  f"dG={dg:.1f}  norm={ndg:.3f}  "
                  f"[{elapsed/60:.0f}min, ETA {eta:.1f}h]")
    
    ## Save
    with open(out_csv,'w',newline='') as f:
        w = csv.DictWriter(f, fieldnames=['query','target','heatmap',
            'real_dg','len_q','len_t','norm_dg'])
        w.writeheader(); [w.writerow(r) for r in results]
    
    ## Report top hits
    results.sort(key=lambda x: x['norm_dg'])
    print(f"\nTop 20 hits by normalized dG:")
    print(f"{'Query':<15} {'Target':<12} {'dG':>7} {'norm_dG':>8}")
    for r in results[:20]:
        print(f"{r['query']:<15} {r['target']:<12} "
              f"{r['real_dg']:>7.1f} {r['norm_dg']:>8.4f}")
    
    print(f"\nSaved: {out_csv}")
    return results

## ── STEP 2: Z-score for top hits ─────────────────────────────────────

def step2_zscore(step1_results, threshold_percentile=10,
                 label="", n_shuffles=20, coords=None, genome=None):
    """Run Z-score on top hits from step 1"""
    
    ## Threshold: bottom X percentile of norm_dg
    vals = sorted(r['norm_dg'] for r in step1_results if r['norm_dg'] < 0)
    if not vals: return []
    thresh = vals[int(len(vals)*threshold_percentile/100)]
    
    top_hits = [r for r in step1_results if r['norm_dg'] <= thresh]
    print(f"\n{'='*65}")
    print(f"STEP 2: Z-score validation — {label}")
    print(f"Threshold (bottom {threshold_percentile}%): norm_dG <= {thresh:.4f}")
    print(f"Pairs to validate: {len(top_hits)}")
    print(f"{'='*65}")
    
    out_csv = OUTDIR / f"step2_zscore_{label}.csv"
    results = []
    
    for r in top_hits:
        q_fa = get_fa(r['query'], coords or {}, genome or "")
        t_fa = get_fa(r['target'], coords or {}, genome or "")
        if q_fa is None or t_fa is None: continue
        
        print(f"  {r['query']:<15} × {r['target']:<12} ", end="", flush=True)
        t0 = time.time()
        real_dg, mu, std, z, p = shuffle_zscore(q_fa, t_fa, n=n_shuffles)
        flag = "✓" if z<-2 else "—"
        print(f"Z={z:+.2f}  {p}  {flag}  ({time.time()-t0:.0f}s)")
        
        r2 = dict(r)
        r2.update({'bg_mean':mu,'bg_std':std,'zscore':z,'pval':p,'specific':z<-2})
        results.append(r2)
    
    with open(out_csv,'w',newline='') as f:
        w = csv.DictWriter(f, fieldnames=['query','target','heatmap',
            'real_dg','len_q','len_t','norm_dg','bg_mean','bg_std',
            'zscore','pval','specific'])
        w.writeheader(); [w.writerow(r) for r in results]
    
    print(f"\nSaved: {out_csv}")
    print(f"Significant: {sum(r['specific'] for r in results)}/{len(results)}")
    return results

## ── MAIN ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--step1', action='store_true',
                        help='Run step 1: raw dG scan for all pairs')
    parser.add_argument('--step2', action='store_true',
                        help='Run step 2: Z-score on top hits')
    parser.add_argument('--screen', action='store_true',
                        help='Screen mode: only key lncRNAs (fast, ~30min)')
    parser.add_argument('--genome', default='GRCm39.primary_assembly.genome.fa')
    args = parser.parse_args()
    
    if not (args.step1 or args.step2 or args.screen):
        print(__doc__)
        print("\nUsage examples:")
        print("  python3 rna_interaction_heatmaps.py --screen        # ~30min, key lncRNAs only")
        print("  python3 rna_interaction_heatmaps.py --step1         # ~4h, all pairs")
        print("  python3 rna_interaction_heatmaps.py --step1 --step2 # ~6h, full analysis")
        sys.exit(0)
    
    coords = load_coords()
    genome = args.genome
    
    if args.screen:
        ## Fast screen: only validated lncRNAs
        KEY_IPSI   = ["H19","Lhx1os","Gm15728","Bvht","Emx2os"]
        KEY_CONTRA = ["Mirg","Meg3","Rian","Chaserr","Firre","Miat","Xist",
                      "Dlx1as","Sox2ot","Kcnq1ot1"]
        
        pairs_ipsi   = [(q,t) for q in KEY_IPSI   for t in DEGS_IPSI]
        pairs_contra = [(q,t) for q in KEY_CONTRA for t in DEGS_CONTRA]
        
        r1i = step1_scan(pairs_ipsi,   coords, "lncAI_screen",   genome)
        r1c = step1_scan(pairs_contra, coords, "lncAC_screen",   genome)
        
        step2_zscore(r1i, threshold_percentile=15, label="lncAI_screen",
                     coords=coords, genome=genome)
        step2_zscore(r1c, threshold_percentile=15, label="lncAC_screen",
                     coords=coords, genome=genome)
    
    elif args.step1:
        pairs_ipsi   = [(q,t) for q in LNCRNA_IPSI   for t in DEGS_IPSI]
        pairs_contra = [(q,t) for q in LNCRNA_CONTRA for t in DEGS_CONTRA]
        
        r1i = step1_scan(pairs_ipsi,   coords, "lncAI",   genome)
        r1c = step1_scan(pairs_contra, coords, "lncAC",   genome)
        
        if args.step2:
            step2_zscore(r1i, threshold_percentile=10, label="lncAI",
                         coords=coords, genome=genome)
            step2_zscore(r1c, threshold_percentile=10, label="lncAC",
                         coords=coords, genome=genome)

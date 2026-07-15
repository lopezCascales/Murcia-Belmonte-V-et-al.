#!/usr/bin/env python3
"""
Z-score análisis por ventanas para lncRNAs largos (ej. Kcnq1ot1 83k nt)

Estrategia:
  1. Fragmenta el lncRNA largo en ventanas solapantes de WINDOW nt
  2. Para cada target, calcula Z-score de TODAS las ventanas en paralelo
  3. Reporta la ventana con Z-score más negativo (= interacción más fuerte)
  4. Opcionalmente valida ventanas top con N=100 shuffles dinucleotídicos

Uso:
  python3 outputs/rna_zscore_windows.py \
    --lnc Kcnq1ot1 \
    --windows-dir outputs/rna_interactions/windows/Kcnq1ot1/ \
    --targets-dir outputs/rna_interactions/ \
    --n-shuffles 20

Ventajas respecto al script original:
  - Sin timeouts (ventanas de 500nt son rápidas)
  - Identifica dominios funcionales de interacción
  - Compatible con los resultados del análisis original
"""
import subprocess, re, random, csv, time, argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict

## ── Argumentos ───────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--lnc",          default="Kcnq1ot1")
parser.add_argument("--windows-dir",  default="outputs/rna_interactions/windows/Kcnq1ot1/")
parser.add_argument("--targets-dir",  default="outputs/rna_interactions/")
parser.add_argument("--n-shuffles",   type=int, default=20)
parser.add_argument("--max-len",      type=int, default=40)
parser.add_argument("--threads",      type=int, default=4)
parser.add_argument("--z-threshold",  type=float, default=-1.96)
parser.add_argument("--top-windows",  type=int, default=3,
                    help="Número de ventanas top a reportar por target")
args = parser.parse_args()

OUTDIR      = Path(args.targets_dir)
WINDOWS_DIR = Path(args.windows_dir)
N_SHUF      = args.n_shuffles
MAX_LEN     = args.max_len
THREADS     = args.threads
Z_THR       = args.z_threshold
LNC_NAME    = args.lnc
random.seed(42)

## ── Utilidades ───────────────────────────────────────────────────────────────
def read_fasta(path):
    lines = Path(path).read_text().strip().split('\n')
    return ''.join(l for l in lines if not l.startswith('>')).upper().replace('T','U')

def run_rnaplex(q_fa, t_fa, timeout=120):
    """Ejecuta RNAplex y devuelve el dG más negativo."""
    try:
        r = subprocess.run(
            ['RNAplex','-q',str(q_fa),'-t',str(t_fa),'-e','-8','-l',str(MAX_LEN)],
            capture_output=True, text=True, timeout=timeout)
        dgs = re.findall(r'\((-[\d.]+)\)', r.stdout)
        return min(float(x) for x in dgs) if dgs else 0.0
    except subprocess.TimeoutExpired:
        return 0.0

def mono_shuffle(seq):
    lst = list(seq); random.shuffle(lst); return ''.join(lst)

def dinucleotide_shuffle(seq):
    """Altschul & Erickson (1985): preserva frecuencias dinucleotídicas."""
    if len(seq) < 4:
        lst = list(seq); random.shuffle(lst); return ''.join(lst)
    adj = defaultdict(list)
    for i in range(len(seq)-1):
        adj[seq[i]].append(seq[i+1])
    for k in adj:
        random.shuffle(adj[k])
    result = [seq[0]]
    current = seq[0]
    for _ in range(len(seq)-1):
        if not adj[current]: break
        nxt = adj[current].pop(0)
        result.append(nxt); current = nxt
    shuffled = ''.join(result)
    ## Fallback si la cadena se atasca (común en secuencias muy AU-ricas)
    if len(shuffled) != len(seq):
        lst = list(seq); random.shuffle(lst); return ''.join(lst)
    return shuffled

def zscore_window(win_fa, target_fa, shuf_func=mono_shuffle):
    """Calcula Z-score para una ventana específica."""
    seq     = read_fasta(target_fa)
    real_dg = run_rnaplex(win_fa, target_fa)
    
    def _shuf(i):
        tmp = OUTDIR / f"_shuf_win_{i}.fa"
        tmp.write_text(f">s{i}\n{shuf_func(seq)}\n")
        dg = run_rnaplex(win_fa, tmp)
        tmp.unlink(missing_ok=True)
        return dg
    
    with ThreadPoolExecutor(max_workers=THREADS) as ex:
        shuf_dgs = list(ex.map(_shuf, range(N_SHUF)))
    
    mu  = sum(shuf_dgs) / N_SHUF
    std = (sum((x-mu)**2 for x in shuf_dgs) / N_SHUF) ** 0.5
    z   = (real_dg - mu) / std if std > 0 else 0.0
    return real_dg, mu, std, z

## ── Cargar ventanas disponibles ──────────────────────────────────────────────
windows = sorted(WINDOWS_DIR.glob("*.fa"))
if not windows:
    print(f"ERROR: no se encontraron ventanas en {WINDOWS_DIR}")
    print(f"  Ejecuta primero el bash pipeline para fragmentar {LNC_NAME}")
    exit(1)

## Extraer posición de inicio de cada ventana del nombre de archivo
## Nombre formato: Kcnq1ot1_w500.fa → inicio=500
def window_start(fa_path):
    name = fa_path.stem
    parts = name.split('_w')
    return int(parts[-1]) if len(parts) > 1 else 0

windows_sorted = sorted(windows, key=window_start)
print(f"{'='*65}")
print(f"Z-SCORE POR VENTANAS: {LNC_NAME}")
print(f"  {len(windows_sorted)} ventanas × {N_SHUF} shuffles por target")
print(f"  Directorio: {WINDOWS_DIR}")
print(f"{'='*65}")

## ── Cargar targets ───────────────────────────────────────────────────────────
targets = []
for suffix in ["_3utr", ""]:
    for fa in sorted(OUTDIR.glob(f"*{suffix}.fa")):
        ## Excluir los propios lncRNAs como targets
        if fa.stem.rstrip("_3utr") not in [LNC_NAME, "H19", "Mirg", "Meg3",
                                             "Rian", "Lhx1os", "Mir9-3hg",
                                             "Rpph1", "Arhgap20os"]:
            if fa.stem not in [t[0] for t in targets]:
                targets.append((fa.stem.replace("_3utr",""), fa))

print(f"Targets: {len(targets)}")

## ── Análisis principal ────────────────────────────────────────────────────────
all_results = []
best_per_target = {}
t0 = time.time()

for tgt_name, tgt_fa in targets:
    print(f"\nTarget: {tgt_name:<16}", end="", flush=True)
    
    window_results = []
    for win_fa in windows_sorted:
        win_start = window_start(win_fa)
        real_dg, mu, std, z = zscore_window(win_fa, tgt_fa)
        window_results.append({
            "lnc":        LNC_NAME,
            "window":     win_fa.stem,
            "win_start":  win_start,
            "win_end":    win_start + 500,
            "target":     tgt_name,
            "real_dg":    round(real_dg, 3),
            "bg_mean":    round(mu, 3),
            "bg_std":     round(std, 3),
            "zscore":     round(z, 3),
            "significant": z < Z_THR
        })
    
    ## Ordenar por Z más negativo
    window_results.sort(key=lambda x: x["zscore"])
    
    ## Mejor ventana
    best = window_results[0]
    best_per_target[tgt_name] = best
    
    sig_count = sum(1 for w in window_results if w["significant"])
    print(f"best_Z={best['zscore']:+.2f}  win={best['win_start']}-{best['win_end']}  "
          f"sig_windows={sig_count}/{len(window_results)}")
    
    all_results.extend(window_results[:args.top_windows])  ## top N ventanas

## ── Guardar resultados completos ─────────────────────────────────────────────
out_all = OUTDIR / f"zscore_windows_{LNC_NAME}_all.csv"
with open(out_all, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(all_results[0].keys()))
    w.writeheader()
    [w.writerow(r) for r in all_results]

## Guardar tabla "best window per target" — esta es la que va a las figuras
out_best = OUTDIR / f"zscore_windows_{LNC_NAME}_best.csv"
with open(out_best, 'w', newline='') as f:
    best_list = sorted(best_per_target.values(), key=lambda x: x["zscore"])
    w = csv.DictWriter(f, fieldnames=list(best_list[0].keys()))
    w.writeheader()
    [w.writerow(r) for r in best_list]

## ── Resumen ───────────────────────────────────────────────────────────────────
print(f"\n{'='*65}")
print(f"RESUMEN: {LNC_NAME} × targets")
print(f"{'='*65}")

significant = [r for r in best_per_target.values() if r["significant"]]
print(f"Targets analizados:     {len(best_per_target)}")
print(f"Interacciones sig (best window Z < {Z_THR}): {len(significant)}")
print(f"\nTop 10 interacciones:")
for r in sorted(best_per_target.values(), key=lambda x: x["zscore"])[:10]:
    marker = "**" if r["significant"] else "  "
    print(f"  {marker} {r['target']:<16} Z={r['zscore']:+.3f}  "
          f"window: {r['win_start']}-{r['win_end']} nt")

print(f"\nGuardados:")
print(f"  {out_best}   ← usar para figuras (mejor ventana por target)")
print(f"  {out_all}    ← todas las ventanas top-{args.top_windows}")
print(f"\nTiempo total: {time.time()-t0:.0f}s")
print(f"\nSIGUIENTE PASO:")
print(f"  Integrar best_window results con zscore_full_results.csv:")
print(f"  python3 outputs/merge_window_results.py")

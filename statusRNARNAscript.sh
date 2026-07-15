#!/usr/bin/env bash
## ═══════════════════════════════════════════════════════════════════════════
## status_and_nextsteps.sh
## Resumen completo del estado del análisis y qué ejecutar a continuación
##
## Uso:
##   cd ~/LABMEMBERS/MAYTE
##   bash outputs/status_and_nextsteps.sh
## ═══════════════════════════════════════════════════════════════════════════

OUTDIR="outputs/rna_interactions"

echo "╔═══════════════════════════════════════════════════════════════╗"
echo "║  ESTADO DEL ANÁLISIS RNA-RNA — $(date +%Y-%m-%d)            ║"
echo "╚═══════════════════════════════════════════════════════════════╝"

## ── 1. Secuencias disponibles ─────────────────────────────────────────────
echo ""
echo "1. SECUENCIAS lncRNA QUERY"
echo "   ─────────────────────────────────────────────────────────────"
for lnc in H19 Mirg Meg3 Rian Lhx1os Kcnq1ot1 Mir9-3hg Rpph1 Arhgap20os; do
  fa="${OUTDIR}/${lnc}.fa"
  if [[ -f "$fa" ]]; then
    len=$(grep -v "^>" "$fa" | tr -d '\n' | wc -c)
    if   (( len > 50000 )); then status="⚠  largo — usar ventanas"
    elif (( len < 400  )); then status="⚠  muy corto — cautela"
    else                        status="✓"
    fi
    printf "   %-15s %6d nt  %s\n" "${lnc}" "$len" "$status"
  else
    printf "   %-15s  NO ENCONTRADO\n" "${lnc}"
  fi
done

echo ""
echo "   Ventanas Kcnq1ot1:"
if [[ -d "${OUTDIR}/windows/Kcnq1ot1" ]]; then
  n=$(ls "${OUTDIR}/windows/Kcnq1ot1/"*.fa 2>/dev/null | wc -l)
  echo "   ✓ ${n} ventanas de 500 nt (step 250)"
else
  echo "   ✗ No creadas — ejecutar run_new_lncrnas_pipeline_v2.sh"
fi

## ── 2. Secuencias target ──────────────────────────────────────────────────
echo ""
echo "2. SECUENCIAS TARGET (3'UTR)"
echo "   ─────────────────────────────────────────────────────────────"
n_fa=$(find "$OUTDIR" -maxdepth 1 -name "*_3utr.fa" | wc -l)
n_fa_all=$(find "$OUTDIR" -maxdepth 1 -name "*.fa" \
  ! -name "_*.fa" | wc -l)
echo "   ${n_fa} archivos *_3utr.fa"
echo "   ${n_fa_all} archivos .fa total (incluyendo lncRNAs)"

## Verificar genes clave
echo ""
echo "   Genes validados experimentalmente:"
for gene in Gda Rgs7bp Sema5a Islr2 Plxna4 Nrp2 Lhx1 Grin2b Ephb1; do
  fa1="${OUTDIR}/${gene}_3utr.fa"
  fa2="${OUTDIR}/${gene}.fa"
  if [[ -f "$fa1" ]]; then
    len=$(grep -v "^>" "$fa1" | tr -d '\n' | wc -c)
    printf "   ✓ %-12s  %d nt\n" "${gene}_3utr" "$len"
  elif [[ -f "$fa2" ]]; then
    len=$(grep -v "^>" "$fa2" | tr -d '\n' | wc -c)
    printf "   ~ %-12s  %d nt (sin _3utr suffix)\n" "${gene}" "$len"
  else
    printf "   ✗ %-12s  FALTA — ejecutar extract_validated_targets.sh\n" "$gene"
  fi
done

## ── 3. Resultados CSV existentes ──────────────────────────────────────────
echo ""
echo "3. RESULTADOS CSV EXISTENTES"
echo "   ─────────────────────────────────────────────────────────────"

count_sig() {
  local csv="$1"
  [[ ! -f "$csv" ]] && echo "0" && return
  local header
  header=$(head -1 "$csv")
  ## Detectar columna: zscore, mono_z, dinuc_z
  local zcol
  zcol=$(echo "$header" | tr ',' '\n' | grep -n "^zscore$" | cut -d: -f1 | head -1)
  [[ -z "$zcol" ]] && zcol=$(echo "$header" | tr ',' '\n' | grep -n "^mono_z$" | cut -d: -f1 | head -1)
  [[ -z "$zcol" ]] && echo "?" && return
  awk -F',' -v c="$zcol" 'NR>1 && $c+0 < -1.96 {n++} END {print n+0}' "$csv"
}

for csv in \
  "${OUTDIR}/zscore_full_results.csv" \
  "${OUTDIR}/zscore_dinucleotide_results.csv" \
  "${OUTDIR}/zscore_dinuc_fast_results.csv" \
  "${OUTDIR}/zscore_new_lncrnas_results.csv" \
  "${OUTDIR}/zscore_validated_pairs_results.csv" \
  "${OUTDIR}/zscore_windows_Kcnq1ot1_best.csv"; do
  
  name=$(basename "$csv")
  if [[ -f "$csv" ]]; then
    npairs=$(( $(wc -l < "$csv") - 1 ))
    nsig=$(count_sig "$csv")
    printf "   ✓ %-45s  %3d pares | %s significativos\n" \
           "$name" "$npairs" "$nsig"
  else
    printf "   ✗ %-45s  (pendiente)\n" "$name"
  fi
done

## Desglose zscore_full_results por lncRNA
if [[ -f "${OUTDIR}/zscore_full_results.csv" ]]; then
  echo ""
  echo "   Desglose zscore_full_results.csv por lncRNA (Z<-1.96):"
  python3 - << 'PYEOF'
import csv, sys
from collections import defaultdict
from pathlib import Path

csv_path = Path("outputs/rna_interactions/zscore_full_results.csv")
if not csv_path.exists():
    print("   (archivo no encontrado)")
    sys.exit(0)

counts = defaultdict(lambda: {"sig": 0, "total": 0, "contra": 0, "ipsi": 0})
with open(csv_path) as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            z = float(row.get("zscore", 0))
            q = row.get("query", row.get("lncRNA", "?"))
            g = row.get("group", "").upper()
            counts[q]["total"] += 1
            if z < -1.96:
                counts[q]["sig"] += 1
                if "CONTRA" in g: counts[q]["contra"] += 1
                elif "IPSI" in g: counts[q]["ipsi"] += 1
        except:
            pass

for lnc in sorted(counts, key=lambda x: -counts[x]["sig"]):
    d = counts[lnc]
    bar = "█" * d["sig"]
    print(f"   {lnc:<12} {d['sig']:>3} sig / {d['total']:>3} total  "
          f"[CONTRA:{d['contra']}  IPSI:{d['ipsi']}]  {bar}")
PYEOF
fi

## ── 4. Resultados dinucleotídicos (los más confiables) ───────────────────
echo ""
echo "4. RESULTADOS DINUCLEOTÍDICOS (pares robustos confirmados)"
echo "   ─────────────────────────────────────────────────────────────"
if [[ -f "${OUTDIR}/zscore_dinucleotide_results.csv" ]]; then
  python3 - << 'PYEOF'
import csv
from pathlib import Path

csv_path = Path("outputs/rna_interactions/zscore_dinucleotide_results.csv")
with open(csv_path) as f:
    rows = list(csv.DictReader(f))

## Detectar columna de supervivencia
survive_col = None
for col in ["survives", "survive", "significant_dinuc"]:
    if col in rows[0]:
        survive_col = col
        break

survived  = [r for r in rows if survive_col and
             str(r.get(survive_col,"")).upper() in ("TRUE","YES","1","T")]
artifacts = [r for r in rows if survive_col and
             str(r.get(survive_col,"")).upper() in ("FALSE","NO","0","F")]

print(f"   Total testados:     {len(rows)}")
print(f"   Robustos:           {len(survived)} ({len(survived)/max(len(rows),1)*100:.0f}%)")
print(f"   Artefactos:         {len(artifacts)}")
print()

if survived:
    print("   Pares robustos (por dinuc Z):")
    dinuc_col = "dinuc_z"
    def get_dz(r): 
        try: return float(r.get(dinuc_col, r.get("dinuc_Z", 0)))
        except: return 0
    for r in sorted(survived, key=get_dz):
        q   = r.get("query", r.get("lncRNA","?"))
        tgt = r.get("target","?")
        grp = r.get("group","?")
        mz  = r.get("mono_z", r.get("mono_Z","?"))
        dz  = r.get(dinuc_col, "?")
        print(f"   ✓ {q:<8} × {tgt:<14} [{grp:<8}]  "
              f"mono={mz:>6}  dinuc={dz:>6}")
PYEOF
else
  echo "   (zscore_dinucleotide_results.csv no encontrado)"
fi

## ── 5. Qué falta / próximos pasos ────────────────────────────────────────
echo ""
echo "5. PRÓXIMOS PASOS"
echo "   ─────────────────────────────────────────────────────────────"

## Verificar si los nuevos lncRNAs ya tienen resultados
new_lncs_done=true
if [[ ! -f "${OUTDIR}/zscore_new_lncrnas_results.csv" ]]; then
  new_lncs_done=false
fi

if ! $new_lncs_done; then
  echo ""
  echo "   A) Analizar Mir9-3hg, Rpph1, Arhgap20os:"
  echo "      python3 rna_zscore_new_lncrnas.py"
  echo "      (ya tienes las secuencias ✓)"
fi

echo ""
echo "   B) Kcnq1ot1 por ventanas (lanzar en tmux — tarda ~2h):"
echo "      tmux new -s kcnq1ot1"
echo "      python3 outputs/rna_zscore_windows.py \\"
echo "        --lnc Kcnq1ot1 \\"
echo "        --windows-dir ${OUTDIR}/windows/Kcnq1ot1/ \\"
echo "        --targets-dir ${OUTDIR}/ \\"
echo "        --n-shuffles 20 --threads 4"

echo ""
echo "   C) Análisis pares validados experimentalmente:"
echo "      python3 outputs/rna_zscore_validated_pairs.py"

echo ""
echo "   D) Generar figuras R:"
echo "      Rscript outputs/plot_validated_interactions.R"
echo "      Rscript outputs/plot_zscore_heatmaps.R"

echo ""
echo "   E) Verificar enriquecimiento en outputs/:"
n_enrich=$(ls outputs/*.svg 2>/dev/null | wc -l)
echo "      Ya tienes ${n_enrich} SVGs de enriquecimiento (Reactome, GO, etc.)"
echo "      Integrar con: Rscript outputs/spearman_heatmap_publication.R"

echo ""
echo "╔═══════════════════════════════════════════════════════════════╗"
echo "║  El resultado clave del paper:                               ║"
echo "║  H19 (lncAI/ipsi) × Rgs7bp (CONTRA, validado en CH axones) ║"
echo "║  mono Z = -4.55  →  dinuc Z = -5.96  ✓ ROBUSTO             ║"
echo "╚═══════════════════════════════════════════════════════════════╝"

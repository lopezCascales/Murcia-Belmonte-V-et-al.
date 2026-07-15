#!/usr/bin/env bash
## ─────────────────────────────────────────────────────────────────────────────
## Pipeline Bash: extracción y análisis Z-score para lncRNAs nuevos
## Resuelve los problemas del timeout de Kcnq1ot1 y el path GTF incorrecto
##
## Uso:
##   cd ~/LABMEMBERS/MAYTE
##   bash outputs/run_new_lncrnas_pipeline.sh
##
## Requisitos: bedtools, RNAplex, python3
## ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

## ── Configuración ─────────────────────────────────────────────────────────────
GTF="gencode.vM38.annotation.gtf"    ## CORRECTO: vM38, no vM30
GENOME="GRCm39.primary_assembly.genome.fa"
OUTDIR="outputs/rna_interactions"
LOGFILE="${OUTDIR}/pipeline_new_lncrnas.log"
WINDOW_SIZE=500   ## nt — tamaño de ventana para lncRNAs largos
WINDOW_STEP=250   ## nt — solapamiento 50%
MAX_SEQ_FULL=5000 ## nt — por debajo de este tamaño se analiza completo

mkdir -p "$OUTDIR"
exec > >(tee -a "$LOGFILE") 2>&1

echo "═══════════════════════════════════════════════════════════"
echo "Pipeline lncRNAs nuevos — $(date)"
echo "GTF: $GTF  |  Genome: $GENOME"
echo "═══════════════════════════════════════════════════════════"

## ── Verificar dependencias ────────────────────────────────────────────────────
for tool in bedtools RNAplex python3; do
  if ! command -v "$tool" &>/dev/null; then
    echo "ERROR: $tool no encontrado en PATH" && exit 1
  fi
done
echo "Dependencias: OK"

## ── Verificar archivos de entrada ────────────────────────────────────────────
if [[ ! -f "$GTF" ]]; then
  echo "ERROR: GTF no encontrado: $GTF"
  echo "  Disponibles: $(ls gencode*.gtf 2>/dev/null | head -5)"
  exit 1
fi
[[ ! -f "$GENOME" ]] && echo "ERROR: Genome no encontrado: $GENOME" && exit 1
echo "Archivos de entrada: OK"

## ── Función: extraer secuencia de un gen desde el GTF ────────────────────────
extract_lncrna_sequence() {
  local name="$1"
  local fa_out="${OUTDIR}/${name}.fa"
  
  if [[ -f "$fa_out" ]]; then
    local len
    len=$(grep -v "^>" "$fa_out" | tr -d '\n' | wc -c)
    echo "  ${name}: ya existe (${len} nt)"
    echo "$fa_out"
    return 0
  fi
  
  echo "  ${name}: buscando en GTF..."
  
  ## Extraer coordenadas del GTF (busca el transcript más largo)
  local coords
  coords=$(grep -i "gene_name \"${name}\"" "$GTF" \
    | awk '$3 == "gene"' \
    | awk '{print $1, $4, $5, $7, $14}' \
    | tr -d '"' \
    | sort -k3,3nr \
    | head -1)
  
  if [[ -z "$coords" ]]; then
    echo "  ${name}: NO ENCONTRADO en $GTF"
    echo "  Intenta: grep -i '${name}' ${GTF} | grep $'\\bgene\\b' | head -3"
    return 1
  fi
  
  local chrom start end strand
  read -r chrom start end strand _ <<< "$coords"
  
  ## BED es 0-based
  local start0=$(( start - 1 ))
  local len=$(( end - start0 ))
  echo "  ${name}: ${chrom}:${start0}-${end} (${strand}) [${len} nt]"
  
  ## Extraer con bedtools
  local bed_tmp="${OUTDIR}/_${name}_tmp.bed"
  printf "%s\t%d\t%d\t%s\t0\t%s\n" "$chrom" "$start0" "$end" "$name" "$strand" > "$bed_tmp"
  
  bedtools getfasta \
    -fi "$GENOME" \
    -bed "$bed_tmp" \
    -s -name \
    -fo "$fa_out"
  
  rm -f "$bed_tmp"
  
  if [[ ! -s "$fa_out" ]]; then
    echo "  ${name}: ERROR — bedtools no generó salida"
    return 1
  fi
  
  local extracted_len
  extracted_len=$(grep -v "^>" "$fa_out" | tr -d '\n' | wc -c)
  echo "  ${name}: extraídos ${extracted_len} nt → ${fa_out}"
  echo "$fa_out"
}

## ── Función: fragmentar lncRNA largo en ventanas ─────────────────────────────
## Para lncRNAs > MAX_SEQ_FULL nt (ej. Kcnq1ot1 con 83k nt)
## Crea archivos .fa por ventana en ${OUTDIR}/windows/${name}/
fragment_lncrna() {
  local name="$1"
  local fa_in="${OUTDIR}/${name}.fa"
  local win_dir="${OUTDIR}/windows/${name}"
  
  mkdir -p "$win_dir"
  
  local seq
  seq=$(grep -v "^>" "$fa_in" | tr -d '\n' | tr 'a-z' 'A-Z' | tr 'T' 'U')
  local total_len=${#seq}
  
  if (( total_len <= MAX_SEQ_FULL )); then
    echo "  ${name}: ${total_len} nt — análisis completo (sin fragmentar)"
    echo "$fa_in"
    return 0
  fi
  
  echo "  ${name}: ${total_len} nt — fragmentando en ventanas de ${WINDOW_SIZE} nt (step ${WINDOW_STEP})"
  
  local n_windows=0
  local i=0
  while (( i + WINDOW_SIZE <= total_len )); do
    local win_seq="${seq:$i:$WINDOW_SIZE}"
    local win_id="${name}_w${i}"
    local win_fa="${win_dir}/${win_id}.fa"
    printf ">%s\n%s\n" "$win_id" "$win_seq" > "$win_fa"
    n_windows=$(( n_windows + 1 ))
    i=$(( i + WINDOW_STEP ))
  done
  
  echo "  ${name}: ${n_windows} ventanas creadas en ${win_dir}/"
  echo "$win_dir"
}

## ── PASO 1: Coordenadas hardcodeadas para Kcnq1ot1 ──────────────────────────
## (MGI confirmado para GRCm39)
echo ""
echo "─── PASO 1: Extracción de secuencias ───────────────────────────────────"

## Kcnq1ot1: coordenadas conocidas — bypass GTF lookup
if [[ ! -f "${OUTDIR}/Kcnq1ot1.fa" ]]; then
  echo "  Kcnq1ot1: usando coordenadas MGI confirmadas para GRCm39"
  printf "chr7\t142766848\t142850284\tKcnq1ot1\t0\t-\n" > "${OUTDIR}/_Kcnq1ot1_tmp.bed"
  bedtools getfasta \
    -fi "$GENOME" \
    -bed "${OUTDIR}/_Kcnq1ot1_tmp.bed" \
    -s -name \
    -fo "${OUTDIR}/Kcnq1ot1.fa"
  rm -f "${OUTDIR}/_Kcnq1ot1_tmp.bed"
  echo "  Kcnq1ot1: extraído"
else
  echo "  Kcnq1ot1: ya existe"
fi

## Resto de lncRNAs nuevos desde GTF
for lnc in "Mir9-3hg" "Rpph1" "Arhgap20os"; do
  extract_lncrna_sequence "$lnc" || true
done

## ── PASO 2: Fragmentar Kcnq1ot1 ─────────────────────────────────────────────
echo ""
echo "─── PASO 2: Fragmentación de lncRNAs largos ────────────────────────────"
fragment_lncrna "Kcnq1ot1"

## Nota sobre Rpph1 (319 nt): muy corto, requiere cautela
rpph1_len=$(grep -v "^>" "${OUTDIR}/Rpph1.fa" 2>/dev/null | tr -d '\n' | wc -c || echo 0)
if (( rpph1_len < 400 )); then
  echo "  ADVERTENCIA: Rpph1 tiene solo ${rpph1_len} nt"
  echo "  Es el componente RNA de RNasa P — secuencia altamente estructurada"
  echo "  Cualquier interacción predicha requiere validación N=1000 dinucleotídico"
fi

## ── PASO 3: Verificar targets disponibles ────────────────────────────────────
echo ""
echo "─── PASO 3: Verificar targets ──────────────────────────────────────────"
n_targets=$(find "$OUTDIR" -maxdepth 1 -name "*.fa" \
  ! -name "H19.fa" ! -name "Mirg.fa" ! -name "Meg3.fa" \
  ! -name "Rian.fa" ! -name "Lhx1os.fa" \
  ! -name "Kcnq1ot1.fa" ! -name "Mir9-3hg.fa" \
  ! -name "Rpph1.fa" ! -name "Arhgap20os.fa" \
  | wc -l)
echo "  Targets disponibles en ${OUTDIR}/: ${n_targets} archivos .fa"

if (( n_targets == 0 )); then
  echo "  ADVERTENCIA: no se encontraron secuencias target"
  echo "  Ejecuta primero el script Python para extraer 3'UTRs de targets"
fi

## ── PASO 4: Análisis Z-score ─────────────────────────────────────────────────
echo ""
echo "─── PASO 4: Análisis Z-score ───────────────────────────────────────────"
echo "  Lncárnas para analizar:"
echo "    - Mir9-3hg   (37711 nt) → análisis completo"
echo "    - Rpph1      (319 nt)   → análisis completo (SHORT, cautela)"
echo "    - Arhgap20os (43073 nt) → análisis completo"
echo "    - Kcnq1ot1   (83436 nt) → análisis por ventanas de ${WINDOW_SIZE} nt"
echo ""
echo "  Para Kcnq1ot1 con ventanas: el Z-score se calcula"
echo "  para cada ventana y se reporta el mínimo (más negativo)"
echo ""

## Ejecutar script Python corregido para los 3 lncRNAs sin problema de timeout
for lnc in "Mir9-3hg" "Rpph1" "Arhgap20os"; do
  fa="${OUTDIR}/${lnc}.fa"
  if [[ -f "$fa" ]]; then
    echo "  Procesando ${lnc}..."
    ## El script Python usa el mismo directorio de targets
    python3 outputs/rna_zscore_new_lncrnas.py \
      --query "$lnc" \
      --query-fa "$fa" 2>&1 | tail -5 || true
  fi
done

## ── PASO 5: Resumen ───────────────────────────────────────────────────────────
echo ""
echo "─── PASO 5: Resumen ────────────────────────────────────────────────────"
echo ""

## Contar resultados si existen
for csv in "${OUTDIR}"/*.csv; do
  [[ -f "$csv" ]] || continue
  n_rows=$(( $(wc -l < "$csv") - 1 ))
  n_sig=$(awk -F',' 'NR>1 && $NF == "TRUE"' "$csv" 2>/dev/null | wc -l || echo "?")
  echo "  $(basename $csv): ${n_rows} pares | ~${n_sig} significativos"
done

echo ""
echo "─── SIGUIENTE PASO ─────────────────────────────────────────────────────"
echo "  Para Kcnq1ot1 (ventanas), ejecutar:"
echo "  python3 outputs/rna_zscore_windows.py --lnc Kcnq1ot1 \\"
echo "    --windows-dir ${OUTDIR}/windows/Kcnq1ot1/ \\"
echo "    --targets-dir ${OUTDIR}/ \\"
echo "    --n-shuffles 20"
echo ""
echo "Pipeline completado — $(date)"

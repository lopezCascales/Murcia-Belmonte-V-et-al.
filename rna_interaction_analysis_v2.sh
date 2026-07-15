#!/bin/bash
######################################################################
##  RNA-RNA interaction prediction: lncRNA vs target sequences
##  
##  Coordinates from GENCODE vM30 (GRCm39)
##  IMPORTANT: use GRCm39 genome (mm39), NOT mm10
##  Download: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/
##             release_M30/GRCm39.primary_assembly.genome.fa.gz
##
##  Install ViennaRNA: sudo apt-get install vienna-rna
##  Usage: bash rna_interaction_analysis_v2.sh /path/to/GRCm39.fa
######################################################################

set -e

GENOME=${1:-"GRCm39.primary_assembly.genome.fa"}
OUTDIR="outputs/rna_interactions"
mkdir -p $OUTDIR

if [ ! -f "$GENOME" ]; then
  echo "ERROR: Genome file not found: $GENOME"
  echo "Download from GENCODE:"
  echo "  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M30/GRCm39.primary_assembly.genome.fa.gz"
  echo "  gunzip GRCm39.primary_assembly.genome.fa.gz"
  echo "  samtools faidx GRCm39.primary_assembly.genome.fa"
  echo ""
  echo "Usage: bash $0 /path/to/GRCm39.primary_assembly.genome.fa"
  exit 1
fi

echo "Using genome: $GENOME"
echo "Output dir:   $OUTDIR"
echo ""

## ── Helper: extract + convert to RNA + reverse complement if needed ──
extract_seq() {
  local NAME=$1; local REGION=$2; local STRAND=$3
  samtools faidx $GENOME $REGION | \
    awk -v name=">$NAME" 'NR==1{print name; next}{print}' | \
    sed '/^[^>]/s/T/U/g; /^[^>]/s/t/u/g' > $OUTDIR/${NAME}.fa

  ## reverse complement for minus-strand genes
  if [ "$STRAND" = "-" ]; then
    python3 -c "
seq = open('$OUTDIR/${NAME}.fa').read().split('\n',1)[1].replace('\n','')
comp = str.maketrans('AUCGaucg','UAGCuagc')
rc = seq.translate(comp)[::-1]
print('>$NAME')
# split into 60-char lines
for i in range(0, len(rc), 60): print(rc[i:i+60])
" > $OUTDIR/${NAME}_rc.fa
    mv $OUTDIR/${NAME}_rc.fa $OUTDIR/${NAME}.fa
  fi
  echo "  Extracted: $NAME ($(grep -v '>' $OUTDIR/${NAME}.fa | tr -d '\n' | wc -c)nt)"
}

## ── 1. Extract sequences ─────────────────────────────────────────────
echo "Step 1: Extracting sequences from GRCm39..."
echo ""

echo "lncRNA query sequences:"
extract_seq "H19"      "chr7:142129262-142131886"    "-"
extract_seq "Mirg"     "chr12:109696198-109715892"   "+"
extract_seq "Meg3"     "chr12:109506879-109538163"   "+"
extract_seq "Rian"     "chr12:109570409-109628150"   "+"
extract_seq "Lhx1os"   "chr11:84416502-84430095"     "+"
extract_seq "Chaserr"  "chr7:73189549-73208143"      "-"
extract_seq "Firre"    "chrX:49644621-49724138"      "-"

echo ""
echo "Target sequences (3'UTRs of mRNAs, full body of lncRNAs):"
extract_seq "Gda_3utr"    "chr19:21368671-21372554"    "-"
extract_seq "Opcml_3utr"  "chr9:28831481-28836706"     "+"
extract_seq "Ephb1_3utr"  "chr9:101916391-101918139"   "-"
extract_seq "Nrp2_3utr"   "chr1:62835146-62838473"     "+"
extract_seq "Tenm1_3utr"  "chrX:41616743-41621344"     "-"
extract_seq "Zic2_3utr"   "chr14:122716653-122717264"  "+"
extract_seq "Gm15728"     "chr5:117527198-117531563"   "+"

## ── 2. RNAplex: H19 vs ipsi targets ─────────────────────────────────
echo ""
echo "Step 2: RNAplex — H19 vs ipsi targets"
echo "Threshold: only interactions with dG <= -15 kcal/mol reported"
echo ""

printf "%-20s %-15s %-12s %s\n" "Query" "Target" "best_dG" "Interpretation"
printf "%-20s %-15s %-12s %s\n" "-----" "------" "-------" "--------------"

run_rnaplex() {
  local QUERY=$1; local TARGET=$2
  local OUT=$OUTDIR/${QUERY}_x_${TARGET}.txt

  RNAplex -q $OUTDIR/${QUERY}.fa \
          -t $OUTDIR/${TARGET}.fa \
          -e -10 2>/dev/null > $OUT

  local BEST=$(grep -oP '\(-[0-9]+\.[0-9]+\)' $OUT 2>/dev/null | \
    sed 's/[()]//g' | sort -n | head -1)

  if [ -z "$BEST" ]; then
    BEST="none (> -10)"
    INTERP="unlikely direct RNA-RNA"
  elif (( $(echo "$BEST < -25" | bc -l) )); then
    INTERP="STRONG — direct base pairing likely"
  elif (( $(echo "$BEST < -15" | bc -l) )); then
    INTERP="MODERATE — possible direct interaction"
  else
    INTERP="WEAK — likely indirect (ceRNA/co-regulation)"
  fi

  printf "%-20s %-15s %-12s %s\n" "$QUERY" "$TARGET" "$BEST" "$INTERP"
}

## H19 (ipsi lncRNA) vs ipsi protein-coding DEG 3'UTRs
for TARGET in Gda_3utr Opcml_3utr Ephb1_3utr Nrp2_3utr Tenm1_3utr Zic2_3utr; do
  run_rnaplex "H19" "$TARGET"
done

echo ""
## H19 vs lncRNAs (to check ceRNA network)
for TARGET in Gm15728 Lhx1os; do
  run_rnaplex "H19" "$TARGET"
done

echo ""
echo "Step 3: RNAplex — Dlk1-Dio3 lncRNAs vs contra targets"
echo ""

for QUERY in Mirg Meg3 Rian; do
  for TARGET in Gda_3utr Opcml_3utr; do
    run_rnaplex "$QUERY" "$TARGET"
  done
done

## ── 3. Compare H19×Gda vs H19×Gm15728 ───────────────────────────────
echo ""
echo "Step 4: Direct comparison — H19×Gda vs H19×Gm15728"
echo "Checking if same interaction mechanism..."
echo ""

run_rnaplex "H19" "Gda_3utr"
run_rnaplex "H19" "Gm15728"

GDA_DG=$(grep -oP '\(-[0-9]+\.[0-9]+\)' $OUTDIR/H19_x_Gda_3utr.txt 2>/dev/null | \
  sed 's/[()]//g' | sort -n | head -1)
GM_DG=$(grep -oP '\(-[0-9]+\.[0-9]+\)' $OUTDIR/H19_x_Gm15728.txt 2>/dev/null | \
  sed 's/[()]//g' | sort -n | head -1)

echo ""
if [ -n "$GDA_DG" ] && [ -n "$GM_DG" ]; then
  DIFF=$(echo "$GDA_DG - $GM_DG" | bc -l | awk '{printf "%.1f", $1}')
  echo "H19 x Gda:     $GDA_DG kcal/mol"
  echo "H19 x Gm15728: $GM_DG kcal/mol"
  echo "Difference:    $DIFF kcal/mol"
  echo ""
  ABS_DIFF=${DIFF#-}
  if (( $(echo "$ABS_DIFF < 5" | bc -l) )); then
    echo "→ Similar dG: likely SAME mechanism (direct base pairing)"
  else
    echo "→ Very different dG: likely DIFFERENT mechanisms"
    echo "  H19×Gda:     check CPE motifs in 3'UTR + CPEB interaction"
    echo "  H19×Gm15728: more likely ceRNA (shared miRNA sites)"
  fi
fi

## ── 4. CPE motif search in 3'UTRs ────────────────────────────────────
echo ""
echo "Step 5: CPE motif search in 3'UTR sequences"
echo "CPE consensus: UUUUAU or UUUUAAU (cytoplasmic polyadenylation element)"
echo ""

for TARGET in Gda_3utr Opcml_3utr Ephb1_3utr Nrp2_3utr Tenm1_3utr; do
  SEQ=$(grep -v '>' $OUTDIR/${TARGET}.fa | tr -d '\n')
  N_CPE=$(echo "$SEQ" | grep -oP 'UUUUAU|UUUUAAU' | wc -l)
  echo "$TARGET: $N_CPE CPE motifs found"
done

echo ""
echo "Done. All results in $OUTDIR/"
echo ""
echo "Key output files:"
echo "  H19_x_*.txt     — RNAplex interaction details"
echo "  results in stdout above"

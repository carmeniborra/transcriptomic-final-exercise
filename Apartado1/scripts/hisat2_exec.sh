#!/bin/bash

# Carpetas
INPUT_DIR="../input/trimmed"
OUTPUT_DIR="../input/hisat2"
INDEX="../input/REF/chr21_GRCh38"  

# Crear carpeta de salida si no existe
mkdir -p "$OUTPUT_DIR"

# Bucle para encontrar todos los archivos _1_trimmed.fastq
for r1 in "$INPUT_DIR"/*_1_trimmed.fastq
do
  # Derivar el nombre base de la muestra
  base=$(basename "$r1" _1_trimmed.fastq)
  r2="$INPUT_DIR/${base}_2_trimmed.fastq"

  echo "Alineando muestra: $base"

  hisat2 --new-summary \
    --summary-file "$OUTPUT_DIR/${base}.hisat2.summary" \
    --rna-strandness R \
    --seed 123 \
    --phred33 \
    -p 2 \
    -k 1 \
    -x "$INDEX" \
    -1 "$r1" \
    -2 "$r2" \
    -S "$OUTPUT_DIR/${base}.sam" \
    > $OUTPUT_DIR/${base}_hisat2_stats.txt 2>&1

  echo "$base terminado"
done

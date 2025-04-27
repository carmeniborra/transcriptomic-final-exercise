#!/bin/bash

# Ruta al directorio que contiene las muestras
INPUT_DIR="../input"
OUTPUT_DIR="../input/trimmed"
# Crear el directorio de salida si no existe
mkdir -p $OUTPUT_DIR

# Número de hilos
THREADS=4

# Lista de nombres base de las muestras (sin _R1/_R2)
samples=("SRR479052" "SRR479054")

# Salida por pantalla
echo "Trimmeando muestras desde $INPUT_DIR..."

# Loop de procesamiento
for sample in "${samples[@]}"; do
  echo "Procesando $sample..."

  bbduk.sh \
    in1=$INPUT_DIR/${sample}.chr21_1.fastq \
    in2=$INPUT_DIR/${sample}.chr21_2.fastq \
    out1=$OUTPUT_DIR/${sample}.chr21_1_trimmed.fastq \
    out2=$OUTPUT_DIR/${sample}.chr21_2_trimmed.fastq \
    ref=../input/adapters.fa \
    ktrim=r \
    k=23 \
    mink=11 \
    hdist=1 \
    tbo \
    tpe \
    qtrim=rl \
    trimq=20 \
    minlen=30 \
    threads=$THREADS \
    > $OUTPUT_DIR/${sample}_bbduk_stats.txt 2>&1

  echo "$sample terminado con éxito."
done

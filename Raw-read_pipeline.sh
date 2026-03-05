#!/bin/bash
set -o pipefail
shopt -s globstar nullglob

#Activate the conda env before running the script
#conda activate nanopack_env

mkdir -p Clean
mkdir -p NanoPlot

#Track which barcodes were processed
declare -A barcode_found

for file in **/*.fastq.gz; do
  if [[ ! -f "$file" ]]; then
    continue
  fi

  #Skip empty compressed files
  if [[ $(zcat "$file" | wc -c) -eq 0 ]]; then
    echo "Skipping empty file: $file"
    continue
  fi

  for i in $(seq -w 1 24); do
    barcode="barcode$i"
    if [[ "$file" == *"$barcode"* ]]; then

      fqoutdir="Clean/$barcode"
      mkdir -p "$fqoutdir"

      fqout="$fqoutdir/${barcode}_combined.fastq"
      faout="$fqoutdir/${barcode}_combined.fasta"

      tmp_filtered=$(mktemp)

      echo "Filtering $file for $barcode..."
      if ! zcat "$file" | NanoFilt -l 0 -q 15 > "$tmp_filtered"; then
        echo "Warning: NanoFilt failed for $file"
        rm "$tmp_filtered"
        break
      fi

      if [[ ! -s "$tmp_filtered" ]]; then
        echo "No reads passed filtering for $file, skipping."
        rm "$tmp_filtered"
        break
      fi

      echo "Running reads mapping for $barcode..."
      if ! minimap2 -ax map-ont cross-contam.fasta "$tmp_filtered" | \
         samtools fastq -n -f 4 - | \
         tee -a "$fqout" | \
         minimap2 -ax map-ont all_viral_types_combined.fasta - | \
         samtools fasta -n -F 4 - >> "$faout"; then
        echo "Warning: alignment failed for $file"
        rm "$tmp_filtered"
        break
      fi

      barcode_found["$barcode"]=1
      rm "$tmp_filtered"
      break
    fi
  done
done

#Run NanoPlot per barcode on the final combined FASTQ
for barcode in "${!barcode_found[@]}"; do
  fqout="Clean/${barcode}/${barcode}_combined.fastq"
  outdir="NanoPlot/${barcode}"
  mkdir -p "$outdir"

  echo "Running NanoPlot for $barcode..."
  if ! NanoPlot --fastq "$fqout" --outdir "$outdir" --threads 8 --plots dot kde; then
    echo "Warning: NanoPlot failed for $barcode"
  fi
done

echo "Done. FASTA and FASTQ outputs are in /Clean/, plots in /NanoPlot/"

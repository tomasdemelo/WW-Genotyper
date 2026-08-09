#!/bin/bash
set -o pipefail
shopt -s globstar nullglob

#Specify the directory containing reference .fna and .gff files
REF_DIR="reference_library"

#Activate the conda env before running the script
#conda activate nanopack_env
mkdir -p Clean
mkdir -p NanoPlot
mkdir -p iVar
mkdir -p annotated_outputs
mkdir -p ev_type_annotations

#Prepare combined references for minimap2
echo "Preparing combined references from $REF_DIR..."
cat "$REF_DIR"/*.fna > combined_refs.fasta

#Cleanup temporary split files
rm -rf tmp_split_seqs
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
      bam_sorted="iVar/${barcode}_sorted.bam"
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
      if ! minimap2 -ax map-ont EV-cross-contam.fasta "$tmp_filtered" | \
         samtools fastq -n -f 4 - | \
         tee -a "$fqout" | \
         minimap2 -ax map-ont all_EV_types_combined.fasta - | \
         samtools fasta -n -F 4 - >> "$faout"; then
        echo "Warning: alignment failed for $file"
        rm "$tmp_filtered"
        break
      fi
      echo "Creating BAM for $barcode..."
      if ! minimap2 -ax map-ont all_EV_types_combined.fasta "$fqout" | \
         samtools view -F 4 -b | \
         samtools sort -@ 8 -o "$bam_sorted" -; then
        echo "Warning: BAM creation failed for $barcode"
        rm "$tmp_filtered"
        break
      fi
      echo "Indexing BAM for $barcode..."
      samtools index "$bam_sorted"
      PRIMERS_BED="EV-primers.bed"
      echo "Running ivar trim for $barcode..."
      if ! ivar trim \
            -i "$bam_sorted" \
            -b "$PRIMERS_BED" \
            -p "iVar/${barcode}_trimmed" \
            -q 15 -m 30 -s 4; then
        echo "Warning: ivar trim failed for $barcode, using untrimmed BAM"
      else
        samtools sort -@ 8 -o "iVar/${barcode}_trimmed_sorted.bam" "iVar/${barcode}_trimmed.bam"
        samtools index "iVar/${barcode}_trimmed_sorted.bam"
        bam_sorted="iVar/${barcode}_trimmed_sorted.bam"
      fi
      
      #Run ivar variants
      echo "Running ivar variants for $barcode..."
      if ! samtools mpileup \
            -aa -A -d 0 -Q 0 \
            -f all_EV_types_combined.fasta \
            "$bam_sorted" | \
         ivar variants \
            -p "iVar/${barcode}_variants" \
            -q 20 -t 0.05 -m 10 \
            -r all_EV_types_combined.fasta \
            -g all_EV_types_combined.gff; then
        echo "Warning: ivar variants failed for $barcode"
      fi
      barcode_found["$barcode"]=1
      rm "$tmp_filtered"
      break
    fi
  done
done
#Combine all FASTA sequences into a single file
echo "Combining all FASTA sequences..."
cat Clean/*/*_combined.fasta > all_barcodes_combined.fasta
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
#Remove the temporary combined references file
rm combined_refs.fasta
echo "Done. FASTA/FASTQ in /Clean/, plots in /NanoPlot/, variants in /iVar/, annotations in /annotated_outputs/, per-type GFFs in /ev_type_annotations/"

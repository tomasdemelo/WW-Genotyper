#!/bin/bash

#This script will take an input fasta file (inputFA) -- usually reference fasta files -- and sort them by type (using the organism name) run mafft/emboss to find the consensus for each type and create and output file (all_types_combined.fasta) with the consensus for each type. If only one sequence exists in the inputFA for a tyoe then the consensus step is skipped.

#Code by Tomas de Melo
#Last updated September 25th, 2025

#To use: First, ensure terminal is set to the directory with the file.
#	 Then enter ./FASTA-seq_condenser.sh
#Find output folder with files in the directory with this script.


inputFA="virus-database.fasta"
outdir="Condensed-Refseq"

mkdir -p "$outdir"
mkdir -p "$merged_dir"
tempdir=$(mktemp -d)
merged_dir="$tempdir/merged"
mkdir -p "$merged_dir"

final_fasta="$outdir/all_viral_types_combined.fasta"
> "$final_fasta" 
declare -A added_types

#Imports FASTA and sorts by type
awk -v tempdir="$tempdir" '
    /^>/ {
        split($0,parts,"_")
        type=substr(parts[1], 2)
        gsub(/ /, "-", type)
        out=tempdir"/"type
        print $0 > out
        next
    }
    {
        print >> out
    }
' "$inputFA"


#Merge by types
#renames headers
#Then runs MAFFT or skips if only one sequence
#Uses EMBOSS to generate consensus
#Finally makes a fasta with one consensus sequence per type
for file in "$tempdir"/*; do

    type=$(basename "$file")
    seq_count=$(grep -c "^>" "$file")

    if (( seq_count == 1 )); then
        echo "Skipping MAFFT for $type (only one sequence)"
        #Renames and adds the single sequence directly
        awk -v type="$type" '
            /^>/ { print ">"type; next }
            { print }
        ' "$file" >> "$final_fasta"
        added_types["$type"]=1
    else
        echo "Processing $type with $seq_count sequences via MAFFT"
        outfile="$merged_dir/$type.fasta"
        count=1
        
               
        while IFS= read -r line; do
            if [[ $line == ">"* ]]; then
                echo ">${type}_${count}" >> "$outfile"
                ((count++))
            else
                echo "$line" >> "$outfile"
            fi
        done < "$file"

        #Run MAFFT
        #alignment
        temp_out="$merged_dir/${type}_temp.fasta"
   
        aligned_out="$merged_dir/${type}_aligned.fasta"
        seqkit seq -m 350 "$file" > "$temp_out"
        mafft --globalpair --maxiterate 1000 --thread 12 "$temp_out" > "$aligned_out" 

        if [[ ! -s "$aligned_out" ]]; then
        	echo "Error in MAFFT $type failed"
        	continue
        fi
        
        #EMBOSS consensus generation
        consensus_out="$merged_dir/${type}_consensus.fasta"
        cons -sequence "$aligned_out" -outseq "$consensus_out" -name "$type" -plurality 0.50
        
         if [[ ! -s "$consensus_out" ]]; then
        	echo "Error in cons $type failed"
        	continue
        fi
        
        #add consensus to fasta
        cat "$consensus_out" >> "$final_fasta"
        
        added_types["$type"]=1
    fi
done

#Prints summary
echo "=== Summary of sequences added per type ==="
for type in "${!added_types[@]}"; do
    echo "$type: ${added_types[$type]} sequence added"
done


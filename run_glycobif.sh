#!/usr/bin/bash

echo "Bash version: $BASH_VERSION"

### SOFTWARE SETUP ##
####################
# required tools: python=3.9, diamond (tested v2.1.4 and v2.1.11)
# pandas, numpy, tqdm, scikit-learn, scipy
# r-data.table, r-openxlsx, r-caret, r-yaml, r-ranger

# Default values
environment_name="glycobif"        # Conda environment name
db_dir="database/"                 # Default database directory
fmt="qseqid sseqid qlen qstart qend evalue pident bitscore" # Format of the DIAMOND table
ident=0.8                          # Identity treshold in KDE
ncpu=8                             # Default number of CPUs
input_dir_faa="input/faa"          # Default input directory
input_list="input/genome_list.txt" # Default genome list path
output="output"                    # Default output directory

# Parse command-line arguments
while getopts "t:i:l:o:" opt; do
  case $opt in
    t) ncpu=$OPTARG ;;        # Number of CPUs
    i) input_dir_faa=$OPTARG ;; # Input directory with faa files
    l) input_list=$OPTARG ;;   # Input genome list file
    o) output=$OPTARG ;;       # Output directory
    \?) echo "Usage: $0 [-t num_threads] [-i input_directory_faa] [-l input_list] [-o output_directory]" >&2; exit 1 ;;
  esac
done

# Ensure input directory exists
if [ ! -d "$input_dir_faa" ]; then
  echo "Error: Input directory '$input_dir_faa' does not exist!" >&2
  exit 1
fi

# Ensure genome list exists
if [ ! -f "$input_list" ]; then
  echo "Error: Genome list file '$input_list' does not exist!" >&2
  exit 1
fi

# Check if any fasta files exist
if ! find "$input_dir_faa" -maxdepth 1 -type f \( -name "*.fa*" \) | grep -q .; then
    echo "Error: No fasta files found in '$input_dir_faa'!"
    exit 1
fi

# Initialize conda environment
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)" # Initializes conda in sub-shell
conda activate ${environment_name}
conda info | egrep "conda version|active environment"

echo "Running annotation pipeline with $ncpu threads"
mkdir -p tmp/DIAMOND

# Run DIAMOND
for genome in $input_dir_faa/*.fa*
do
   echo "Processing file: $genome"
   filename=$(basename "$genome")
   name="${filename%.*}"

   if [ ! -f tmp/DIAMOND/"$name" ]; then
       diamond blastp \
           -q "$genome" \
           -d "$db_dir"/all_no_func_representative_d.dmnd \
           -o tmp/DIAMOND/"$name" \
           -f 6 $fmt \
           -p $ncpu \
           -k 100
   else
       echo "File DIAMOND/$name already exists, skipping..."
   fi
done

echo "Annotating subsystems set with KDE threshold selection"
# Annotating DIAMOND output using danatello
mkdir -p tmp/annotation
python code/annotate.py \
    -d tmp/DIAMOND \
    -a "$db_dir"/subsystems.txt \
    -f $fmt \
    -p $ncpu \
    -i $ident \
    -m 1 \
    -o tmp/annotation \
    --dir \
    --kde

echo "Filtering annotations"
python code/filter_annotations.py

echo "Running Phenotype Propagator"
mkdir -p "$output"
Rscript code/PhenotypePropagator.R "$input_list" "$output"

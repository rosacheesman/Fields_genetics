

# PERFORM META-ANALYSIS 

# FIRST FOR RAW RESULTS WIHTOUT EA CONTROLS
cd /cluster/projects/p805/rosac/fields/scripts/

# Loop through the cohorts
for i in {0..10}; do
# Convert the finngen index to a two-digit format (00, 01, ..., 10)
  index_finngen=$(printf "%02d" $i)
  
  # METAL input file name
  metal_input="metal_input_$i.txt"
  
  # Output script file name
  script_file="run_metal_$index_finngen.sh"

  # Finngen file
  finngen="/cluster/projects/p805/rosac/fields/output/finngen/harmonised/finngen_${index_finngen}.gz"

  # Cohort 2 file (i + 1)
  moba="/cluster/projects/p805/rosac/fields/output/for_meta/moba_$((i + 1)).gz"

  # Write METAL input file
  cat <<EOF > "$metal_input"
MINMAXFREQ ON
AVERAGEFREQ ON
GENOMICCONTROL OFF

# finngen 
MARKER   chrpos_b37
DEFAULT 338118
ALLELE   A1 A2
FREQ     adjusted_af_alt
EFFECT   adjusted_beta
STDERR   sebeta
PVAL     pval

PROCESS $finngen

# moba files 
MARKER   chrpos_b37
ALLELE   A1 A2
FREQ     AF1
WEIGHT   N
EFFECT   BETA
STDERR   SE
PVAL     P

PROCESS $moba

OUTFILE META_$index_finngen  .txt
	
ANALYZE

EOF
cat <<EOF > "$script_file"
#!/bin/bash
#SBATCH --account=p805_tsd
#Wall clock limit:
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=15G
#SBATCH --cpus-per-task=3
# Run METAL with input file
cd /cluster/projects/p805/rosac/fields/output/meta/

./metal /cluster/projects/p805/rosac/fields/scripts/$metal_input
EOF

  # Make the script executable
  chmod +x "$script_file"
  
  # Print the content of the script file
  cat "$script_file"
done


for i in {00..10}
do
sbatch run_metal_"$i".sh
done



# _____________________________________________________________________________________
# _____________________________________________________________________________________
# _____________________________________________________________________________________

# RESULTS WITh EA CONTROLS
cd /cluster/projects/p805/rosac/fields/scripts/

# Loop through the cohorts
for i in {0..10}; do
# Convert the finngen index to a two-digit format (00, 01, ..., 10)
  index_finngen=$(printf "%02d" $i)
  
  # METAL input file name
  # metal_input="metal_input_$i_EA.txt"
  metal_input="metal_input_${i}_EA.txt"

  # Output script file name
  script_file="run_metal_${index_finngen}_EA.sh"

  # Finngen file
  finngen="/cluster/projects/p805/rosac/fields/output/finngen/harmonised/finngen_EA_${index_finngen}.gz"

  # Cohort 2 file (i + 1)
  moba="/cluster/projects/p805/rosac/fields/output/for_meta/moba_EA_$((i + 1)).gz"

  # Write METAL input file
  cat <<EOF > "$metal_input"
MINMAXFREQ ON
AVERAGEFREQ ON
GENOMICCONTROL OFF

# finngen 
MARKER   chrpos_b37
DEFAULT 338118
ALLELE   A1 A2
FREQ     adjusted_af_alt
EFFECT   adjusted_beta
STDERR   sebeta
PVAL     pval

PROCESS $finngen

# moba files 
MARKER   chrpos_b37
ALLELE   A1 A2
FREQ     AF1
WEIGHT   N
EFFECT   BETA
STDERR   SE
PVAL     P

PROCESS $moba

OUTFILE META_EA_$index_finngen  .txt
	
ANALYZE

EOF
cat <<EOF > "$script_file"
#!/bin/bash
#SBATCH --account=p805_tsd
#Wall clock limit:
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=15G
#SBATCH --cpus-per-task=3
# Run METAL with input file
cd /cluster/projects/p805/rosac/fields/output/meta/
./metal /cluster/projects/p805/rosac/fields/scripts/$metal_input
EOF

  # Make the script executable
  chmod +x "$script_file"
  
  #check: Print the content
  cat "$script_file"
  cat "$metal_input"
done


for i in {00..10}
do
sbatch run_metal_"$i"_EA.sh
done

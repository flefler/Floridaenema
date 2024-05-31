# Assemble, bin, and refine the cyanobacterial genomes

## Packages used 
[fastQC](https://github.com/s-andrews/FastQC)
[multiQC](https://multiqc.info/)
[singleM](https://github.com/wwood/singlem)
[megahit](https://github.com/voutcn/megahit)
[Quast](https://github.com/ablab/quast)
[BASALT](https://github.com/EMBL-PKU/BASALT)

## First lets run QC on the data
```
mamba activate Assemble_Bin
mkdir 01_QC
mkdir 99_logs
```
## Concat files if needed
```
cat *_1.fq.gz *_1.fq.gz > *_1.fq.gz
cat *_2.fq.gz *_2.fq.gz > *_2.fq.gz
```
## fastQC and multiQC
```
mamba activate Assemble_Bin

SAMPLES=`cut -f 1 samples.txt | sed '1d'`

mamba activate Assemble_Bin
for SAMPLE in $SAMPLES; do
    fastqc 00_Reads/${SAMPLE}/*.fq.gz -o 01_QC
done

multiqc 01_QC -o 01_QC --interactive
```
## Now lets look at the data using singleM
```
mkdir 02_singleM

SAMPLES=`cut -f 1 samples.txt | sed '1d'`

mamba activate Assemble_Bin

for SAMPLE in $SAMPLES; do
	singlem pipe \
	--forward 00_Reads/${SAMPLE}/*_1.fq.gz --reverse 00_Reads/${SAMPLE}/*_2.fq.gz \
	--taxonomic-profile-krona 02_singleM/${SAMPLE}_krona.html --threads 16 --quiet
done
```
## Now lets assemble with megahit
```
mkdir 03_ASSEMBLIES

SAMPLES=`cut -f 1 samples3.txt | sed '1d'`
for SAMPLE in $SAMPLES; do
  megahit -1 00_Reads/${SAMPLE}/*_1.fq.gz \
          -2 00_Reads/${SAMPLE}/*_2.fq.gz \
          --out-dir 03_ASSEMBLIES/${SAMPLE} \
          --min-contig-len 1000 \
          -m 0.99
done
```
## Assess Assemblies with Quast
```
mkdir 04_ASSEMBLIES_QC

SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in $SAMPLES; do
  metaquast 03_ASSEMBLIES/${SAMPLE}/final.contigs.fa \
        --output-dir 04_ASSEMBLIES_QC/${SAMPLE} \
        --max-ref-number 0 \
        --threads 16
done
```
## BASALT
This bins genomes and then refines them
```
conda activate BASALT
mkdir 04_BASALT
```
You must move the files into a directory. The commands are NOT path friendly.
```
SAMPLES=`cut -f 1 samples_6.txt | sed '1d'`
for SAMPLE in $SAMPLES; do
mkdir 04_BASALT/${SAMPLE}
cp 03_ASSEMBLIES/${SAMPLE}/final.contigs.fa 04_BASALT/${SAMPLE}/${SAMPLE}_final.contigs.fa
cp 00_Reads/${SAMPLE}/*_1.fq.gz 04_BASALT/${SAMPLE}/${SAMPLE}_1.fq.gz
cp 00_Reads/${SAMPLE}/*_2.fq.gz 04_BASALT/${SAMPLE}/${SAMPLE}_2.fq.gz
done

SAMPLES=`cut -f 1 /blue/hlaughinghouse/flefler/BLCC_Genomes/samples_2_1.txt | sed '1d'`
for SAMPLE in $SAMPLES; do
	
	# Set a useful name for log files
	N="${SAMPLE}";

	# Set up command
	CMD="cd ${SAMPLE} && BASALT -a ${SAMPLE}_final.contigs.fa -s ${SAMPLE}_1.fq.gz,${SAMPLE}_2.fq.gz -t 28 -m 218 --min-cpn 80 --max-ctn 20 && cd .."
	
	# If you have SLURM engine, do something like this:
	sbatch -A hlaughinghouse -J ${N} -c 28 --mem=218G -o ${N}.o -e ${N}.e --export=ALL --mail-type=ALL --mail-user=flefler@ufl.edu -t 196:00:00  --wrap="${CMD}"
	
done
```
## BASALT for F46 with long read
```
cd ${SAMPLE} && BASALT -a ${SAMPLE}_contigs.fa -s ${SAMPLE}_1.fq,${SAMPLE}_2.fq -l ${SAMPLE}_lr.fq -t 8 -m 100 --min-cpn 80 --max-ctn 20 && cd ..
```
# Genome quality and other stastics
## Packages used 
[gtdbtk](https://github.com/Ecogenomics/GTDBTk)
[CheckM2](https://github.com/chklovski/CheckM2)
[CoverM](https://github.com/wwood/CoverM)
[SeqTK](https://github.com/lh3/seqtk)

## Rename Bins and move them
```
main_directory="/blue/hlaughinghouse/flefler/BLCC_Genomes/04_BASALT"

for directory in "$main_directory"/*/Final_bestbinset; do
    prefix="$(basename "$(dirname "$directory")")_"
    cd "$directory" || continue
    
    for file in *; do
        if [[ "$file" == "$prefix"* ]]; then
            continue  # Skip files that already have the prefix
        fi
        new_name="${prefix}${file}"
        mv -- "$file" "$new_name"
        echo "Renamed $file to $new_name"
    done
done

mkdir 05_GENOMES
SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in $SAMPLES; do
    mv 04_BASALT/${SAMPLE}/Final_bestbinset/*.fa 05_GENOMES
done
```
## gzip files
```
find /blue/hlaughinghouse/flefler/BLCC_Genomes/05_GENOMES -type f -print0 | xargs -0 gzip
```
## CheckM2
Used to assess completion and contamination
```
conda activate checkm2
mkdir 06_checkM

#Set job name
N="checkm2"
checkm2 predict --threads 28 --input 05_GENOMES --output-directory 06_checkM -x .fa.gz --quiet --remove_intermediates
```
## CoverM
Used to determine mean coverage
```
conda activate reassemble
mkdir 07_BAMFILES

for SAMPLE in $SAMPLES; do
    N=${SAMPLE}_coverm
    R1="/blue/hlaughinghouse/flefler/BLCC_Genomes/00_Reads/${SAMPLE}/*_1.fq.gz"
    R2="/blue/hlaughinghouse/flefler/BLCC_Genomes/00_Reads/${SAMPLE}/*_2.fq.gz"
    CMD="coverm genome -1 ${R1} -2 ${R2} --genome-fasta-files 05_GENOMES/${SAMPLE}_* --output-file 09_COVERM/${SAMPLE}_output.tsv \
    -x .gz --threads 8 --methods mean relative_abundance --bam-file-cache-directory process/bamcache"
    sbatch -A hlaughinghouse -J ${SAMPLE} -c 8 --mem=100G -o 99_logs/${SAMPLE}_h_o -e 99_logs/${SAMPLE}_h_e --export=ALL --mail-type=ALL --mail-user=flefler@ufl.edu -t 196:00:00  --wrap="${CMD}"
done
```
## GTDB
Used to assign taxonomy
```
conda activate gtdbtk-2.4.0
mkdir 08_gtdbtk
N="gtdb"
gtdbtk classify_wf --genome_dir 05_GENOMES --mash_db /blue/hlaughinghouse/flefler --out_dir 08_gtdbtk -x .fa.gz --cpus 28 --pplacer_cpus 28 --force
```

## SeqTK
Used to gather other stats, e.g., N50, length
```
conda activate reassemble

#general QC file
seqkit stats -Ta *.fa.gz | csvtk tab2csv -o SeqTK_Output.csv

#CheckM file
csvtk tab2csv quality_report.tsv | csvtk rename -f 1 -n file -o quality_report.csv

#GTDB file
csvtk tab2csv gtdb/gtdbtk.bac120.summary.tsv | csvtk rename -f 1 -n file -o gtdb/gtdbtk.bac120.summary.csv

#coverM file
SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in $SAMPLES; do
    cat ${SAMPLE}_output.tsv | csvtk rename -t -f 1,2,3 -n file,meancoverage,relativeabundance -o ${SAMPLE}_output2.csv
done

csvtk concat -t *2.csv | csvtk tab2csv | csvtk filter2 -f '$file!="unmapped"' -o coverMoutput.csv

csvtk join -f 1 SeqTK_Output.csv coverm/coverMoutput.csv checkm/quality_report.csv gtdb/gtdbtk.bac120.summary.csv -o merged.csv

csvtk concat -k merged.csv process/selectedmags.csv -o genomeinfo_230524.csv
```
# Analyses
## Packages used 
[skani](https://github.com/bluenote-1577/skani)
[ezaai](https://github.com/endixk/ezaai)
[gtdbtk](https://github.com/Ecogenomics/GTDBTk)
[ModelTest-NG](https://github.com/ddarriba/modeltest)
[RAxML-NG](https://github.com/amkozlov/raxml-ng)
[GToTree](https://github.com/AstrobioMike/GToTree)
[IQ-TREE](https://github.com/iqtree/iqtree2)

## ANI
We used skani to determine ANI and AF between genomes
```
skani dist -t 3 -q Genomes/*.fa -r Genomes/*.fa -s 70 --medium -o new_ANI/all-to-all_results.txt
```

## Create a large phylogenomic tree
We used GTDB to to create a concatenated alignment of the filamentous cyanobacterial families and orders in GTDB with *Gloeobacter* as the outgroup
```
gtdbtk de_novo_wf --genome_dir /media/HDD_4T/Forrest/floridaenema/stuff/Genomes/BLCC_genomes -x fa --bacteria --outgroup_taxon g__Gloeobacter \
--taxa_filter f__Coleofasciculaceae,f__Desertifilaceae,f__Microcoleaceae,f__Oscillatoriaceae,f__Phormidiaceae_A,f__PCC-6304,f__Geitlerinemaceae,f__Geitlerinemaceae_A,f__Spirulinaceae,g__Oculatella,g__Elainella,f__FACHB-T130,o__Phormidesmiales,g__Gloeobacter \
--write_single_copy_genes --cpus 4 --out_dir /media/HDD_4T/Forrest/floridaenema/stuff/new_GTDB
```
Using modeltest-ng to detemine the evolutionary model
```
modeltest-ng -i Floridanema_.fasta -d aa -p 10 -t ml
```
Run the phylogenomic tree with raxml-NG
```
raxml-ng --all --msa Floridanema_.fasta --threads auto{10} --workers auto --model LG+I+G4+F --bs-trees autoMRE{1000} --tree Floridanema_.fasta.tree --prefix Floridanema_tree
```

## Creating a more refined phylogenomic tree of the Aerosakkonematales 
Using GToTree and iqtree
fasta_files.txt is a file which contains the paths to the genomes of interest
```
GToTree -f fasta_files.txt -H Cyanobacteria -N -n 16 -j 4 -o /media/HDD_4T/Forrest/floridaenema/stuff/new_GToTree/GToTree_Floridaenema -F
```
feed the GToTree output to iqtree, runs a model on each partition
```
iqtree -s /media/HDD_4T/Forrest/floridaenema/stuff/new_GToTree/GToTree_Floridaenema/Aligned_SCGs.faa \
       -p /media/HDD_4T/Forrest/floridaenema/stuff/new_GToTree/GToTree_Floridaenema/run_files/Partitions.txt \
       -m MFP -B 1000 -pre iqtree_out
```

# Making ANI and AAI figures in ```R```

## Install and load packages 
```
library(devtools)
install_github("jokergoo/ComplexHeatmap")
install_github("GuangchuangYu/ggtree")
install.packages("tidyverse")

library(tidyverse)
library(ggtree)
library(ComplexHeatmap)
```
## Upload and edit the phylogenomic tree of the Aerosakkonematales 
There is probably a more effective way to do this, but if it aint broke dont fix it
```
tree = ggtree::read.tree("Floridaenema_GToTree.nwk")

dendrogram <- ape::chronos(tree2)
row_cluster <- as.hclust(dendrogram)
col_cluster <- as.hclust(dendrogram)

ordered_data <- result_df[tree2$tip.label, tree2$tip.label]

row.names(ordered_data)

row.names(ordered_data)[row.names(ordered_data) == "BLCC_F154"] <- "F. flavialum BLCC-F154"
row.names(ordered_data)[row.names(ordered_data) == "BLCC_F50"] <- "F. flaviceps BLCC-F50"
row.names(ordered_data)[row.names(ordered_data) == "BLCC_F46"] <- "F. aerugineus BLCC-F46"
row.names(ordered_data)[row.names(ordered_data) == "BLCC_F167"] <- "F. evergladium BLCC-F167"
row.names(ordered_data)[row.names(ordered_data) == "BLCC_F43"] <- "Microseira sp. BLCC-F43"
row.names(ordered_data)[row.names(ordered_data) == "BLCC_F2"] <- "Ae. funiforme BLCC-F2"
row.names(ordered_data)[row.names(ordered_data) == "BLCC_F183"] <- "Ae. funiforme BLCC-F183"
row.names(ordered_data)[row.names(ordered_data) == "GCA_014696265.1"] <- "Aerosakkonema sp. FACHB-1375"
row.names(ordered_data)[row.names(ordered_data) == "GCA_949128025.1_PMM_0008_genomic"] <- "Argonema sp. MAG PMM_0008"
row.names(ordered_data)[row.names(ordered_data) == "GCA_949127755.1_PMM_0001_genomic"] <- "Argonema sp. MAG PMM_0001"
row.names(ordered_data)[row.names(ordered_data) == "GCA_023333585.1"] <- "Ar. antarcticum A004/B2"
row.names(ordered_data)[row.names(ordered_data) == "GCA_023333595.1"] <- "Ar. galeatum A003/A1"
row.names(ordered_data)[row.names(ordered_data) == "GCA_003486305.1"] <- "Microseira sp. UBA11371" #Cyanobacteria bacterium UBA11371 
row.names(ordered_data)[row.names(ordered_data) == "GCA_003486675.1"] <- "Microseira sp. UBA11372" #    Cyanobacteria bacterium UBA11372 
row.names(ordered_data)[row.names(ordered_data) == "GCA_001904725.1"] <- "Floridaenema sp. IAM M-71" #Phormidium ambiguum IAM M-71
row.names(ordered_data)[row.names(ordered_data) == "GCA_015207735.1"] <- "Floridaenema sp. LEGE 05292" #Phormidium sp. LEGE 05292
row.names(ordered_data)[row.names(ordered_data) == "GCA_020521235.1"] <- "M. wollei NIES-4236 "

colnames(ordered_data)[colnames(ordered_data) == "BLCC_F154"] <- "F. flavialum BLCC-F154"
colnames(ordered_data)[colnames(ordered_data) == "BLCC_F50"] <- "F. flaviceps BLCC-F50"
colnames(ordered_data)[colnames(ordered_data) == "BLCC_F46"] <- "F. aerugineus BLCC-F46"
colnames(ordered_data)[colnames(ordered_data) == "BLCC_F167"] <- "F. evergladium BLCC-F167"
colnames(ordered_data)[colnames(ordered_data) == "BLCC_F43"] <- "Microseira sp. BLCC-F43"
colnames(ordered_data)[colnames(ordered_data) == "BLCC_F2"] <- "Ae. funiforme BLCC-F2"
colnames(ordered_data)[colnames(ordered_data) == "BLCC_F183"] <- "Ae. funiforme BLCC-F183"
colnames(ordered_data)[colnames(ordered_data) == "GCA_014696265.1"] <- "Aerosakkonema sp. FACHB-1375"
colnames(ordered_data)[colnames(ordered_data) == "GCA_949128025.1_PMM_0008_genomic"] <- "Argonema sp. MAG PMM_0008"
colnames(ordered_data)[colnames(ordered_data) == "GCA_949127755.1_PMM_0001_genomic"] <- "Argonema sp. MAG PMM_0001"
colnames(ordered_data)[colnames(ordered_data) == "GCA_023333585.1"] <- "Ar. antarcticum A004/B2"
colnames(ordered_data)[colnames(ordered_data) == "GCA_023333595.1"] <- "Ar. galeatum A003/A1"
colnames(ordered_data)[colnames(ordered_data) == "GCA_003486305.1"] <- "Microseira sp. UBA11371" #Cyanobacteria bacterium UBA11371 
colnames(ordered_data)[colnames(ordered_data) == "GCA_003486675.1"] <- "Microseira sp. UBA11372" #    Cyanobacteria bacterium UBA11372 
colnames(ordered_data)[colnames(ordered_data) == "GCA_001904725.1"] <- "Floridaenema sp. IAM M-71" #Phormidium ambiguum IAM M-71
colnames(ordered_data)[colnames(ordered_data) == "GCA_015207735.1"] <- "Floridaenema sp. LEGE 05292" #Phormidium sp. LEGE 05292
colnames(ordered_data)[colnames(ordered_data) == "GCA_020521235.1"] <- "M. wollei NIES-4236 "
```

## Load ANI and AAI results
```
ani <- read.delim("all-to-all_results.txt") %>% select(c(Ref_file, Query_file, ANI))

ani_result_df <- ani %>% pivot_wider(names_from = Ref_file, values_from = ANI, values_fn = mean) %>% column_to_rownames(var = "Query_file")

aai <- read.delim("AAI_output.tsv") %>% select(c(Label.1, Label.2, AAI))

aai_result_df <- aai %>% pivot_wider(names_from = Label.1, values_from = AAI, values_fn = mean) %>% column_to_rownames(var = "Label.2")
```

## Create and plot the heatmaps
```
p1 = ComplexHeatmap::pheatmap(as.matrix(ordered_data), legend_breaks = c(70,100), display_numbers = TRUE,
                              number_color = "black", cluster_rows = row_cluster, cluster_cols = col_cluster, fontsize_col = 8,
                              fontsize_number = 6, name = "ANI", angle_col = "45", number_format = "%.1f", column_title = "ANI")

p2 = ComplexHeatmap::pheatmap(as.matrix(ordered_data), legend_breaks = c(70,100), display_numbers = TRUE,
                              number_color = "black", cluster_rows = row_cluster, cluster_cols = col_cluster, fontsize_col = 8,
                              fontsize_number = 6, name = "AAI", angle_col = "45", number_format = "%.1f", column_title = "AAI")

p1 + p2
```



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

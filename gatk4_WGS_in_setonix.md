# User guide
**Guide on setting up and running GATK4 pipeline in Setonix**
## useful resources
#### https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/
#### https://github.com/gencorefacility/variant-calling-pipeline-gatk4
#### https://eriqande.github.io/eca-bioinf-handbook/
#### https://www.jianshu.com/p/aefd9a0eb67d/
#### https://jnmaloof.github.io/BIS180L_web/slides/
#### https://github.com/ShifengCHENG-Laboratory/WWWG2B/tree/main/

## install required tools from conda env file
```bash
/scratch/pawsey0399/cxiao/tools
git clone https://github.com/gencorefacility/variant-calling-pipeline-gatk4
cp variant-calling-pipeline-gatk4/bin/parse_metrics.sh your/active/environment/path/

conda env create -f gatk4_seto_nf_env.yaml
conda activate nf-env
```
## fastqc check
```bash
## fastqc.sh
find /scratch/pawsey0399/yjia/skylar/WGS/ -name "*fq.gz"|grep trimmed| while read R;
do
	fastqc $R -o raw_fastqc_output -t 128
done

## multiqc.conf
conda activate bio
srun --export=all -n 1 -c 64 multiqc -i skylar_raw_qc -n raw_qc raw_fastqc_output &> multiqc_log.txt
```
## trim fastq using fastp
```bash
## workdir
/scratch/pawsey0399/cxiao/WGS/01.RawData

## generate fastq pairs
find . -iname "*.fq.gz" | sort | paste - - > SPLIT

## trim fastq using fastp tool
while read R1 R2;
do
	sample="trimmed_"$(basename $R1|cut -d '_' -f1-4)
                fastp \
-i $R1 \
-I $R2 \
-o ./fastq_trimmed/$sample"_1.fq.gz" \
-O ./fastq_trimmed/$sample"_2.fa.gz" \
--qualified_quality_phred 20 \
--cut_front \
--cut_tail \
--length_required 30 \
--detect_adapter_for_pe \
--thread 128
done < SPLIT
```

## Solution 1 = Run GATK step-by-step
#### index reference genome
```bash
#!/bin/bash --login

#SBATCH --job-name=picard
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
#SBATCH --export=NONE

conda activate nf-env
REF=/scratch/pawsey0399/yjia/wheat/WGS/201216_Fielder_pseudomolecules_V1+unanchored_contigs.fasta.gz
PJAR=/scratch/pawsey0399/yjia/tools/miniconda3/envs/nf-env/share/picard-2.18.29-0/picard.jar

srun --export=all -n 1 -c 15 bwa index $REF
srun --export=all -n 1 -c 64 samtools faidx $REF
srun --export=all -n 1 -c 64 java -jar $PJAR CreateSequenceDictionary \
   R=201216_Fielder_pseudomolecules_V1+unanchored_contigs.fasta.gz \
   O=201216_Fielder_pseudomolecules_V1+unanchored_contigs.dict
```
#### https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently\
```bash
#!/bin/bash --login

#SBATCH --job-name=S100
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --account=ubuntu
#SBATCH --mem=50G
#SBATCH --export=NONE

source /data/tools/miniconda3/bin/activate nf-env
export REF=/data/skylar/Reference_genome_old/NLL_v2.fa
export PJAR=/data/tools/miniconda3/envs/nf-env/share/picard-2.18.29-0/picard.jar
srun --export=all -n 1 -c 15 bwa mem -M -t 15 -R '@RG\tID:NLL100\tLB:NLL100\tPL:ILLUMINA\tPM:HISEQ\tSM:NLL100' $REF NLL100_1.fq.gz NLL100_2.fq.gz |samtools view -@ 6 -Sb - | samtools sort -@ 6 -o NLL100.bam ## read alignment
srun --export=all -n 1 -c 15 gatk MarkDuplicatesSpark -I NLL100.bam -O mark_fix_NLL100.bam ## mark duplicates
srun --export=all -n 1 -c 15 rm NLL100.bam
srun --export=all -n 1 -c 15 gatk FixMateInformation -MC true -I mark_fix_NLL100.bam ## fixmate
srun --export=all -n 1 -c 15 samtools index -@ 15 mark_fix_NLL100.bam ## index bam, use -c for long chromosome
srun --export=all -n 1 -c 15 rm mark_fix_NLL100.bam.sbi ## this index file interferes with downstream
srun --export=all -n 1 -c 15 gatk HaplotypeCallerSpark -R $REF -I mark_fix_NLL100.bam -ERC GVCF -O mark_fix_NLL100.bam.g.vcf ## haplotype call
## use HaplotypeCaller, some bam csi index could not recognized by HaplotypeCallerSpark
srun --export=all -n 1 -c 15 rm mark_fix_NLL100.bam
srun --export=all -n 1 -c 15 bgzip -@ 15 mark_fix_NLL100.bam.g.vcf ## compress gvcf file
srun --export=all -n 1 -c 5 gatk IndexFeatureFile -I mark_fix_NLL100.bam.g.vcf.gz ## index gvcf file, produce tbi index file, or idx file for long chromosome which does not work downstream

## combine gvcfs for multiple samples
srun --export=all -n 1 -c 15 gatk --java-options "-Xmx60g -Xms60g" GenomicsDBImport \ ## use ***GenomicsDBImport*** generate genomicsdb for each chromosome downstream genotype call
	--genomicsdb-workspace-path CHROM_NAME \
       -V mark_fix_NLL001.bam.g.vcf.gz \
       -V mark_fix_NLL002.bam.g.vcf.gz \
       -V mark_fix_NLL003.bam.g.vcf.gz \
	...
	-L CHROM_NAME
Or srun --export=all -n 1 -c 15 gatk CombineGVCFs \ ## use ***CombineGVCFs***
	   -R ../Reference_genome/NLL_v2.fa \
       -V mark_fix_NLL001.bam.g.vcf.gz \
       -V mark_fix_NLL002.bam.g.vcf.gz \
       -V mark_fix_NLL003.bam.g.vcf.gz \
	...
	-O merged.g.vcf.gz

## generate genotype vcf
srun --export=all -n 1 -c 15 gatk --java-options "-Xmx40g" GenotypeGVCFs \
	-R /data/skylar/Reference_genome/NLL_v2.fa \
	-V NLL_old.g.vcf.gz \
	-O NLL_old.raw.vcf.gz
Or srun --export=all -n 1 -c 15 gatk --java-options "-Xmx64g" GenotypeGVCFs \
	-R /scratch/pawsey0399/yjia/skylar/Reference_genome/NLL_v2.fa \
	-V gendb://HiC_scaffold_10 \
	-O HiC_scaffold_10.raw.vcf.gz

## concat multiple chromosomes vcf
srun --export=all -n 1 -c 64 bcftools concat -f vcf_file_list \
	-o 300NLL_merged.vcf.gz \
	--threads 64

## filter vcf
srun --export=all -n 1 -c 64 gatk VariantFiltration \
	-R /scratch/pawsey0399/yjia/skylar/Reference_genome/NLL_v2.fa \
	-V HiC_scaffold_20.raw.vcf.gz \
	-O filtered_HiC_scaffold_20.raw.vcf.gz \
	--filter-name "QD_filter" -filter "QD < 2.0" \
	--filter-name "FS_filter" -filter "FS > 60.0" \
	--filter-name "MQ_filter" -filter "MQ < 40.0"
```
## merge gvcf files without indexing
```bash
## https://github.com/dnanexus-rnd/GLnexus/wiki/Getting-Started
### static executable: for modern Linux x86-64 hosts, download glnexus_cli from the Releases page and chmod +x glnexus_cli

## add sample header to bam if missed
https://github.com/IARCbioinfo/addreplacerg-nf
samtools addreplacerg -r "@RG\tID:file_name\tPG:samtools addreplacerg\tSM:file_name}"
#!/bin/bash --login
 
#SBATCH --job-name=gln
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=10:00:00
#SBATCH --account=pawsey0399
#SBATCH --mem=980G
#SBATCH --export=NONE
conda activate nf-env
export PATH=$PATH:/scratch/pawsey0399/yjia/tools/
srun --export=all -n 1 -c 128 glnexus_cli --config gatk --bed TraesFLD2D01G513900.bed -m 980 --threads 128 2Dvcf/*.g.vcf.gz > TraesFLD2D01G513900.bcf
srun --export=all -n 1 -c 128 bcftools view --threads 64 TraesFLD2D01G513900.bcf | bgzip -@ 4 -c > TraesFLD2D01G513900.vcf.gz

## snp annotation
### build snpeff database
cp cds.fa genes.gff protein.fa sequences.fa sequences.fa.fai /data/tools/snpEff/data/fielder
echo "fielder.genome : fielder" >> snpEff.config
java -jar snpEff.jar build -gff3 -v fielder
java -Xmx8g -jar snpEff.jar fielder /data/wheat/fielder/TraesFLD2D01G513900.vcf.gz >  /data/wheat/fielder/TraesFLD2D01G513900.annotated.vcf
bgzip TraesFLD2D01G513900.annotated.vcf
tabix -C TraesFLD2D01G513900.annotated.vcf.gz
bcftools annotate -x ^FORMAT/GT TraesFLD2D01G513900.annotated.vcf.gz > TraesFLD2D01G513900.annotated.simple.vcf
```

## Solution 2 = Run GATK nextflow
```bash
## split jobs for setonix
find fastq_trimmed -name "*.fq.gz" | sort | paste - - > new_pairs
split -a 3 -d -l 1 new_pairs input_

## mv.sh
ls input*|while read R;
do
	mkdir $R"_files"
	cat $R| while read P1 P2; do mv $P1 $P2 $R"_files";done
done

## prepare nextflow
for i in $(ls --color=never -d input*/);do cp nextflow.config $i;done
for i in $(ls --color=never -d input*/);do VAR=$(basename $i);sed -i "s/IN_DIR/$VAR/" $i"nextflow.config";done

for i in $(ls --color=never -d input*/);do cp nextflow.conf $i;done
for i in $(ls --color=never -d input*/);do VAR=$(basename $i);sed -i "s/INPUT/$VAR/" $i"nextflow.conf";done

for i in $(ls --color=never -d input*/);do cp main.nf $i;done

## submit jobs
for i in $(ls --color=never -d input*/);do cd $i;sbatch nextflow.conf;cd -;done

## nextflow.conf
#!/bin/bash --login

#SBATCH --job-name=input_004_files
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=60:00:00
#SBATCH --account=pawsey0399
###SBATCH --mem=64G
#SBATCH --exclusive
#SBATCH --export=NONE

##activate env
conda activate nf-env

export SNPEFF_JAR=/scratch/pawsey0399/yjia/tools/snpEff/snpEff.jar
export PICARD_JAR=/scratch/pawsey0399/yjia/tools/miniconda3/envs/nf-env/share/picard-2.18.29-0/picard.jar

srun --export=all -n 1 -c 128  nextflow run main.nf -c nextflow.config
```
## process results
```bash
## find succeeded runs
for file in ./input_*_files/slurm*;do if grep -q "Succeeded" $file; then echo $file" succeeded";fi;done > TMP_succeeded

## find resumed runs
find input_*_files/ -name "nextflow.conf" -type f |while read R;do if grep -q "resume" $R;then echo $R;fi;done > TMP_resumed_run

cat TMP_resumed_run|cut -d '/' -f1|while read R;do grep $R TMP_succeeded ;done
grep -f running_jobs -v TMP_resumed_run > resume_submit_again
cat resume_submit_again |cut -d '/' -f1|while read R;do cd $R;sbatch nextflow.conf;cd -;done

## check which jobs with errors
for file in ./input_*_files/slurm*;do if grep -q "error" $file; then echo $file" contains error";fi;done > TMP

## save deduplicated bam for later use
find . -type d -name "dedup_sorted"|while read R;do mv $R/* dedup_bam_sorted/;done

## save recallibriated bam for later use
find . -type d -name "bqsr"|while read R;do mv $R/* recal_bams/;done

## copy annotated snp vcf files
find . input_*_output/out/snpeff/ -type f -name "*ann.vcf" | grep -v nextflow_work_dir | while read R;do cp $R /scratch/pawsey0399/yjia/skylar/GATK_output2;done

## compress, index, and merge vcf
htsfile compressed.vcf.gz ## check file type
bcftools view file.vcf -Oz --threads 30 -o file.vcf.gz or bgzip file.vcf ## compress vcf
bcftools index file.vcf.gz --threads 30 or tabix file.vcf.gz ## index vcf
bcftools merge file*.vcf.gz -Oz --threads 30 -o merged.files.vcf.gz
bcftools concat file*.vcf.gz -Oz --threads 64 -o merged.snps.vcf.gz

##variant filteration
bcftools view -i 'INFO/QD > 2 & MQ > 40 & INFO/SOR < 4 & INFO/FS < 60 & INFO/MQRankSum > -12.5 & INFO/ReadPosRankSum > -8'
bcftools view -i 'QUAL > 200 & MQ > 50 & INFO/DP > 5 & INFO/QD > 2 & INFO/SOR < 3 & INFO/FS < 60 & INFO/MQRankSum > -12.5 & INFO/ReadPosRankSum > -8'
https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants \

```

## SNP marker filteration using plink
#### https://zzz.bwh.harvard.edu/plink/tutorial.shtml
#### https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/
#### https://speciationgenomics.github.io/pca/

```bash
## Convert SNP vcf file into plink binary file
plink --vcf input.vcf --make-bed --allow-extra-chr --double-id --out output_prefix

## filter snp by MAF and missingness
plink --bfile output_prefix --maf 0.05 --geno 0.2 --mind 0.2 --make-bed --out filtered_output_prefix --allow-extra-chr ## geno missing, mind individual missing, maf minor allele missing

## filter snp by chromosomes CHROM=HiC_scaffold_1 HiC_scaffold_10 HiC_scaffold_11 .......
plink --bfile filtered_merged_snps --chr $CHROM --out chromosome_only_genotype --make-bed --allow-extra-chr

## filter genotype data by samples
plink --bfile chromosome_only_genotype --keep sample_list.txt --make-bed --out test --allow-extra-chr ## sample_list.txt needs to have two columns

## population stratification analyses
plink --bfile filtered_output_prefix --genome --out filtered_output_prefix --allow-no-sex --allow-extra-chr
plink --bfile filtered_output_prefix --read-genome filtered_output_prefix.genome --cluster --ppc 0.05 --allow-no-sex --allow-extra-chr

## association anlayses with phenotype and population clustering; --mh requires control/case information
plink --bfile filtered_output_prefix --pheno your.phe --within plink.cluster --assoc --adjust --allow-no-sex --allow-extra-chr --out GWAS_output

## filter genotype for target chromosome id only
plink --bfile filtered_updated_large_snp --chr $(cat chromosome_id.txt) --make-bed --out chrosome_only_genotype --allow-no-sex --allow-extra-chr

## prune snp based on LD calculation
plink2 --bfile chromosome_only_genotype --set-all-var-ids @:# --make-bed --out test --allow-extra-chr ## name snp if not
plink --bfile test --indep-pairwise 50 5 0.5 --out LD_pruned --allow-extra-chr ## prune snp based LD, per 50 snp
plink --bfile test --extract LD_pruned.prune.in --out LD_pruned --make-bed --allow-extra-chr ## filter SNP genotype data

## modify/offset snp positions
### fix_chrpos.sh
SNP="QUAL200MQ50DP6_more_HaplotypeCaller.chr1A_part2.vcf.gz.snp.vcf.gz"
zcat $SNP | awk -F '\t' -v OFS='\t' 'BEGIN {
first_file_data["chr1A_part2"] = 424817035
first_file_data["chr2C_part2"] = 426814308
first_file_data["chr3C_part2"] = 476205391
first_file_data["chr4C_part2"] = 360718161
first_file_data["chr5C_part2"] = 293874722
first_file_data["chr6C_part2"] = 297374895
first_file_data["chr7C_part2"] = 442998614
first_file_data["chr2D_part2"] = 199338498
first_file_data["chr7D_part2"] = 469736193
} {if ($1 in first_file_data) $2 += first_file_data[$1]; print}' | bgzip > modified_snp.vcf.gz

## recode and modify chromosome ID
plink --bfile LD_pruned --recode --out changeChr_LD_pruned --allow-extra-chr ## generate map file
sed -i 's/_part1//g;s/_part2//g' changeChr_LD_pruned.map ## modify chromosome IDs
plink --file changeChr_LD_pruned --make-bed --out changeChr_LD_pruned --allow-extra-chr ## convert back to binary files

## notes:
1. no LD pruning for gwas
2. plink --keep-allele-order ## plink auto assign allele: 0- 1-
3. plink --biallelic-only strict ## remove multi alleles
4. plink --set-missing-var-ids @:#
```
## bcftools cheatsheet
#### https://gist.github.com/elowy01/93922762e131d7abd3c7e8e166a74a0b

## create phylogeny tree from vcf
```bash
### solution 1. iqtree
## use plink to filter your genotype file, see above
## convert plink binary to vcf
plink2 --bfile LD_pruned --recode vcf id-paste=iid --out pruned --allow-extra-chr
## convert vcf to fasta or phylip using vcf2phylip at https://github.com/edgardomortiz/vcf2phylip
git clone https://github.com/edgardomortiz/vcf2phylip
python vcf2phylip/vcf2phylip.py -i LD_pruned.vcf --fasta --min-samples-locus 60
iqtree -s SNP_data.fasta -m GTR+ASC

### solution 2. vcf2pop webpage tool
https://github.com/sansubs/vcf2pop

### snphylo
https://github.com/thlee/SNPhylo/blob/master/docs/install_on_linux.rst
```

## GWAS using rMVP
### https://github.com/xiaolei-lab/rMVP
### https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2800123/
```bash
## install rMVP
conda create -n rMVP r=4.2
conda activate rMVP
conda install -c conda-forge r-rmvp
```
## using rMVP in jupyter lab in my nimbus server
```bash
# install rMVP in conda, from scratch
conda create --name rMVP
conda activate rMVP
conda install -c conda-forge r-base ## current default version 4.3.1
R
# within R, install IRkernel
install.packages("IRkernel")
IRkernel::installspec()
# within R, install rMVP
install.packages("rMVP")
q()
# launch jupyter lab
jupyter-lab --ip 0.0.0.0 --port 8888 --no-browser ## note the token key on your screen

# open your browser, paste the address below, and input the token:
http://146.118.64.19:8888/lab/

```
## running rMVP
```R
## https://github.com/yongjiam/chickpea_acid_soil/blob/main/rMVP_gwas.ipynb
library(rMVP)

## transform genotype vcf file, Full-featured function (Recommended)
MVP.Data(fileVCF="input.genotype.vcf",
         fileKin=FALSE,
         filePC=FALSE,
         out="mvp.genotype"
         )
## or plink binary genotype data
# Full-featured function (Recommended)
MVP.Data(fileBed="LD_pruned",
         filePhe=NULL,
         fileKin=FALSE,
         filePC=FALSE,       
         #priority="speed",
         #maxLine=10000,
         out="mvp.plink"
         )

## read genotye and phenotype data
genotype <- attach.big.matrix("mvp.genotype.desc")
phenotype <- read.table("input.phenotype",head=TRUE)
map <- read.table("mvp.genotype.map" , head = TRUE)

## calculate kinship and principle components from genotype data
MVP.Data.Kin(TRUE, mvp_prefix='mvp.plink', out='mvp')
Kinship <- attach.big.matrix("mvp.kin.desc")

MVP.Data.PC(TRUE, mvp_prefix='mvp.plink', out='mvp', pcs.keep=10)
Covariates_PC <- bigmemory::as.matrix(attach.big.matrix("mvp.pc.desc"))

## run GWAS
for(i in 2:ncol(phenotype)){
  imMVP <- MVP(
    phe=phenotype[, c(1, i)],
    geno=genotype,
    map=map,
    K=Kinship,
    CV.GLM=Covariates_PC,
    CV.MLM=Covariates_PC,
    CV.FarmCPU=Covariates_PC,
    nPC.GLM=5,
    nPC.MLM=3,
    nPC.FarmCPU=3,
    priority="speed",
    #ncpus=10,
    vc.method="BRENT",
    maxLoop=10,
    method.bin="static",
    #permutation.threshold=TRUE,
    #permutation.rep=100,
    threshold=0.05,
    method=c("GLM", "MLM", "FarmCPU"),
    #file.output=c("pmap", "pmap.signal", "plot", "log")
    file.output=TRUE
  )
  gc()
}
```

## understanding QQ plots
https://jnmaloof.github.io/BIS180L_web/slides/11_QQPlots.html#1



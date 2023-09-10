# User guide
**Guide on setting up and running GATK4 pipeline in Setonix using lupin WGS data as an example**
## useful resources
#### https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/
#### https://github.com/gencorefacility/variant-calling-pipeline-gatk4
#### https://eriqande.github.io/eca-bioinf-handbook/
#### https://www.jianshu.com/p/aefd9a0eb67d/

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
## prepare input files
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
## Run GATK nextflow
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
bcftools view file.vcf -Oz -o file.vcf.gz
bcftools index file.vcf.gz
bcftools merge file*.vcf.gz -Oz > merged.files.vcf.gz
```

## SNP marker filteration using plink
#### https://zzz.bwh.harvard.edu/plink/tutorial.shtml
#### https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/

```bash
## Convert SNP vcf file into plink binary file
plink --vcf input.vcf --make-bed --allow-extra-chr --double-id --out output_prefix

## filter snp by MAF and missingness
plink --bfile output_prefix --maf 0.05 --geno 0.2 --make-bed --out filtered_output_prefix

## population stratification analyses
plink --bfile filtered_output_prefix --genome --out filtered_output_prefix --allow-no-sex --allow-extra-chr
plink --bfile filtered_output_prefix --read-genome filtered_output_prefix.genome --cluster --ppc 0.05 --allow-no-sex --allow-extra-chr

## association anlayses with phenotype and population clustering; --mh requires control/case information
plink --bfile filtered_output_prefix --pheno your.phe --within plink.cluster --assoc --adjust --allow-no-sex --allow-extra-chr --out GWAS_output

## filter genotype for target chromosome id only
plink --bfile filtered_updated_large_snp --chr $(cat chromosome_id.txt) --make-bed --out chrosome_only_genotype --allow-no-sex --allow-extra-chr
```
## bcftools cheatsheet
#### https://gist.github.com/elowy01/93922762e131d7abd3c7e8e166a74a0b

## GWAS using rMVP
### https://github.com/xiaolei-lab/rMVP
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

## read genotye and phenotype data
genotype <- attach.big.matrix("mvp.genotype.desc")
phenotype <- read.table("input.phenotype",head=TRUE)
map <- read.table("mvp.genotype.map" , head = TRUE)

## run GWAS
for(i in 2:ncol(phenotype)){
  imMVP <- MVP(
    phe=phenotype[, c(1, i)],
    geno=genotype,
    map=map,
    #K=Kinship,
    #CV.GLM=Covariates,
    #CV.MLM=Covariates,
    #CV.FarmCPU=Covariates,
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


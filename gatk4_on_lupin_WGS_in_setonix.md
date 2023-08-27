# User guide
**Guide on setting up and running GATK4 pipeline in Setonix using lupin WGS data as an example**

## install required tools from conda env file
```bash
/scratch/pawsey0399/cxiao/tools
conda env create -f gatk4_seto_nf_env.yaml
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
## check which jobs with errors
for file in ./input_*_files/slurm*;do if grep -q "error" $file; then echo $file" contains error";fi;done > TMP

## copy annotated snp vcf files
find . input_*_output/out/snpeff/ -type f -name "*ann.vcf" | grep -v nextflow_work_dir | while read R;do cp $R /scratch/pawsey0399/yjia/skylar/GATK_output2;done

## compress, index, and merge vcf
bcftools view file.vcf -Oz -o file.vcf.gz
bcftools index file.vcf.gz
bcftools merge file*.vcf.gz -Ov | gzip > merged.files.vcf.gz
```

## SNP marker filteration using plink
```bash
## Convert SNP vcf file into plink binary file
plink --vcf input.vcf --make-bed --allow-extra-chr --double-id --out output_prefix

## filter snp by MAF and missingness
plink --bfile output_prefix --maf 0.05 --geno 0.2 --make-bed --out filtered_output_prefix
```


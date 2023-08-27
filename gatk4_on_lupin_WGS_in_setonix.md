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
## Run GATK nextflow pipeline
```bash
## split jobs for setonix

```

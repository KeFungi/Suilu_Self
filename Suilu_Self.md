# genome annotation
## funannotate
```bash
mkdir -p funannotate_suilu funannotate_aphal 
module purge
conda activate /nfs/turbo/lsa-tyjames/funannotate/funannotate
source /nfs/turbo/lsa-tyjames/funannotate/envvars.txt
module load openjdk

# normalize genome
funannotate clean -i ref/Suilu4.fasta -o funannotate_suilu/assembly.clean.fasta # no change
funannotate clean -i Aphal/ge/Aphal.fasta -o funannotate_aphal/assembly.clean.fasta # no change

funannotate sort -i funannotate_suilu/assembly.clean.fasta -o funannotate_suilu/assembly.clean.sort.fasta
funannotate sort -i funannotate_aphal/assembly.clean.fasta -o funannotate_aphal/assembly.clean.sort.fasta

funannotate mask -i funannotate_suilu/assembly.clean.sort.fasta -o funannotate_suilu/assembly.clean.sort.masked.fasta
funannotate mask -i funannotate_aphal/assembly.clean.sort.fasta -o funannotate_aphal/assembly.clean.sort.masked.fasta

# gene prediction
sbatch -J funannotate_predict -A tyjames1 -t 48:00:00 -c 36 --mem=24g  --wrap="funannotate predict -i funannotate_suilu/assembly.clean.sort.masked.fasta -o funannotate_suilu/ -s 'Suillus luteus' --strain "UH-Slu-Lm8-n1" --cpus 36 --optimize_augustus"
sbatch -J funannotate_predict -A tyjames1 -t 48:00:00 -c 36 --mem=24g  --wrap="funannotate predict -i funannotate_aphal/assembly.clean.sort.masked.fasta -o funannotate_aphal/ -s 'Amanita phalloides' --strain Alameda01 --cpus 36 --optimize_augustus"

# InterProScan
sbatch -J iprscan -A tyjames1 -t 48:00:00 -c 4 --mem=64g --wrap="funannotate iprscan -i funannotate_suilu/ --cpus 4 -m local --iprscan_path /nfs/turbo/lsa-tyjames/funannotate/interproscan-5.60-92.0/interproscan.sh"
sbatch -J iprscan -A tyjames1 -t 48:00:00 -c 4 --mem=64g --wrap="funannotate iprscan -i funannotate_aphal/ --cpus 4 -m local --iprscan_path /nfs/turbo/lsa-tyjames/funannotate/interproscan-5.60-92.0/interproscan.sh"

# annotation results
sbatch -J annotate -d afterok:63444576 -A tyjames1 -t 48:00:00 -c 20 --mem=64g --wrap="funannotate annotate -i funannotate_suilu/ --cpus 20 --iprscan funannotate_suilu/annotate_misc/iprscan.xml"
sbatch -J annotate -d afterok:63444577 -A tyjames1 -t 48:00:00 -c 20 --mem=64g --wrap="funannotate annotate -i funannotate_aphal/ --cpus 20 --iprscan funannotate_aphal/annotate_misc/iprscan.xml"

# extract results
mkdir funannotate
cp funannotate_suilu/assembly.clean.sort.masked.fasta funannotate/suilu.assembly.fasta
cp funannotate_suilu/annotate_results/Suillus_luteus_UH-Slu-Lm8-n1.gff3 funannotate_suilu/annotate_results/Suillus_luteus_UH-Slu-Lm8-n1.annotations.txt funannotate_suilu/annotate_results/Suillus_luteus_UH-Slu-Lm8-n1.cds-transcripts.fa funannotate_suilu/annotate_results/Suillus_luteus_UH-Slu-Lm8-n1.proteins.fa   funannotate

cp funannotate_aphal/assembly.clean.sort.masked.fasta funannotate/aphal.assembly.fasta
cp funannotate_aphal/annotate_results/Amanita_phalloides_Alameda01.gff3 funannotate_aphal/annotate_results/Amanita_phalloides_Alameda01.annotations.txt funannotate_aphal/annotate_results/Amanita_phalloides_Alameda01.cds-transcripts.fa funannotate_aphal/annotate_results/Amanita_phalloides_Alameda01.proteins.fa funannotate
```

# input SNP
```bash
vcftools --gzvcf bwa_haplotypecaller_finalvcf/s208_SNPfiltered.vcf.gz --keep metadata/suilu_s189.samples.txt --recode --recode-INFO-all --stdout | bgzip > bwa_haplotypecaller_finalvcf/s189_SNPfiltered.vcf.gz

# filter SNP 1.5 DP
vcftools --gzvcf bwa_haplotypecaller_finalvcf/s189_SNPfiltered.vcf.gz --site-mean-depth --out bwa_haplotypecaller_finalvcf/suilu
vcftools --gzvcf Aphal/Aphal_SNPfiltered.vcf.gz --site-mean-depth --out Aphal/aphal

MED=$(awk 'NR>1{print $3}' Aphal/aphal.ldepth.mean | sort -g | awk '{a[NR]=$1} END{print a[int(NR/2)]}')
MAXDP=$(echo "1.5*$MED" | bc -l)
vcftools --gzvcf Aphal/Aphal_SNPfiltered.vcf.gz --max-meanDP $MAXDP --max-missing 0.9 --min-alleles 2 --max-alleles 2 --recode -c | bgzip > bwa_haplotypecaller_finalvcf/Aphal_SNPDPfiltered.vcf.gz

MED=$(awk 'NR>1{print $3}' bwa_haplotypecaller_finalvcf/suilu.ldepth.mean | sort -g | awk '{a[NR]=$1} END{print a[int(NR/2)]}')
MAXDP=$(echo "1.5*$MED" | bc -l)
vcftools --gzvcf bwa_haplotypecaller_finalvcf/s189_SNPfiltered.vcf.gz  --max-meanDP $MAXDP --max-missing 0.9 --min-alleles 2 --max-alleles 2 --recode -c | bgzip > bwa_haplotypecaller_finalvcf/s189_SNPDPfiltered.vcf.gz
```


## plink files
```bash
# s208 (no clone, no mixed, no haploid; Slu1, Slu2): pca, admixture
mkdir plink
plink --allow-extra-chr --vcf bwa_haplotypecaller_finalvcf/s189_SNPDPfiltered.vcf.gz  --out plink/s189 --make-bed --set-missing-var-ids @:# --double-id --update-ids metadata/pop.fam.txt 
plink --allow-extra-chr --bfile plink/s189 --out plink/s189  --het

plink --allow-extra-chr --vcf bwa_haplotypecaller_finalvcf/s189_SNPDPfiltered.vcf.gz --out plink/s132 --make-bed --set-missing-var-ids @:# --double-id --update-ids metadata/s132.rename.fam.txt --keep metadata/s132.fam.txt
plink --allow-extra-chr --bfile plink/s132 --out plink/s132_admixture --max-alleles 2 --mac 2 --indep-pairwise 50 10 0.1 --geno 0.1
plink --allow-extra-chr --bfile plink/s132 --out plink/s132_admixture --extract plink/s132_admixture.prune.in --make-bed # biallelic
plink --allow-extra-chr --bfile plink/s132_admixture --out plink/s132_admixture --het
plink --allow-extra-chr --bfile plink/s132_admixture --out plink/s132_admixture --recode A

plink --allow-extra-chr --bfile plink/s132_admixture --out plink/s132o_admixture --keep metadata/s132o.fam.txt --update-ids metadata/s132o.rename.fam.txt --make-bed
plink --allow-extra-chr --bfile plink/s132o_admixture --out plink/s132o_admixture --recode A
plink --allow-extra-chr --bfile plink/s132o_admixture --out plink/s132o_admixture --het

plink --allow-extra-chr --vcf bwa_haplotypecaller_finalvcf/Aphal_SNPDPfiltered.vcf.gz --out plink/Aphalraw --make-bed --set-missing-var-ids @:# --double-id
plink --allow-extra-chr --bfile plink/Aphalraw --out plink/Aphal --mac 1 --make-bed 
```

# ROH
```bash
mkdir bcfroh_input bcfroh
plink --allow-extra-chr --bfile plink/s189 --out bcfroh_input/s189 --recode vcf bgz
plink --allow-extra-chr --bfile plink/s189 --out bcfroh_input/s189 --freq
tail -n +2 bcfroh_input/s189.afreq | sed -E 's/scaffold_[0-9]+://g' | awk 'BEGIN{OFS="\t"};{print $1,$2,$5}' | bgzip -c > bcfroh_input/AFs.tab.gz
tabix -s1 -b2 -e2 bcfroh_input/AFs.tab.gz
echo '##INFO=<ID=AF1KG,Number=1,Type=Float,Description="Allele Frequency">' > bcfroh_input/AFs.tab.gz.hdr

for file in bcfroh_input/s189.vcf.gz; do
  for sample in `bcftools query -l $file`; do
    bcftools view -Ou -s $sample $file | bcftools annotate -c CHROM,POS,AF1KG -h bcfroh_input/AFs.tab.gz.hdr -a bcfroh_input/AFs.tab.gz | bcftools roh --AF-tag AF1KG -M 0.00003 -o bcfroh/${sample}.roh.txt -G30
  done
done

for file in bcfroh/*.roh.txt
do grep '^RG' $file > ${file%.roh.txt}.roh_RG.txt; done
```

# MAT allele reconstruction
## bowtie approach
```bash
mkdir MATA_bowtie_raw
bowtie2-build MATA_ref/Suilu_MATA_refl8.fasta MATA_ref/Suilu_MATA_refl

seql=$(cat metadata/MATA_run.txt)
for seq in $seql
do 
  R1path=raw_fastq/${seq}_R1all_trim.fastq 
  R2path=raw_fastq/${seq}_R2all_trim.fastq 
  outpath=MATA_bowtie_raw/${seq}.fastq
  bowtie2 --local --very-sensitive-local -N 1 -k 1 -q --no-unal -x MATA_ref/Suilu_MATA_refl -1 $R1path -2 $R2path | samtools fastq - > $outpath
done

for seq in $seql
  do
    scripts/grep_allfastq.sh $seq
  done
```

# SnpEff
## Build SnpEff databases
```bash
mkdir snpeff snpeff/data/
mkdir snpeff/data/aphal snpeff/data/suilu

samtools faidx ./funannotate/suilu.assembly.fasta
samtools faidx ./funannotate/aphal.assembly.fasta

cp ./funannotate/Suillus_luteus_UH-Slu-Lm8-n1.gff3 ./snpeff/data/suilu/genes.gff
cp ./funannotate/Amanita_phalloides_Alameda01.gff3 ./snpeff/data/aphal/genes.gff

cp ./funannotate/suilu.assembly.fasta ./snpeff/data/suilu/sequences.fa
cp ./funannotate/aphal.assembly.fasta ./snpeff/data/aphal/sequences.fa

cp ./funannotate/suilu.assembly.fasta.fai ./snpeff/data/suilu/sequences.fa.fai
cp ./funannotate/aphal.assembly.fasta.fai ./snpeff/data/aphal/sequences.fa.fai

cp ./funannotate/Amanita_phalloides_Alameda01.cds-transcripts.fa ./snpeff/data/aphal/cds.fa
cp ./funannotate/Amanita_phalloides_Alameda01.proteins.fa ./snpeff/data/aphal/protein.fa

cp ./funannotate/Suillus_luteus_UH-Slu-Lm8-n1.cds-transcripts.fa ./snpeff/data/suilu/cds.fa
cp ./funannotate/Suillus_luteus_UH-Slu-Lm8-n1.proteins.fa ./snpeff/data/suilu/protein.fa

cat > snpeff/snpeff.config <<EOF
data.dir = ./snpeff/data/

aphal.genome : Amanita phalloides 10511
suilu.genome : Suillus luteus Suilu4
EOF

snpEff build -gff3 -v -c snpeff/snpeff.config  aphal > snpeff/snpeff_build_aphal.log 2>&1
snpEff build -gff3 -v -c snpeff/snpeff.config  suilu > snpeff/snpeff_build_suilu.log 2>&1
```

## Build VCF from PLINK
```bash
plink --allow-extra-chr --bfile plink/s189 --out snpeff/data/s189 --recode vcf bgz
tabix -f -p vcf snpeff/data/s189.vcf.gz

plink --allow-extra-chr --bfile plink/Aphal --out snpeff/data/aphal --recode vcf bgz
tabix -f -p vcf snpeff/data/aphal.vcf.gz
```

## annotate
```bash
  ann() {  # $1 in.vcf.gz  $2 contig map (vcf_name ... funannotate_name)  $3 fasta  $4 db  $5 out.vcf.gz
    bcftools view -t $(cut -f1 $2 | paste -sd,) $1 -Ou \
    | bcftools annotate --rename-chrs <(awk -v OFS='\t' '{print $1,$NF}' $2) -Ou \
    | bcftools norm --check-ref e -f $3 -Ov \
    | snpEff -c snpeff/snpeff.config -dataDir $PWD/snpeff/data $4 2> snpeff/snpeff_ann_$4.log \
    | bcftools annotate --rename-chrs <(awk -v OFS='\t' '{print $NF,$1}' $2) -Oz -o $5
    tabix -f -p vcf $5
  }

  ann snpeff/data/aphal.vcf.gz funannotate/aphal_contig_map.tsv ./snpeff/data/aphal/sequences.fa aphal snpeff/aphal.ann.vcf.gz
  ann snpeff/data/s189.vcf.gz  funannotate/suilu_contig_map.tsv ./snpeff/data/suilu/sequences.fa suilu snpeff/suilu.ann.vcf.gz

```

## HIGH and LOW counts
```bash
# HIGH / LOW site lists (first ANN entry)
for sp in suilu aphal; do
  bcftools query -f '%CHROM\t%POS\t%INFO/ANN\n' snpeff/$sp.ann.vcf.gz \
  | awk -F'\t' -v o=snpeff/$sp '{split($3,e,","); split(e[1],f,"|"); if(f[3]=="HIGH"||f[3]=="LOW") print $1"\t"$2 > (o"_"f[3]".sites")}' &
done; wait

# sample -> population (FID) maps: Suilu s132, Aphal 17
awk 'NR==FNR{p[$2]=$1; next} {i=$0; sub(/^[^_]*_/,"",i); if(i in p) print $0"\t"$0"\t"p[i]"\t"i}' metadata/s132.fam.txt <(bcftools query -l snpeff/suilu.ann.vcf.gz) > snpeff/suilu.update_ids.txt
bcftools query -l snpeff/aphal.ann.vcf.gz | awk '{print $1"\t"$1"\tAphal-Drake\t"$1}' > snpeff/aphal.update_ids.txt

# HIGH / LOW VCFs (biallelic SNPs) -> plink, variable in the sample set
for sp in suilu aphal; do for imp in HIGH LOW; do
  bcftools view -T snpeff/${sp}_${imp}.sites -S <(cut -f1 snpeff/$sp.update_ids.txt) -m2 -M2 -v snps snpeff/$sp.ann.vcf.gz -Oz -o snpeff/${sp}_${imp}.vcf.gz
  plink --allow-extra-chr --vcf snpeff/${sp}_${imp}.vcf.gz --double-id --set-missing-var-ids @:# --update-ids snpeff/$sp.update_ids.txt --mac 1 --make-bed --out plink/${sp}_${imp}
done; done
```

## Genome denominators
```bash
{ printf "species\tgenes\tprotein_coding_genes\tCDS_bp\n"
for sp in aphal:Amanita_phalloides_Alameda01 suilu:Suillus_luteus_UH-Slu-Lm8-n1; do
  awk -F'\t' -v sp=${sp%%:*} 'BEGIN{OFS="\t"} $3=="gene"{g++} $3=="mRNA"{m++} $3=="CDS"{c+=$5-$4+1} END{print sp, g, m, c}' funannotate/${sp#*:}.gff3
done; } > funannotate/genome_denominators.tsv
```

# Amanita phalloides data
```bash
source ~/miniconda3/etc/profile.d/conda.sh; conda activate genomics
module load Bioinformatics gatk/4.3.0.0
mkdir -p ge metadata sra_cache raw_fastq cutadapt_reports \
         bwa_mapping bwa_depth \
         bwa_haplotypecaller_gvcf bwa_haplotypecaller_log bwa_haplotypecaller_finalvcf \
         scripts
```

## 0. Index reference + soft-mask BED
```bash
ln -sf ../GCA_032173475.1_ASM3217347v1_genomic.fna ge/Aphal.fasta
bwa index ge/Aphal.fasta
samtools faidx ge/Aphal.fasta
picard CreateSequenceDictionary R=ge/Aphal.fasta O=ge/Aphal.dict

# soft-mask BED (lowercase regions in the assembly)
cat ge/Aphal.fasta | perl -lne 'if(/^(>.*)/){ $head=$1 } else { $fa{$head} .= $_ } END{ foreach $s (sort(keys(%fa))){ print "$s\n$fa{$s}\n" }}' \
  | perl -lne 'if(/^>(\S+)/){ $n=$1 } else { while(/([a-z]+)/g){ printf("%s\t%d\t%d\n",$n,pos($_)-length($1),pos($_)) } }' \
  > ge/Aphal.mask.bed
awk -F'\t' '{SUM+=$3-$2} END{print SUM}' ge/Aphal.mask.bed   # 5,464,260 bp masked
```

```bash
tail -n +2 Aphal_DrakeCA_HDnot1_dedup.sra.csv | awk -F, 'NF>=2 && $1!=""{print $2}' > metadata/accession_list.txt
tail -n +2 Aphal_DrakeCA_HDnot1_dedup.sra.csv | awk -F, 'NF>=2 && $1!=""{split($1,a,"_"); print a[1]"\t"$1"\t"$2}' > metadata/lane_map.tsv
cut -f1 metadata/lane_map.tsv | sort -u > metadata/sample_list.txt   # 17 samples
cut -f1 ge/Aphal.fasta.fai > metadata/contig_list.txt                # 605 contigs
```

## 1. SRA download
```bash
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" metadata/lane_map.tsv)
sample=$(echo "$line" | cut -f1); lane_label=$(echo "$line" | cut -f2); srr=$(echo "$line" | cut -f3)
[[ -s raw_fastq/${lane_label}_R1.fastq.gz && -s raw_fastq/${lane_label}_R2.fastq.gz ]] && exit 0

prefetch --max-size 100g -O sra_cache/ "$srr"
TMP=${WD}/sra_cache/.tmp_${srr}; mkdir -p "$TMP"
fasterq-dump --split-files -e 4 -t "$TMP" -O raw_fastq/ "sra_cache/${srr}/${srr}.sra"
rm -rf "$TMP"

mv raw_fastq/${srr}_1.fastq raw_fastq/${lane_label}_R1.fastq
mv raw_fastq/${srr}_2.fastq raw_fastq/${lane_label}_R2.fastq
gzip -f raw_fastq/${lane_label}_R1.fastq raw_fastq/${lane_label}_R2.fastq
rm -rf sra_cache/${srr}
```

## 2. Lane merge → per-sample fastq
```bash
merge_one() {
    local s="$1"
    [[ -s raw_fastq/${s}_R1.fastq.gz && -s raw_fastq/${s}_R2.fastq.gz ]] && return 0
    mapfile -t lanes < <(awk -F'\t' -v x="$s" '$1==x{print $2}' metadata/lane_map.tsv)
    cat $(printf "raw_fastq/%s_R1.fastq.gz " "${lanes[@]}") > raw_fastq/${s}_R1.fastq.gz.tmp
    cat $(printf "raw_fastq/%s_R2.fastq.gz " "${lanes[@]}") > raw_fastq/${s}_R2.fastq.gz.tmp
    mv raw_fastq/${s}_R1.fastq.gz.tmp raw_fastq/${s}_R1.fastq.gz
    mv raw_fastq/${s}_R2.fastq.gz.tmp raw_fastq/${s}_R2.fastq.gz
}
export -f merge_one
xargs -a metadata/sample_list.txt -P 4 -I {} bash -c 'set -euo pipefail; merge_one "$@"' _ {}
```

## 3. Conditional read-level subsample (gated on >30× expected coverage)
```bash
s=$(sed -n "${SLURM_ARRAY_TASK_ID}p" metadata/sample_list.txt)
R1in=raw_fastq/${s}_R1.fastq.gz; R2in=raw_fastq/${s}_R2.fastq.gz
R1out=raw_fastq/${s}_R1p.fastq.gz; R2out=raw_fastq/${s}_R2p.fastq.gz

[[ -e $R1out && -e $R2out && -s $R1out && -s $R2out ]] && exit 0
sz=$(stat -c %s $R1in); THRESH=$((1700*1024*1024))
if (( sz < THRESH )); then
    ln -sf ${s}_R1.fastq.gz $R1out; ln -sf ${s}_R2.fastq.gz $R2out
else
    seqtk sample -2 -s 100 $R1in 10000000 | pigz -p 4 -c > ${R1out}.tmp
    seqtk sample -2 -s 100 $R2in 10000000 | pigz -p 4 -c > ${R2out}.tmp
    mv ${R1out}.tmp $R1out; mv ${R2out}.tmp $R2out
fi
```

## 4. Trim adaptors (cutadapt q25, m80)
```bash
s=$(sed -n "${SLURM_ARRAY_TASK_ID}p" metadata/sample_list.txt)
cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -q 25 -m 80 -j 4 \
    -o raw_fastq/${s}_R1_trim.fastq.gz -p raw_fastq/${s}_R2_trim.fastq.gz \
    raw_fastq/${s}_R1p.fastq.gz raw_fastq/${s}_R2p.fastq.gz \
    > cutadapt_reports/${s}.log
```

## 5. BWA mapping + sort
```bash
s=$(sed -n "${SLURM_ARRAY_TASK_ID}p" metadata/sample_list.txt)
BAM=bwa_mapping/${s}_sorted.bam
bwa mem -t 8 ge/Aphal.fasta raw_fastq/${s}_R1_trim.fastq.gz raw_fastq/${s}_R2_trim.fastq.gz \
  | samtools view -F 4 -h -@ 2 - \
  | samtools sort -@ 4 -m 1g -O BAM -o ${BAM}.tmp -
mv ${BAM}.tmp $BAM
samtools index -@ 4 $BAM
```

## 6. Per-sample depth
```bash
s=$(sed -n "${SLURM_ARRAY_TASK_ID}p" metadata/sample_list.txt)
samtools depth bwa_mapping/${s}_sorted.bam > bwa_depth/${s}_depth.txt
awk -v x="$s" '{sum+=$3; n++} END{print x"\t"(n>0?sum/n:0)}' \
    bwa_depth/${s}_depth.txt > bwa_depth/${s}_depth_avg.txt
```

```bash
cat bwa_depth/*_depth_avg.txt > bwa_depth/all_depth_avg.txt
```

## 7. BAM thin to 30x
```bash
s=$(sed -n "${SLURM_ARRAY_TASK_ID}p" metadata/sample_list.txt)
INBAM=bwa_mapping/${s}_sorted.bam; OUTBAM=bwa_mapping/${s}_thin.bam
cov=$(grep -P "^${s}\t" bwa_depth/all_depth_avg.txt | cut -f2)
frac=$(awk -v c=$cov 'BEGIN{printf "%.6f", 30/c}')
if awk -v f=$frac 'BEGIN{exit !(f<1)}'; then
    samtools view -h -s $frac -b -@ 2 $INBAM > ${OUTBAM}.tmp
else
    cp $INBAM ${OUTBAM}.tmp
fi
mv ${OUTBAM}.tmp $OUTBAM
samtools index -@ 2 $OUTBAM
```

## 8. Read groups + index
```bash
s=$(sed -n "${SLURM_ARRAY_TASK_ID}p" metadata/sample_list.txt)
picard AddOrReplaceReadGroups \
    I=bwa_mapping/${s}_thin.bam O=bwa_mapping/${s}_info.bam \
    RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${s}
picard BuildBamIndex INPUT=bwa_mapping/${s}_info.bam
```

## 9. HaplotypeCaller to per-sample GVCF
```bash
s=$(sed -n "${SLURM_ARRAY_TASK_ID}p" metadata/sample_list.txt)
gatk --java-options '-Xmx10g' HaplotypeCaller \
    -R ge/Aphal.fasta \
    --native-pair-hmm-threads 4 \
    -I bwa_mapping/${s}_info.bam \
    -O bwa_haplotypecaller_gvcf/${s}.g.vcf.gz \
    -ploidy 2 -ERC GVCF
```

## 10a. Joint genotyping per contig (parallel)
```bash
contig=$(sed -n "${SLURM_ARRAY_TASK_ID}p" metadata/contig_list.txt)
GVCF=bwa_haplotypecaller_gvcf/Aphal_${contig}.g.vcf.gz
VCF=bwa_haplotypecaller_gvcf/Aphal_${contig}.vcf.gz
[[ -s $VCF && -s ${VCF}.tbi ]] && exit 0

vflag=$(awk '{print "--variant bwa_haplotypecaller_gvcf/"$1".g.vcf.gz"}' metadata/sample_list.txt | xargs)

gatk --java-options '-Xmx6g' CombineGVCFs \
    -R ge/Aphal.fasta $vflag \
    -L $contig -XL ge/Aphal.mask.bed \
    -O ${GVCF}.partial.g.vcf.gz
mv ${GVCF}.partial.g.vcf.gz $GVCF
[[ -f ${GVCF}.partial.g.vcf.gz.tbi ]] && mv ${GVCF}.partial.g.vcf.gz.tbi ${GVCF}.tbi

gatk --java-options '-Xmx6g' GenotypeGVCFs \
    -R ge/Aphal.fasta -V $GVCF \
    -O ${VCF}.partial.vcf.gz
mv ${VCF}.partial.vcf.gz $VCF
[[ -f ${VCF}.partial.vcf.gz.tbi ]] && mv ${VCF}.partial.vcf.gz.tbi ${VCF}.tbi
```

## 10b. Merge per-contig VCFs
```bash
> bwa_haplotypecaller_gvcf/Aphal_genotyped.list
while read contig; do
    echo bwa_haplotypecaller_gvcf/Aphal_${contig}.vcf.gz \
      >> bwa_haplotypecaller_gvcf/Aphal_genotyped.list
done < metadata/contig_list.txt

gatk --java-options '-Xmx12g' MergeVcfs \
    -I bwa_haplotypecaller_gvcf/Aphal_genotyped.list \
    -O bwa_haplotypecaller_gvcf/Aphal.vcf.gz
```

## 11. Hard filter + SNP-only + strict_SNP
```bash
IN=bwa_haplotypecaller_gvcf/Aphal.vcf.gz
HARD=bwa_haplotypecaller_gvcf/Aphal_hardfilter.vcf.gz
FILT=bwa_haplotypecaller_gvcf/Aphal_filtered.vcf.gz
FINAL=Aphal_SNPfiltered.vcf.gz
mkdir -p bwa_haplotypecaller_finalvcf

gatk --java-options '-Xmx6g' VariantFiltration \
    -R ge/Aphal.fasta -V $IN -O $HARD \
    --filter-name hard_filters \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"

gatk --java-options '-Xmx6g' SelectVariants \
    -R ge/Aphal.fasta -V $HARD -O $FILT \
    --exclude-filtered --select-type-to-include SNP -restrict-alleles-to BIALLELIC

# strict_SNP.py uses PyVCF3 (`pip install --user PyVCF3`)
python vcf_minitools/strict_SNP.py -i $FILT -o $FINAL
gatk --java-options '-Xmx2g' IndexFeatureFile -I $FINAL
```

## 12. Verify
```bash
bcftools stats Aphal_SNPfiltered.vcf.gz | grep -E '^SN'
# number of samples:    17
# number of records:    407558
# number of SNPs:       407558
# number of multiallelic SNP sites: 0
```

## Per-sample mean mapped depth (pre-thin) and final SNP-set depth

  | sample | mapped × | SNP-set × |
  |---|---:|---:|
  | 10221 | 55.6 | 24.6 |
  | 10222 | 91.4 | 25.0 |
  | 10223 | 73.0 | 25.6 |
  | 10225 | 63.7 | 24.4 |
  | 10229 | 31.0 | 24.3 |
  | 10231 | 26.5 | 20.9 |
  | 10232 | 36.1 | 24.5 |
  | 10281 | 106.6 | 24.1 |
  | 10294 | 103.1 | 23.5 |
  | 10334 | 53.0 | 21.9 |
  | 10708 | 103.2 | 22.8 |
  | 10709 | 97.0 | 22.5 |
  | 10715 | 128.3 | 22.5 |
  | 10716 | 80.5 | 22.3 |
  | 10717 | 87.5 | 22.3 |
  | 10718 | 67.4 | 21.9 |
  | 10721 | 127.4 | 22.3 |



#Determine number of Sequences 
grep -c "^>" signalp_results_mature.fasta




#Extracting sequences
/home/muzhinjin/tikafinal/ncbi-blast-2.15.0+/bin/makeblastdb -in genome.fasta -dbtype nucl -out genome_db
/home/muzhinjin/tikafinal/ncbi-blast-2.15.0+/bin/blastn -query SeptoriaITS.fasta -db genome_db -outfmt 6 -out SeptriaITS_hits.txt
/home/muzhinjin/tikafinal/samtools-1.19.2/samtools faidx Illuminaseptoriacontigs.fasta  NODE_259_length_7729_cov_2354.519937:3121-3611 > ITS1.fasta

# Septoria-PacBio-illumina-Hybrid-
Complete genome and comparative genomics
#FastQC Command:  
fastqc SLMR.L350_FDSW250082146-1r_1.fq.clean.gz SLMR.L350_FDSW250082146-1r_2.fq.clean.gz SLMR_350_1.fq.gz SLMR_350_2.fq.gz -o hastqc
3MultiQC Command:  
multiqc hastqc/ -o hastqc/

ern jobs submit --name=fastqc_multiqc --threads=4 --memory=16gb --hours=24 --module='fastqc multiqc' --command="fastqc SLMR.L350_FDSW250082146-1r_1.fq.clean.gz SLMR.L350_FDSW250082146-1r_2.fq.clean.gz SLMR_350_1.fq.gz SLMR_350_2.fq.gz -o hastqc/ && multiqc hastqc/ -o hastqc/"
#Trimmomatic
trimmomatic PE -threads 4 -trimlog NameLog SLMR.L350_FDSW250082146-1r_1.fq.clean.gz SLMR.L350_FDSW250082146-1r_2.fq.clean.gz SLMR_1_trimmed_paired.fq.gz SLMR_1_trimmed_unpaired.fq.gz SLMR_2_trimmed_paired.fq.gz SLMR_2_trimmed_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Extract PacBio Reads
# If PacBio reads are in a BAM file (SLMR.bam), convert to FASTQ:
samtools fastq SLMR.bam > SLMR_pacbio.fastq

# Hybrid Assembly Using Flye + Illumina Polishing
# Use long-read assembly with Flye:
flye --pacbio-raw SLMR_pacbio.fastq --out-dir flye_assembly --threads 16
ern jobs submit --name=Septoria --threads=8 --memory=128gb  --hours=48  --input="SLMR_pacbio.fastq" --module="flye" --command=flye -- --nano-raw SLMR_pacbio.fastq --out-dir flye_assembly --threads 16 --genome size 
Output: flye_assembly/assembly.fasta

# Polish Assembly with Illumina Reads
a. Index the assembly
bwa index flye_assembly/assembly.fasta

b. Map Illumina reads

ern jobs submit --name=Septoria --threads=8 --memory=128gb  --hours=48  --input="flye_assembly/assembly.fasta,SLMR_1_trimmed_paired.fq.gz,SLMR_2_trimmed_paired.fq.gz" --module="bwa" --command=bwa -- mem -t 16 flye_assembly/assembly.fasta SLMR_1_trimmed_paired.fq.gz SLMR_2_trimmed_paired.fq.gz > illumina_mapped.sam

samtools index illumina_mapped.bam
 /home/muzhinjin/tikafinal/samtools-1.19.2/samtools view -H illumina_mapped.sorted.bam view names of contigs
 grep "^>" polished.fasta

c. Polish with Pilon
pilon --genome flye_assembly/assembly.fasta --frags illumina_mapped.bam  --output polished --threads 8 --fix all
java -Xmx64G -jar pilon-1.24.jar --genome flye_assembly/assembly.fasta --frags illumina_mapped.sorted.bam --output polished

Run Pilon 2–3 times to improve qualityusing output from previous round

# Scaffold to Chromosome-Level Using Reference Genome
a. Align to reference genome
nucmer --prefix=align septoriarefgenome.fna polished.fasta
show-coords -rcl align.delta > align.coords

# Scaffold with RagTag (recommended for reference-guided scaffolding)
b. ragtag.py scaffold septoriarefgenome.fna polished.fasta -o ragtag_output

# Extract Mitochondria genome
/home/muzhinjin/tikafinal/ncbi-blast-2.15.0+/bin/makeblastdb -in ragtag_output/ragtagscaffold1.fasta -dbtype nucl -out genome_db
/home/muzhinjin/tikafinal/ncbi-blast-2.15.0+/bin/blastn -query MitochondriaCP099434Septorialinola.fasta -db genome_db -outfmt 6 -out mito_hits.tsv

# Extract unallignedsequences
awk '/^>/ {f=0} /^>contig_14_pilon|^>contig_21_pilon|^>contig_48_pilon|^>contig_50_pilon|^>contig_57_pilon/ {f=1} f' ragtag_output/ragtagscaffold1.fasta > uallignedcontigs.fasta

# Evaluate Final Assembly
quast ragtag_output/ragtag.scaffold.fasta -r septoriarefgenome.fna -o quast_report

# Repeat Annotation

RepeatModeler -database ragtag_output/ragtag.scaffold.fasta -pa 16 -LTRStruct

RepeatMasker -pa 16 -lib repeatmodeler.lib ragtag_output/ragtag.scaffold.fasta

#  Clean and Rename Chromosomes
seqkit rename ragtag_output/ragtag.scaffold.fasta > septoria_final.fasta

# Mask repeats (Funannotate can do this too)
funannotate mask -i septoria_final.fasta -o septoria_masked.fasta

# Prepare genome for annotation
funannotate clean -i septoria_final.fasta -o septoria.cleaned.fasta
funannotate sort -i septoria.cleaned.fasta -o septoria.sorted.fasta

# What to do with unaligned scaffolds?
✅ 1. Identify unaligned scaffolds
RagTag places these in a separate FASTA file:
ragtag_output/ragtag.scaffold.unplaced.fasta

🔎 2. Analyze unplaced scaffolds
a. Length & GC content
seqkit stats ragtag_output/ragtag.scaffold.unplaced.fasta
b. BLAST search
To determine similarity to known sequences:
blastn -query ragtag_output/ragtag.scaffold.unplaced.fasta \
       -db nt -outfmt 6 -max_target_seqs 1 -evalue 1e-5 -out unplaced_blast.txt
c. BUSCO analysis
Check if they contain core fungal genes:
busco -i ragtag_output/ragtag.scaffold.unplaced.fasta -l fungi_odb10 -m genome -o busco_unplaced
d. Repeat content
Repeatmasker ragtag_output/ragtag.scaffold.unplaced.fasta

# Use august 
augustus --species=aspergillus_nidulans ragtagscaffold1flteredfinalrenamed1.fasta > augustus_output.gff


🧩 3. Assign unplaced scaffolds to pseudochromosomes
If you want to include unplaced scaffolds in your final assembly for annotation:

Option 1: Concatenate them as "ChrUn"
Append to your scaffolded genome:
cat ragtag_output/ragtag.scaffold.fasta ragtag_output/ragtag.scaffold.unplaced.fasta > final_with_unplaced.fasta
Rename unplaced scaffolds to ChrUn_1, ChrUn_2, etc., using seqkit rename.

Option 2: Cluster using synteny or Hi-C (if available)
If you have Hi-C or linkage mapping data, use tools like:

SALSA, 3D-DNA, or ALLHiC to place them

Or use synteny-based clustering from related species

🧪 Optional: Tidy chromosomes before annotation
Use funannotate sort to reorder and relabel contigs in a consistent fashion (e.g., Chr1, Chr2, … ChrUn1, …).



1. Data Quality Control
bash
# For Illumina data
fastqc illumina_*.fastq.gz -o fastqc_results
multiqc fastqc_results -o multiqc_report

# For PacBio data
pbqc SLRM.bam LMR.bam -o pacbio_quality
2. Hybrid Genome Assembly
Option A
unicycler -1 illumina_R1.fastq.gz -2 illumina_R2.fastq.gz \
          -l pacbio.fastq.gz -o hybrid_assembly

Option B: Using SPAdes with hybrid option

spades.py --pe1-1 illumina_R1.fastq.gz --pe1-2 illumina_R2.fastq.gz \
          --pacbio pacbio.fastq.gz -o spades_hybrid_assembly
Option C: Using Flye + Polishing (for larger genomes)


# Polish with Illumina
bwa index flye_output/assembly.fasta
bwa mem flye_output/assembly.fasta illumina_R1.fastq.gz illumina_R2.fastq.gz | \
    samtools sort -o aligned.bam
pilon --genome flye_output/assembly.fasta --frags aligned.bam --output pilon_polished
3. Assembly Evaluation

quast.py hybrid_assembly/assembly.fasta -o assembly_quality
busco -i hybrid_assembly/assembly.fasta -l basidiomycota_odb10 -o busco_results -m genome
4. Chromosome-level Scaffolding



# Or using RagTag with reference
ragtag.py scaffold reference_genome.fasta hybrid_assembly/assembly.fasta -o ragtag_output
5. Annotation Pipeline (Enhanced with RNA-seq if available)
bash
funannotate predict -i ragtag_output/ragtag.scaffolds.fasta \
                    -o annotation_output \
                    -s "Species name" \
                    --busco_db basidiomycota \
                    --rna_bam rna_seq.bam \
                    --pasa_gff pasa_annotations.gff
6. Specialized Analyses
Effector prediction (with additional filters):
bash
# First run SignalP and TMHMM
signalp -f short -m signalp_results annotation_output/proteins.fa
tmhmm annotation_output/proteins.fa > tmhmm_results.txt

# Then run EffectorP with filtered set
EffectorP.py -i secreted_proteins.fa -o effector_results
Comprehensive CAZyme analysis:
bash
# Run dbCAN2
run_dbcan.py annotation_output/proteins.fa protein --out_dir dbcan_results --tools all

# Add CUPP classification
cupp.py -i dbcan_results/overview.txt -o cupp_results
Secondary metabolite clusters:
bash
# Run antiSMASH with full options
antismash --cb-general --cb-knownclusters --cb-subclusters --asf \
          --pfam2go --smcog-trees --genefinding-tool prodigal \
          -c 20 ragtag_output/ragtag.scaffolds.fasta \
          -o antismash_results
7. Integrated Visualization
bash
# Create JBrowse2 instance with all annotations
jbrowse create -i ragtag_output/ragtag.scaffolds.fasta \
               -g annotation_output/predict_results/genes.gff3 \
               -t annotation_output/predict_results/genes.fasta \
               --out genome_browser

# Add additional tracks
jbrowse add-track effector_results/effectors.gff --load copy
jbrowse add-track dbcan_results/overview.gff --load copy







Effector prediction
wget https://services.healthtech.dtu.dk/download/9782ed7b-1e4f-4227-9a28-b9abb9a0684e/signalp-5.0b.Linux.tar.gz
tar -xzf signalp-6.0.fast.tar.gz
signalp -fasta SLMR.pep.fa -format short -prefix SLMR_signalp
ls
cd signalp-5.0b
cd bin
echo 'export PATH=$PATH:/home/muzhinjin/Septoria/signalp-5.0b/bin' >> ~/.bashrc
source ~/.bashrc
/home/muzhinjin/Septoria/signalp-5.0b/bin/signalp -fasta SLMR.pep.fa -format short -prefix SLMR_signalp
signalp -fasta /home/muzhinjin/Septoria/SLMR.pep.fa -format short -prefix SLMR_signal
tmhmm SLMR.pep.fa > SLMR.tmhmm
which perl
nano /home/muzhinjin/Septoria/tmhmm-2.0c/bin/tmhmm
chmod +x /home/muzhinjin/Septoria/tmhmm-2.0c/bin/tmhmm
 tmhmm SLMR.pep.fa > SLMR.tmhmm
 echo 'export PATH=$PATH:/home/muzhinjin/Septoria/tmhmm-2.0c/bin' >> ~/.bashrc
 source ~/.bashrc
 tmhmm SLMR.pep.fa > SLMR.tmhmm
awk '!/^>/ {if (length($0) < 300) print prev; } {prev=$0;}' SLMR.pep.fa > SLMR.small_secreted.fa
grep -v '#' SLMR_signalp_summary.signalp5 | awk '$2 == "SP(Sec/SPI)" {print $1}' > SLMR.signalp_ids.txt
seqtk subseq SLMR.pep.fa SLMR.signalp_ids.txt > SLMR.signalp_proteins.fa
wget wget https://github.com/lh3/seqtk/archive/refs/tags/v1.4.tar.gz
tar -xzf v1.4.tar.gz
ls
 cd seqtk-1.4
 ls
  make
 ls
 pwd
 echo 'export PATH=$PATH:'$/home/muzhinjin/Septoria/seqtk-1.4 >> ~/.bashrc
 source ~/.bashrc
seqtk seq
/home/muzhinjin/Septoria/seqtk-1.4/seqtk seq
 /home/muzhinjin/Septoria/seqtk-1.4/seqtk subseq SLMR.pep.fa SLMR.signalp_ids.txt > SLMR.signalp_proteins.fa
  ls
   awk 'BEGIN {RS=">"; ORS=""} length($2) < 300 {print ">"$0}' SLMR.signalp_proteins.fa > SLMR.small_secreted1.fa


# Chromosome-Level Assembly Without Reference Genome or RNA-seq

Since you want to achieve chromosome-level assembly without a reference genome or RNA-seq data, we'll need to rely on alternative scaffolding methods. Here are the best approaches given your available data:

## 1. Initial High-Quality Hybrid Assembly

First, create the best possible contig assembly using both data types:

```bash
# Using MaSuRCA (excellent for fungal genomes)
masurca -g 50m -t 16 \
        -i clean_R1.fq.gz,clean_R2.fq.gz \
        -r pacbio_SLRM.fq,pacbio_LMR.fq \
        -o masurca_config.txt
./assemble.sh masurca_config.txt
```

## 2. Chromosome-Scale Scaffolding Options

### Option A: Hi-C Scaffolding (if you can generate Hi-C data)
```bash
# Process Hi-C data
juicer.sh -g 50m -d juicer_out -s MboI \
          -p chrom.sizes \
          -y restriction_sites.txt \
          masurca_assembly.fasta \
          hic_R1.fq.gz hic_R2.fq.gz

# Create chromosome scaffolds
3d-dna masurca_assembly.fasta juicer_out/merged_nodups.txt
```

### Option B: Optical Mapping (if you can generate Bionano data)
```bash
# Run Bionano Solve
solve -T 16 -f 0.3 -N 2 \
      -R bionano.cmap \
      -G masurca_assembly.fasta \
      -o bionano_scaffolds
```

### Option C: Using Your PacBio BAM Files for Ultra-Long Scaffolding
```bash
# Extract ultra-long reads from BAM files
samtools view -F 3844 SLRM.bam | awk 'length($10) > 50000' > ultra_long.fq

# Scaffold with LINKS
LINKS -f masurca_assembly.fasta -s empty.fof -k 20 -b ultra_long.fq
```

### Option D: De Novo Chromosome Reconstruction (for fungi)
```bash
# Using Salsa2 with PacBio data
samtools faidx masurca_assembly.fasta
bwa index masurca_assembly.fasta
bwa mem -x pacbio masurca_assembly.fasta pacbio_combined.fq | \
    samtools view -Sb - > pacbio_aligned.bam
samtools sort -o pacbio_aligned.sorted.bam pacbio_aligned.bam
samtools index pacbio_aligned.sorted.bam

run_pipeline.pl -Xmx100G -deNovoPipeline \
    -overlapTrimmer KEEP_INTEGRITY=false \
    -fastq1 pacbio_combined.fq \
    -output salsa_scaffolds \
    -iF 10 -iC 50
```

## 3. Chromosome Validation

```bash
# Check chromosome-length scaffolds
assembly-stats chromosome_assembly.fasta

# Validate with BUSCO
busco -i chromosome_assembly.fasta -l basidiomycota_odb10 -o busco_chrom -m genome

# Check synteny with related species (if any genes available)
mummer -mum -b -c ref_genome.fasta chromosome_assembly.fasta > mummer_out.txt
```

## 4. Chromosome Annotation

```bash
# Using Funannotate with chromosome-aware settings
funannotate predict -i chromosome_assembly.fasta \
                    -o chrom_annotation \
                    --busco_db basidiomycota \
                    --repeats2evm \
                    --optimize_augustus \
                    --ploidy 1 \
                    --organism fungus
```

## 5. Visualization of Chromosomes

```bash
# Create a Circos plot
circos -conf circos.conf -outputdir circos_output

# Generate JBrowse with chromosome view
jbrowse create -i chromosome_assembly.fasta \
               -g chrom_annotation/predict_results/genes.gff3 \
               --out chrom_browser \
               --tracks chromosomes
```

## Key Considerations Without Reference:

1. **Assembly Metrics**:
   - N50 should be dramatically higher than contig assembly
   - Expected chromosome count based on literature for your organism
   - Telomere-to-telomere completeness (look for telomeric repeats)

2. **Validation Methods**:
   - K-mer analysis: `kat comp -o kat_results -t 16 -m 21 reads.fq assembly.fasta`
   - Read mapping: `bowtie2-build assembly.fasta; bowtie2 -x index -1 R1.fq -2 R2.fq | samtools view -Sb - > mapped.bam`

3. **Chromosome Confirmation**:
   - Look for centromeric repeats
   - Check gene density distribution
   - Validate synteny blocks with related species

## Alternative Approach: Reference-Free Chromosome Reconstruction

If you can't generate additional data (Hi-C/Bionano), try this computational approach:

```bash
# 1. Identify core eukaryotic genes
cegma -g contigs.fasta -o cegma_out -T 16

# 2. Cluster contigs by sequence composition
sibeliaz -f contigs.fasta -o sibeliaz_out -n 16

# 3. Order contigs using long reads
ragtag.py order reference-free -u pacbio_aligned.sorted.bam -o ragtag_order

# 4. Merge into chromosomes
merge_contigs.py -i ordered_contigs.fasta -o chromosomes.fasta -c 8
```

Would you like me to:
1. Provide specific telomere identification methods for your organism group?
2. Explain how to estimate expected chromosome count without a reference?
3. Detail how to perform manual chromosome curation?
4. Suggest validation experiments to confirm chromosome structure?

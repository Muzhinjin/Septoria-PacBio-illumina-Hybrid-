# Loking for a file 
find / -name "getAnnoFasta.pl" 2>/dev/null
#Extract effectors frm Blast heat file
 awk -F'\t' '
NR==1 { print; next }
{
    if ($1 != prev) {
        if (NR > 2) print line
        line = $0
        prev = $1
        best = $3
    } else if ($3 > best) {
        line = $0
        best = $3
    }
}
END { print line }
' SlinoclaeffectrorsnoTHMM.tabular > best_hits.tabular

 pyani-plus list-runs --database Neo.db

faSize -detailed Septorialinicola.fna  > septorialinolica.tsv
 grep '>' ragtag.patch.fasta
 java -Xmx64G -jar pilon-1.24.jar --genome assembly_pilon1_round1.fasta  --frags illumina_paired_round2.bam --output assembly_pilon2 --threads 32
  /home/muzhinjin/tikafinal/bwa-0.7.17/bwa index illumina_paired_round2.bam
  /home/muzhinjin/tikafinal/samtools-1.19.2/samtools sort -@ 8 -o illumina_paired_round2.bam illumina_paired_round2.sam
  quast.py ragtag_output/ragtag.scaffold.fasta -r septoriarefgenome.fna -o quast_results
  ern jobs submit --name=Septoriabusco --threads=32 --memory=128gb  --hours=48  --input="Septoriagenomeassempledfinal_sorted.fasta" --module="busco/1.0_88de6b8" --command=busco -- -i Septoriagenomeassempledfinal_sorted.fasta  -l dothideomycetes_odb10 -m genome -o busco_out -c 32
  
# Effectors
python /home/muzhinjin/Septoria/septriagenomes/EffectorP/Scripts/EffectorP.py -i SLML2signalp_only.fa -o effectorp_results.txt
python EffectorP.py -i /home/muzhinjin/Septoria/SgnalP/Septoriaglycinesignalp_only.fa -o effectorp_results.txt
# extract the IDs predicted as effectors using awk or Python. For example:
awk '$2=="Effector" {print $1}' Septorialincolaeffectorp_results.txt > effector_ids.txt
grep -E "Apoplastic|Cytoplasmic" Slinicola_results.txt | awk '{print $1}' > Slinolaeeffector_IDs.txt
# Extract sequences from your FASTA
seqkit grep -f effector_ids.txt Septorialincolasignalp_only.fasta > Septorialincola_effectors.fa
seqkit grep -f  Slinolaeeffector_IDs.txt /home/muzhinjin/Septoria/SgnalP/Septorialincolasignalp_only.fasta > Slinicolaeffectorp_only.fasta
seqkit seq -m 1 -M 300 -g Neopestalotiopsis_rosae_1902predictedeffectors.fasta -o Neopestalotiopsis_rosae_1902predictedeffectors_under300.fa


# Extractthe top 4 contigs
seqkit sort -l Finalassemplyragtag.scaffold.fasta | head -n 4 > top_contigs.fasta
# SYN
ern jobs submit --name=Septorisyny --threads=32 --memory=128gb  --hours=48  --input="*.gbff" --module="syny/1.0_1294505" --command=run_syny.pl -- -a *.gbff -o finaloutput directory


Module load cluster/hpc
augustus --species=botrytis_cinerea --protein=on --gff3=on --stopCodonExcludedFromCDS=false combinedallclassfiedsorted2ffinnaotatecleaned.fasta > augustus_septoria.gff3
# Mantain the names
augustus --species=botrytis_cinerea --protein=on --gff3=on --stopCodonExcludedFromCDS=false --uniqueGeneId=true combinedallclassfiedsorted2ffinnaotatecleaned.fasta > augustus_septoria.gff3
/home/muzhinjin/miniconda3/envs/genome_assembly/bin/getAnnoFasta.pl SLM2augustus_septoria.gff3


minimap2 -ax map-pb -t 8 pilonround2_polished.fasta SLMR_pacbio.fasta > polishedalr.sam
/home/muzhinjin/tikafinal/samtools-1.19.2/samtools view -@ 8 -bS polishedalr.sam > polishedlr_unsorted.bam
 /home/muzhinjin/tikafinal/samtools-1.19.2/samtools sort -@ 8 -o polishedlr_sorted.bam > polishedlr_unsorted.bam
  /home/muzhinjin/tikafinal/samtools-1.19.2/samtools index polishedlr_unsorted.bam
  /home/muzhinjin/tikafinal/samtools-1.19.2/samtools sort -@ 8 -o polishedlr_sorted.bam > polishedlr_unsorted.bam
  /home/muzhinjin/tikafinal/samtools-1.19.2/samtools view -@ 8 -bS polishedalr.sam > polishedlr_unsorted.bam
  /home/muzhinjin/tikafinal/samtools-1.19.2/samtools sort -@ 8 -o polishedlr.bam polishedlr_unsorted.bam
  /home/muzhinjin/tikafinal/samtools-1.19.2/samtools index polishedlr.bam
  /home/muzhinjin/tikafinal/samtools-1.19.2/samtools idxstats polishedlr.bam > polisheddepth.tsv
  cat polisheddepth.tsv
#Rename
sed -i 's/Chr\([0-9]\+\)/Chr0\1/g' combinedallclassfiedsorted2f.fasta.gbff

# Augustus
# Extraxct prtein sequences

/home/muzhinjin/miniconda3/envs/genome_assembly/bin/getAnnoFasta.pl SL2augustus_septoria.gff
# make sure they are the same name 
awk 'BEGIN{FS=OFS="\t"} !/^#/{$1=prevseq} /^>/ {sub(/^>/,"",$1); prevseq=$1; next} 1' input.fasta augustus_output.gff > restored_names.gff


 # SignalP
 /home/muzhinjin/miniconda3/envs/genome_assembly/bin/getAnnoFasta.pl Septoriliniclaaugustus_septoria.gff3
signalp -fasta Septoriaglycinesaugustus_septoria3.aa -format short -org euk -prefix Septoriaglycine_signalp
awk '$2 ~ /^SP/ {print $1}' Septoriaglycine_signalp_summary.signalp5 > Septoriaglycinessp_ids.txt
seqkit grep -f Septoriaglycinessp_ids.txt Septoriaglycinesaugustus_septoria3.aa > Septoriaglycinesignalp_only.fa

# Effectors
python /home/muzhinjin/Septoria/septriagenomes/EffectorP/Scripts/EffectorP.py -i SLML2signalp_only.fa -o effectorp_results.txt
 
 git clone https://github.com/JanaSperschneider/EffectorP-3.0.git
 cd EffectorP-3.0
 unzip weka-3-8-4.zip
 python EffectorP.py -i /home/muzhinjin/Septoria/SgnalP/SL2signalp_only.fa -o SL2effectorp_results.txt
# extract the IDs predicted as effectors using awk or Python. For example:
awk '$2=="Effector" {print $1}' Septorialincolaeffectorp_results.txt > effector_ids.txt
grep -E "Apoplastic|Cytoplasmic" effectorp_results.txt | awk '{print $1}' > Septoriaglycineeeffector_IDs.txt
# Extract sequences from your FASTA
seqkit grep -f effector_ids.txt Septorialincolasignalp_only.fasta > Septorialincola_effectors.fa
seqkit grep -f Septoriaglycineeeffector_IDs.txt /home/muzhinjin/Septoria/SgnalP/Septoriaglycinesignalp_only.fa > Septoriaglycineeffectorp_only.fasta

seqkit stats Septoriaglycinesignalp_only.fa
# Determine number of Sequences 
grep -c "^>" signalp_results_mature.fasta

# Extracting sequences
/home/muzhinjin/tikafinal/ncbi-blast-2.15.0+/bin/makeblastdb -in genome.fasta -dbtype nucl -out genome_db
/home/muzhinjin/tikafinal/ncbi-blast-2.15.0+/bin/blastn -query SeptoriaITS.fasta -db genome_db -outfmt 6 -out SeptriaITS_hits.txt
/home/muzhinjin/tikafinal/samtools-1.19.2/samtools faidx Illuminaseptoriacontigs.fasta  NODE_259_length_7729_cov_2354.519937:3121-3611 > ITS1.fasta
cat 
awk '/^>/ {print; next} {seq=$0; rev=""; for(i=length(seq);i!=0;i--) { base=substr(seq,i,1); comp=toupper(base); if(comp=="A") comp="T"; else if(comp=="T") comp="A"; else if(comp=="C") comp="G"; else if(comp=="G") comp="C"; rev=rev comp } print rev }' sepextractedtelf.fasta > sepextractedtelfrc.fasta

# Septoria-PacBio-illumina-Hybrid-
Complete genome and comparative genomics
# FastQC Command:  
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
ern jobs submit --name=Septoriabusco --threads=32 --memory=128gb  --hours=48  --input="Finalassemplyragtag.scaffold.fasta" --module="busco/1.0_88de6b8" --command=busco -- Finalassemplyragtag.scaffold.fasta -l dothideomycetes_odb10 -m genome -o busco_out -c 32

# TEs
conda install -c bioconda trnascan-se
tRNAscan-SE -o trnas.out -f trnas_struct.out $FASTA

# Polish Assembly with Illumina Reads
a. Index the assembly
bwa index flye_assembly/assembly.fasta

b. Map Illumina reads

ern jobs submit --name=Septoria --threads=8 --memory=128gb  --hours=48  --input="flye_assembly/assembly.fasta,SLMR_1_trimmed_paired.fq.gz,SLMR_2_trimmed_paired.fq.gz" --module="bwa" --command=bwa -- mem -t 16 flye_assembly/assembly.fasta SLMR_1_trimmed_paired.fq.gz SLMR_2_trimmed_paired.fq.gz > illumina_mapped.sam

samtools index illumina_mapped.bam
 /home/muzhinjin/tikafinal/samtools-1.19.2/samtools view -H illumina_mapped.sorted.bam view names of contigs
 grep "^>" polished.fasta
 awk '{print "chr -",$1,$1,"0",$2,"chr"$1}' finalscaffolds.sizes > finalkaryotype.txt

c. Polish with Pilon
pilon --genome flye_assembly/assembly.fasta --frags illumina_mapped.bam  --output polished --threads 8 --fix all
java -Xmx64G -jar pilon-1.24.jar --genome flye_assembly/assembly.fasta --frags illumina_mapped.sorted.bam --output polished

Run Pilon 2–3 times to improve qualityusing output from previous round

/home/muzhinjin/tikafinal/bwa-0.7.17/bwa mem -t 16 pillon_round1.fasta SLMR_1_trimmed_paired.fq.gz SLMR_350_2.fq | samtools sort -o illumina_mapped_round2.sorted.bam


# Scaffold to Chromosome-Level Using Reference Genome
a. Align to reference genome
nucmer --prefix=align septoriarefgenome.fna polished.fasta
show-coords -rcl align.delta > align.coords

# Scaffold with RagTag (recommended for reference-guided scaffolding)
b. ragtag.py scaffold septoriarefgenome.fna polished.fasta -o ragtag_output

/home/muzhinjin/tikafinal/bwa-0.7.17/bwa index pilon_polished.fasta
 /home/muzhinjin/tikafinal/bwa-0.7.17/bwa bwa mem -t 16 pilon_polished.fasta SLMR.L350_FDSW250082146-1r_1clean.fq reads_R2.fastq.gz > pilonround2.sam
  /home/muzhinjin/tikafinal/bwa-0.7.17/bwa mem -t 16 pilon_polished.fasta SLMR.L350_FDSW250082146-1r_1clean.fq reads_R2.fastq.gz > pilonround2.sam
  /home/muzhinjin/tikafinal/samtools-1.19.2/samtools view -bS pilonround2.sam > pilonround2.bam
  /home/muzhinjin/tikafinal/samtools-1.19.2/samtools sort -@ 16 -o pilonround2.sorted.bam pilonround2.bam 
  /home/muzhinjin/tikafinal/samtools-1.19.2/samtools index pilonround2.sorted.bam
  java -Xmx64G -jar pilon-1.24.jar --genome pilon_polished.fasta --frags pilonround2.sorted.bam   --output pilonround2_polished --threads 16

# Extract Mitochondria genome
/home/muzhinjin/tikafinal/ncbi-blast-2.15.0+/bin/makeblastdb -in ragtag_output/ragtagscaffold1.fasta -dbtype nucl -out genome_db
/home/muzhinjin/tikafinal/ncbi-blast-2.15.0+/bin/blastn -query MitochondriaCP099434Septorialinola.fasta -db genome_db -outfmt 6 -out mito_hits.tsv

# Extract unallignedsequences
awk '/^>/ {f=0} /^>contig_14_pilon|^>contig_21_pilon|^>contig_48_pilon|^>contig_50_pilon|^>contig_57_pilon/ {f=1} f' ragtag_output/ragtagscaffold1.fasta > uallignedcontigs.fasta
awk '/^>/{f=($0 ~ /CP099434.1_RagTag_pilon_pilon/)} f' pilonround2_polished.fasta > CP099434mitochonria.fasta

# Evaluate Final Assembly
quast ragtag_output/ragtag.scaffold.fasta -r septoriarefgenome.fna -o quast_report

# Repeat Annotation
 RepeatMasker -pa 8 -species "aspergillus_fumigatus" sortedragtagscaffold1flteredfinalrenamed1.fasta

RepeatModeler -database ragtag_output/ragtag.scaffold.fasta -pa 16 -LTRStruct

RepeatMasker -pa 16 -lib repeatmodeler.lib ragtag_output/ragtag.scaffold.fasta

# Determine the sizes
faSize -detailed scaffolds_100kb.fasta > scaffolds.sizes
cat polished2scaffolds.sizes
python chrom_stats.py > chrom_stats.tsv
# Repeat ends
 grep -B1 -A1 -E "TTAGGG|CCCTAA" ragtagscaffold1flteredfinalrenamed1.fasta
 ./telomere_check.py

#  Clean and Rename Chromosomes
seqkit rename ragtag_output/ragtag.scaffold.fasta > septoria_final.fasta

# Remove less than 10000 
seqkit seq -m 1000 ragtagscaffold1flteredfinalrenamed1.fasta > scaffolds_1kb.fasta

# Mask repeats (Funannotate can do this too)
funannotate mask -i septoria_final.fasta -o septoria_masked.fasta

# Chromsme stats
python chrom_stats.py > chrom_stats.tsv

# Prepare genome for annotation
funannotate clean -i septoria_final.fasta -o septoria.cleaned.fasta
ern jobs submit --name=Septoriafunantatecleanglycine --threads=32 --memory=128gb  --hours=48  --input="Septoriaglycines.fna" --module="funannotate/1.0_1e48052" --command=funannotate -- clean i Septoriaglycines.fna -o funnantatecleanSeptoriaglycines.fasta
funannotate sort -i septoria.cleaned.fasta -o septoria.sorted.fasta

# Length & GC content
seqkit stats ragtag_output/ragtag.scaffold.unplaced.fasta

# Antismash
ern jobs submit --name=ProkkaSeptoriaglycine --threads=32 --memory=128gb  --hours=48  --input="Septoriaglycines.fna" --module="prokka/1.0_d8de48f" --command=prokka -- Septoriaglycines.fna --outdir prokka_out --prefix mygenome

# BLAST search
To determine similarity to known sequences:
blastn -query ragtag_output/ragtag.scaffold.unplaced.fasta \
       -db nt -outfmt 6 -max_target_seqs 1 -evalue 1e-5 -out unplaced_blast.txt
# BUSCO analysis
Check if they contain core fungal genes:
busco -i ragtag_output/ragtag.scaffold.unplaced.fasta -l fungi_odb10 -m genome -o busco_unplaced
# Repeat content
Repeatmasker ragtag_output/ragtag.scaffold.unplaced.fasta

# Use august 
augustus --species=aspergillus_nidulans ragtagscaffold1flteredfinalrenamed1.fasta > augustus_output.gff
# Assign unplaced scaffolds to pseudochromosomes
include unplaced scaffolds in your final assembly for annotation:

# Concatenate them as "ChrUn"
Append to your scaffolded genome:
cat ragtag_output/ragtag.scaffold.fasta ragtag_output/ragtag.scaffold.unplaced.fasta > final_with_unplaced.fasta
Rename unplaced scaffolds to ChrUn_1, ChrUn_2, etc., using seqkit rename.

# Tidy chromosomes before annotation
Use funannotate sort to reorder and relabel contigs in a consistent fashion (e.g., Chr1, Chr2, … ChrUn1, …).



1. Data Quality Control

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

Prep: convert BAM → FASTQ (only if you don’t really have SLMR_pacbio.fastq)
# only run if SLMR_pacbio.fastq does not already exist
samtools fastq -@ 16 SLMR.bam > SLMR_pacbio.fastq

#  map long reads to contigs
minimap2 -x map-pb -t 16 Illuminaseptoriacontigs.fasta SLMR_pacbio.fastq > long_vs_contigs.paf

# round 1 racon
# Install racon using conda
conda create -n racon_env -c bioconda racon
conda activate racon_env

racon -t 16 SLMR_pacbio.fastq long_vs_contigs.paf Illuminaseptoriacontigs.fasta > contigs.racon1.fasta
#  round 2 racon
minimap2 -x map-pb -t 16 contigs.racon1.fasta SLMR_pacbio.fastq > 2long_vs_contigs.paf
racon -t 16 SLMR_pacbio.fastq 2long_vs_contigs.paf contigs.racon1.fasta > contigs.racon2.fast

#  correct misassemblies relative to reference
ragtag.py correct -t 16 septoriarefgenome.fna contigs.racon2.fasta -o ragtag_correct

# 2scaffold corrected assembly to reference chromosomes
ragtag.py scaffold -t 16 septoriarefgenome.fna ragtag_correct/ragtag.correct.fasta -o ragtag_scaffold
# outputs:
# ragtag_scaffold/ragtag.scaffold.fasta  -> main scaffolded assembly
# ragtag_scaffold/ragtag.scaffold.agp
# ragtag_scaffold/ (logs)

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

pyani-plus anim files/ --database Septoria.db --create-db --name "Septria_ANIm"
pyani-plus list-runs --database Septoria.db
 pyani-plus export-run --database Neo.db --outdir walkthrough_output/ --run-id 1
pyani-plus plot-run --database Neo.db --outdir walkthrough_output/ --run-id 1

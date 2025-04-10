# Septoria-PacBio-illumina-Hybrid-
Complete genome and comparative genomics

#FastQC Command:  
fastqc SLMR.L350_FDSW250082146-1r_1.fq.clean.gz SLMR.L350_FDSW250082146-1r_2.fq.clean.gz SLMR_350_1.fq.gz SLMR_350_2.fq.gz -o hastqc
3MultiQC Command:  
multiqc hastqc/ -o hastqc/

ern jobs submit --name=fastqc_multiqc --threads=4 --memory=16gb --hours=24 --module='fastqc multiqc' --command="fastqc SLMR.L350_FDSW250082146-1r_1.fq.clean.gz SLMR.L350_FDSW250082146-1r_2.fq.clean.gz SLMR_350_1.fq.gz SLMR_350_2.fq.gz -o hastqc/ && multiqc hastqc/ -o hastqc/"
#Trimmomatic
trimmomatic PE -threads 4 -trimlog NameLog SLMR.L350_FDSW250082146-1r_1.fq.clean.gz SLMR.L350_FDSW250082146-1r_2.fq.clean.gz SLMR_1_trimmed_paired.fq.gz SLMR_1_trimmed_unpaired.fq.gz SLMR_2_trimmed_paired.fq.gz SLMR_2_trimmed_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

ragtag.py patch SLMRcontigs.fasta Illuminaseptoriacontigs.fasta









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



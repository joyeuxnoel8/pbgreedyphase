#!/usr/bin/env bash

samtools mpileup -q0 -Q 0 $1 | /net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/Phasing/MpileupToFreq.py  /dev/stdin | $PBS/Phasing/PrintHetFreq.py 10 > assembly.consensus.nuqfreq

samtools view -h $1|  /net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/Phasing/readToSNVList  --nft assembly.consensus.nucfreq --sam /dev/stdin  --ref assembly.consensus.fasta --out assembly.consensus.fragments --minFraction 0.01 --minCoverage 5

/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/Phasing/FragmentsToPhaseGraph.py assembly.consensus.fragments --vcf assembly.consensus.fragments.vcf --contig assembly.consensus.fasta  --minAlleleCov 5  --minOverlap 2


samtools view -h $1 | /net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/Phasing/partitionByPhasedSNVs  --vcf assembly.consensus.fragments.1.vcf  --ref assembly.consensus.fasta --h1 h1.sam --h2 h2.sam --sam /dev/stdin 


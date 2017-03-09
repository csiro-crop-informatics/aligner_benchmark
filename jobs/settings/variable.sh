#!/bin/bash

# Output root directory
OUTPUT_ROOT_PATH="/project/itmatlab/aligner_benchmark/tool_results"

# Scripts directory
SCRIPTS_PATH="/project/itmatlab/aligner_benchmark/scripts"

#
BIOKANGA_PATH="/home/ste69r/bin"
BIOKANGA_ALIGN_CMD="biokanga align"
BIOKANGA_INDEX_CMD="biokanga index"

# CUTADAPT
CUTADAPT_PATH="/home/baruzzog/.local/bin"
CUTADAPT_CMD="cutadapt"

# BLAST
BLAST_PATH="/home/baruzzog/tools/aligner/blast/ncbi-blast-2.2.31+/bin"
BLAST_ALIGN_CMD="blastn"
BLAST_INDEX_CMD="makeblastdb"

# Bowtie
BOWTIE_PATH="/home/baruzzog/tools/aligner/bowtie/bowtie-1.1.2"
BOWTIE_ALIGN_CMD="bowtie"
BOWTIE_INDEX_CMD="bowtie-build"


# Bowtie2
BOWTIE2_PATH="/home/baruzzog/tools/aligner/bowtie2/bowtie2-2.2.5"
BOWTIE2_ALIGN_CMD="bowtie2"
BOWTIE2_INDEX_CMD="bowtie2-build"


# BWA
BWA_PATH="/home/baruzzog/tools/aligner/bwa/bwa-0.7.12"
BWA_ALIGN_CMD="bwa aln"
BWAMEM_ALIGN_CMD="bwa mem"
BWA_INDEX_CMD="bwa index"


# Contextmap 2
CONTEXTMAP2_PATH="/home/baruzzog/tools/aligner/contextmap2/ContextMap_v2.6.0"
LASTVERSION_CONTEXTMAP2_PATH="/home/baruzzog/tools/aligner-lastversion/contextmap2/ContextMap_v2.7.6"
CONTEXTMAP2_ALIGN_CMD="ContextMap_v2.6.0.jar mapper"
LASTVERSION_CONTEXTMAP2_ALIGN_CMD="ContextMap_v2.7.6.jar mapper"


# CRAC
CRAC_PATH="/home/baruzzog/test_installs/bin"
LASTVERSION_CRAC_PATH="/home/baruzzog/test_installs_2/bin"
CRAC_LIB="/home/baruzzog/test_installs/usr/local/lib"
LASTVERSION_CRAC_LIB="/home/baruzzog/test_installs_2/usr/local/lib"
CRAC_INDEX_CMD="crac-index"
LASTVERSION_CRAC_INDEX_CMD="crac-index"
CRAC_ALIGN_CMD="crac"
LASTVERSION_CRAC_ALIGN_CMD="crac"



# GSNAP
GSNAP_PATH="/home/baruzzog/tools/aligner/gsnap/gmap-2015-09-29/bin"
LASTVERSION_GSNAP_PATH="/home/baruzzog/tools/aligner-lastversion/gsnap/gmap-2016-06-30/bin"
GSNAP_INDEX_CMD="gmap_build"
GSNAP_ALIGN_CMD="gsnap"


# Hisat
HISAT_PATH="/home/baruzzog/tools/aligner/hisat/hisat-master"
HISAT_INDEX_CMD="hisat-build"
HISAT_ALIGN_CMD="hisat"


# Hisat2
HISAT2_PATH="/home/ste69r/hisat2/hisat2-2.0.5"
LASTVERSION_HISAT2_PATH="/home/ste69r/hisat2/hisat2-2.0.5"
HISAT2_INDEX_CMD="hisat2-build"
HISAT2_ALIGN_CMD="hisat2"


# Mapsplice 2
MAPSPLICE2_PATH="/home/baruzzog/tools/aligner/mapsplice2/MapSplice-v2.2.0"
LASTVERSION_MAPSPLICE2_PATH="/home/baruzzog/tools/aligner-lastversion/mapsplice2/MapSplice_3.0.1_beta-master"
MAPSPLICE2_INDEX_CMD=""
MAPSPLICE2_ALIGN_CMD="mapsplice.py"


# Novoalign
NOVOALIGN_PATH="/home/baruzzog/tools/aligner/novoalign/novocraft"
LASTVERSION_NOVOALIGN_PATH="/home/baruzzog/tools/aligner-lastversion/novoalign/novocraft"
NOVOALIGN_INDEX_CMD="novoindex"
NOVOALIGN_ALIGN_CMD="novoalign"

# Useq (used by novoalign)
USEQ_PATH="/home/baruzzog/tools/aligner/useq/USeq_8.9.5/Apps"
USEQ_MAKE_TRANSCRIPTOME_CMD="MakeTranscriptome"
USEQ_MASK_EXONS_IN_FASTA_FILES_CMD="MaskExonsInFastaFiles"
USEQ_SAM_TRANSCRIPTOME_PARSER_CMD="SamTranscriptomeParser"


# olego
OLEGO_PATH="/home/baruzzog/tools/aligner/olego/olego.src.v1.1.6"
OLEGO_INDEX_CMD="olegoindex"
OLEGO_ALIGN_CMD="olego"

# Mono (used to run osa)
MONO_PATH="/home/baruzzog/tools/aligner/mono/mono/bin"
MONO_CMD="mono"


# OSA
OSA_PATH="/home/baruzzog/tools/aligner/osa/OSAv4.1.1.1"
OSA_INDEX_CMD="osa.exe"
OSA_ALIGN_CMD="osa.exe"
OSA_CONF_PATH="/project/itmatlab/aligner_benchmark/jobs/osa"
OSA_CONF_FILE="osa_parameters.ini"


# RNAsequel
RNASEQUEL_PATH="/home/baruzzog/tools/aligner/rnasequel/RNASequel-0.9b/src"
RNASEQUEL_INDEX_CMD="rnasequel"
RNASEQUEL_ALIGN_CMD="rnasequel"


# RUM
RUM_PATH="/home/baruzzog/tools/aligner/rum/rum/bin"
RUM_INDEX_CMD="rum_indexes"
RUM_ALIGN_CMD="rum_runner"


# SOAP
SOAP_PATH="/home/baruzzog/tools/aligner/soap/soap2.21release"
SOAP_INDEX_CMD="2bwt-builder"
SOAP_ALIGN_CMD="soap"
SOAP_SCRIPT_PATH="/home/baruzzog/itmat/aligner_benchmark/scripts"
SOAP_SCRIPT_CONVERT_TO_SAM="soap2sam.pl"


# SOAPsplice
SOAPSPLICE_PATH="/home/baruzzog/tools/aligner/soapsplice/SOAPsplice-v1.10/bin"
SOAPSPLICE_INDEX_CMD="2bwt-builder"
SOAPSPLICE_ALIGN_CMD="soapsplice"


# Spliceseq
SPLICESEQ_PATH="/home/baruzzog/tools/aligner/spliceseq/spliceseq"
SPLICESEQ_INDEX_CMD=""
SPLICESEQ_ALIGN_CMD="SpliceSeqAnalyze.jar"


# STAR
STAR_PATH="/home/baruzzog/tools/aligner/star/STAR-STAR_2.5.0a/source"
LASTVERSION_STAR_PATH="/home/baruzzog/tools/aligner-lastversion/star/STAR-2.5.2a/source"
STAR_INDEX_CMD="STAR"
STAR_ALIGN_CMD="STAR"


# subread-subjunc
#SUBREAD_PATH="/home/baruzzog/tools/aligner/subread/subread-1.5.0-source/bin"
SUBREAD_PATH="/home/baruzzog/tools/aligner/subread/subread-1.5.0-Linux-x86_64/bin"
LASTVERSION_SUBREAD_PATH="/home/baruzzog/tools/aligner-lastversion/subread/subread-1.5.0-p3-source/bin"
SUBREAD_INDEX_CMD="subread-buildindex"
SUBREAD_ALIGN_CMD="subread-align"
SUBJUNC_ALIGN_CMD="subjunc"


# Tophat2
TOPHAT2_PATH="/home/baruzzog/tools/aligner/tophat2/tophat-2.1.0.Linux_x86_64"
LASTVERSION_TOPHAT2_PATH="/home/baruzzog/tools/aligner-lastversion/tophat2/tophat-2.1.1.Linux_x86_64"
TOPHAT2_INDEX_CMD=""
TOPHAT2_ALIGN_CMD="tophat2"





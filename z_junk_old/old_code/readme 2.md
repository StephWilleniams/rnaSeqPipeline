
# Processing fastQ data

## The kit is a SMART-Seq mRNA LP (with UMIs)

## 0. Download the rnaSeq data from UCR HPC

`wget -r -e robots=off --no-parent --http-user='xcge' --http-password='Hd4ej$ghGwo!jd$G!' http://cluster.hpcc.ucr.edu/~genomics/xcge`

Warning: the above is heavily monitored/rate-limited. Repeated use will break the link.
It is also date-sensitive (expected deletion date mid-July 2026)

Data from this wget needs to be stored in `d_data/`

## 0. Alternative gDrive download method

The following can be executed

`FILEID="1lsUj0Mlbs9mGuoLaoAo6LQXD4CkQcqev"`
`FILENAME="fastq_data"`
`wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id='$FILEID -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id="$FILEID -O $FILENAME && rm -rf /tmp/cookies.txt`

Data from this wget should be stored in d_data/

## 0. Download the reference genome + annotation

Download the complete genomic FASTA sequence:
Download the corresponding RefSeq gene annotations (GFF format):

`wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz`

`wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gff.gz`

Data from these wget are stored in `d_data/refGenome/`

## 0. Set up the Conda virtual environment to get the libraries ect

The following sequence of commands should be run to set up the workspace for running the scripts.

`conda create -n bio_qc`
`conda activate bio_qc`

`conda config --add channels defaults`
`conda config --add channels conda-forge`
`conda config --add channels bioconda`
`conda config --set channel_priority strict`

`conda install fastqc multiqc fastp STAR`

to export run:
`conda env export --no-builds > bioqc_environment.yml`

## 0. ALTERNATIVE set up for env W/ bioqc_environment.yml

conda env create -f bioqc_environment.yml

## 1. Run fastqc and multiqc on the files to generate the reports

In this step we check the quality of the rnaSeq data using fastqc.

Run command:
`chmod +x c_code/1_fastqc_check_raw.sh`
to convert "1_fastqc_check_raw.sh" into an executable shell (.sh) file
then run:
`./c_code/1_fastqc_check_raw.sh`

The results of this should produce a folder `o_outputs/`
in this there will be two subfolders `qc_reports` and `qc_reports_multi` with the information on the fastqc and multiqc checks.

## 2. Run fastp for the Takara SMART-seq mRNA LP (UNIs) kit

In this step, we are trimming the fastq files from the rnaSeq to remove the initial strings (attached during processing) as well as removing erroneous sections using a quality trim.

Run command:
`chmod +x c_code/trim_raw.sh`
to convert "chmod +x c_code/2_trim_raw.sh" into an executable shell (.sh) file
then run:
`./c_code/2_trim_raw.sh`

This process also uses the umi_tools library to extract the umis at the starts of the sequences, which will be needed later on for deduplification.

Notes on possible issues:

- If used with a large number of threads per job this can lead to memory issues which will cause the code to cease without warning (it will appear to continue to run).
- It seems like smaller numbers of threads per job partially solves the issue (i.e., 2).
- The issue isn't repeatable (i.e., different files will/won't process).
- The job can also be done piecewise if the logic in the script to remove prior processing runs is removed (it currently stores old versions).
- The job is expected to be deterministic, so reruns should generate the same trims, as such the job performed piecewise should be sufficient.

This job has a bit of a runtime, in a range of 300-1100 seconds for each file to be trimmed.
The first output of this job is stored in `o_outputs/extracted_umis`.
The final output of this job is stored in `o_outputs/trimmed_fastqs`, it contains html/json reports on the timming process, as the trimmed fastq files.

## 3. Rerun the fastQC on the trimmed data to check

To check the trim data has worked correctly we rerun FastQC on the trimmed data, this step is identical to step 1, except using the trimmed data.

This is done by running:
`./c_code/3_fastqc_check_trim.sh`

the output of this will be stored in another pair of subfolders in `o_outputs/`, labelled by the additional indicator `_trimmed` containing the fastqc data on the trimmed fastq files.

## 4. Generate the STAR index library for the sequence alignment

We need to generate the STAR index from the reference genome data to do the alignment.
To do this first the ref files downloaded from step 0 need to be uncompressed.
To do this run the command:

`gunzip d_data/refGenome/GCF_000001635.27_GRCm39_genomic.fna.gz`
`gunzip d_data/refGenome/GCF_000001635.27_GRCm39_genomic.gff.gz`

Note: the above removes the .gz files, and so only needs to be done once.

After this change permissions to execute run the script:

`chmod +x c_code/4_make_star_index.sh`
`./c_code/4_make_star_index.sh`

Note:

- This job has a bit of a runtime (~40 minutes), and only needs to be performed once.
- Additionally, this job is computationally intensive and so should be expected to initate the computers fan.
- Expect issues if you have limited access to memory, the process requires 30-32GB (since the whole seq is stored). It can be done with HISAT2 using ~5GB of memory, however this is expected to slow down the alignment and possibly lead to unmapped regions where there are long intron distances involved.

The output is placed in the `d_data/refGenome/mouse_star_index` for use in step 5, as well as log data about the generation process.

## 5. Perform the sequence alignment

CONDA_SUBDIR=osx-64 mamba install -y bioconda::star=2.7.10b

After this we align the sequences to the reference genome.

to do this, as before:

`chmod +x c_code/5_align_to_STAR_index.sh`

The outputs are then stored in `o_outputs/star_alignments`

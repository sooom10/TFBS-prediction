# TFBS-prediction

# Download data

    # Experiment bams
    wget --no-check-certificate https://www.encodeproject.org/files/ENCFF698HBZ/@@download/ENCFF698HBZ.bam -O rep1.bam
    wget --no-check-certificate https://www.encodeproject.org/files/ENCFF778JLA/@@download/ENCFF778JLA.bam -O rep2.bam

    # Control bam
    wget --no-check-certificate https://www.encodeproject.org/files/ENCFF948ROL/@@download/ENCFF948ROL.bam -O control.bam

# Download reference files and blacklist regions

    # download genome refrence
    wget --no-check-certificate https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz \
    -O hg38.genome.fa.gz
    gunzip hg38.genome.fa.gz

    # index genome reference
    samtools faidx hg38.genome.fa

    # download chrom sizes
    wget --no-check-certificate https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv

    # exclude alt contigs and chrEBV
    grep -v -e '_' -e 'chrEBV' GRCh38_EBV.chrom.sizes.tsv > hg38.chrom.sizes
    rm GRCh38_EBV.chrom.sizes.tsv

    # make file with chromosomes only
    awk '{print $1}' hg38.chrom.sizes > chroms.txt

    # download blacklist
    wget --no-check-certificate https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz -O blacklist.bed.gz
    gunzip blacklist.bed.gz

# Preprocessing steps to generate bigwig counts tracks
## Merge the two replicates and create and index

    # Experiment
    samtools merge -f merged.bam rep1.bam rep2.bam
    samtools index merged.bam

    # Control
    samtools index control.bam

## Create bigwig files using bedtools via intermediate bedGraph files

    # Experiment
    # get coverage of 5’ positions of the plus strand
    samtools view -b merged.bam $(cut -f 1 hg38.chrom.sizes) | \
    	bedtools genomecov -5 -bg -strand + -ibam stdin | \
    	sort -k1,1 -k2,2n > plus.bedGraph

    # get coverage of 5’ positions of the minus strand
    samtools view -b merged.bam $(cut -f 1 hg38.chrom.sizes) | \
            bedtools genomecov -5 -bg -strand - -ibam stdin | \
            sort -k1,1 -k2,2n > minus.bedGraph

    # Convert bedGraph files to bigWig files
    bedGraphToBigWig plus.bedGraph hg38.chrom.sizes plus.bw
    bedGraphToBigWig minus.bedGraph hg38.chrom.sizes minus.bw


    # Control
    # get coverage of 5’ positions of the control plus strand
    samtools view -b control.bam $(cut -f 1 hg38.chrom.sizes) | \
            bedtools genomecov -5 -bg -strand + -ibam stdin | \
            sort -k1,1 -k2,2n > control_plus.bedGraph

    # get coverage of 5' positions of the control minus strand
    samtools view -b control.bam $(cut -f 1 hg38.chrom.sizes) | \
            bedtools genomecov -5 -bg -strand - -ibam stdin | \
            sort -k1,1 -k2,2n > control_minus.bedGraph

    # Convert bedGraph files to bigWig files
    bedGraphToBigWig control_plus.bedGraph hg38.chrom.sizes control_plus.bw
    bedGraphToBigWig control_minus.bedGraph hg38.chrom.sizes control_minus.bw

## Download identify peaks

    wget https://www.encodeproject.org/files/ENCFF869DAZ/@@download/ENCFF869DAZ.bed.gz -O peaks.bed.gz
    gunzip peaks.bed.gz

## Organize data

    mkdir GR
    mkdir GR/data
    mv *.bw GR/data
    mv peaks.bed GR/data

    mkdir GR/reference
    mv hg38.genome.fa* GR/reference
    mv hg38.chrom.sizes GR/reference
    mv chroms.txt GR/reference
    mv blacklist.bed GR/reference

# Outlier removal
    {
        "0": {
            "signal": {
                "source": ["GR/data/plus.bw",
                           "GR/data/minus.bw"]
            },
            "loci": {
                "source": ["GR/data/peaks.bed"]
            },
            "bias": {
                "source": ["GR/data/control_plus.bw",
                           "GR/data/control_minus.bw"],
                "smoothing": [null, null]
            }
        }
    }
make json file like this as input_outliers.json

    bpnet-outliers \
        --input-data input_outliers.json  \
        --quantile 0.99 \
        --quantile-value-scale-factor 1.2 \
        --task 0 \
        --chrom-sizes GR/reference/hg38.chrom.sizes \
        --chroms $(paste -s -d ' ' GR/reference/chroms.txt) \
        --sequence-len 1000 \
        --blacklist GR/reference/blacklist.bed \
        --global-sample-weight 1.0 \
        --output-bed GR/data/peaks_inliers.bed

and then run this code

# Citation

This project is based on [BPNet-refactor](https://github.com/kundajelab/bpnet-refactor.git) by the Kundaje Lab.

If you use this repository, please cite the original BPNet paper:

> Avsec, Ž., Weilert, M., Shrikumar, A. et al.  
> *Base-resolution models of transcription-factor binding reveal soft motif syntax.*  
> Nature Genetics 53, 354–366 (2021). https://doi.org/10.1038/s41588-021-00782-6

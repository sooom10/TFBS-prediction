# TFBS-prediction

# 환경설정

    conda create --name bpnet python=3.7
    conda activate bpnet
    pip install git+https://github.com/kundajelab/bpnet-refactor.git

# 필요 라이브러리 설치(여기서 파이썬 버전 호환 건으로 만약 설치되지 않는 library가 있다면 그 라이브러리만 새로운 환경 만들어서 설치)

    conda install -y -c bioconda samtools=1.1 bedtools ucsc-bedgraphtobigwig

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

# Background_generation
## gc matched negatives

      # bpnet-gc-reference - get gc content after binning the entire genome into bins - You might be choose to run just once for a genome for a specific input sequence length and reuse the genomewide_gc_stride_flank_size.gc.bed output for other datasets

    bpnet-gc-reference \
            --ref_fasta GR/reference/hg38.genome.fa \
            --chrom_sizes GR/reference/hg38.chrom.sizes \
            --output_prefix GR/reference/genomewide_gc_stride_1000_flank_size_1057.gc \
            --inputlen 2114 \
            --stride 1000
        
    
    bpnet-gc-background \
            --ref_fasta GR/reference/hg38.genome.fa \
            --peaks_bed GR/data/peaks_inliers.bed \
            --out_dir GR/data/ \
            --ref_gc_bed GR/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed \
            --output_prefix gc_negatives \
            --flank_size 1057 \
            --neg_to_pos_ratio_train 4

# Train the model
make these two json file

## input_data.json

    {
        "0": {
            "signal": {
                "source": ["GR/data/plus.bw", 
                           "GR/data/minus.bw"]
            },
            "loci": {
                "source": ["GR/data/peaks_inliers.bed"]
            },
            "background_loci": {
                "source": ["GR/data/gc_negatives.bed"],
                "ratio": [0.33]
            },
            "bias": {
                "source": ["GR/data/control_plus.bw",
                           "GR/data/control_minus.bw"],
                "smoothing": [null, null]
            }
        }
    }

## best_params.json

    {
        "input_len": 2114,
        "output_profile_len": 1000,
        "weight_decay": 0.0001,
        "motif_module_params": {
            "filters": [
                128
            ],
            "kernel_sizes": [
                21
            ],
            "padding": "valid"
        },
    
        "transformer_params": {
            "head_size": 32,
            "num_heads": 4,
            "ff_dim": 256,
            "num_blocks": 3,
            "dropout": 0.15
        },

        "profile_head_params": {
            "filters": -1,  
            "kernel_size":  75,
            "padding": "valid"
        },
        "counts_head_params": {
            "units": [
                -1
            ],
            "dropouts": [0.0],
            "activations": ["linear"]
        },
        "profile_bias_module_params": {
            "kernel_sizes": [
                1
            ]
        },
        "counts_bias_module_params": {},
        "use_attribution_prior": false,
        "attribution_prior_params": {
            "frequency_limit": 150,
            "limit_softness": 0.2,
            "grad_smooth_sigma": 3,
            "profile_grad_loss_weight": 200,
            "counts_grad_loss_weight": 100        
        },
        "loss_weights": [
            1,
            42
        ],
        "counts_loss": "MSE"
    }

## Run this code to find optimal loss weight (위에 최적 loss_weight 찾기 위해(나온 결과로 42로 적혀있는 부분 대체하면 됨))

    bpnet-counts-loss-weight --input-data input_data.json

# Splits datasets
## splits.json

    {
        "0": {
            "test":
                ["chr7", "chr13", "chr17", "chr19", "chr21", "chrX"],
            "val":
                ["chr10", "chr18"],
            "train":
                ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr8", "chr9", "chr11", "chr12", "chr14", "chr15", "chr16", "chr21", "chr22", "chrY"]
        }
    }

## Train command

    BASE_DIR=GR
    DATA_DIR=$BASE_DIR/data
    MODEL_DIR=$BASE_DIR/models
    REFERENCE_DIR=reference
    CHROM_SIZES=$BASE_DIR/$REFERENCE_DIR/hg38.chrom.sizes
    REFERENCE_GENOME=$BASE_DIR/$REFERENCE_DIR/hg38.genome.fa
    CV_SPLITS=$BASE_DIR/splits.json
    INPUT_DATA=$BASE_DIR/input_data.json
    MODEL_PARAMS=$BASE_DIR/bpnet_params.json
    
    # 파일 생성 안한 경우만
    # mkdir $MODEL_DIR


    bpnet-train \
            --input-data $INPUT_DATA \
            --output-dir $MODEL_DIR \
            --reference-genome $REFERENCE_GENOME \
            --chroms $(paste -s -d ' ' $BASE_DIR/$REFERENCE_DIR/chroms.txt) \
            --chrom-sizes $CHROM_SIZES \
            --splits $CV_SPLITS \
            --model-arch-name BPNet \
            --model-arch-params-json $MODEL_PARAMS \
            --sequence-generator-name BPNet \
            --model-output-filename GR_model \
            --input-seq-len 2114 \
            --output-len 1000 \
            --shuffle \
            --threads 10 \
            --epochs 100 \
    	--batch-size 32 \
    	--reverse-complement-augmentation \
    	--early-stopping-patience 10 \
    	--reduce-lr-on-plateau-patience 5 \
            --learning-rate 0.001

## Prediction command

    PREDICTIONS_DIR=$BASE_DIR/predictions_and_metrics2

    # 이것도 최초 한 번만
    # mkdir $PREDICTIONS_DIR

    bpnet-predict \
            --model $MODEL_DIR/GR_model_split000 \
            --chrom-sizes $CHROM_SIZES \
            --chroms chr7 chr13 chr17 chr19 chr21 chrX \
            --test-indices-file None \
            --reference-genome $REFERENCE_GENOME \
            --output-dir $PREDICTIONS_DIR \
            --input-data $INPUT_DATA \
            --sequence-generator-name BPNet \
            --input-seq-len 2114 \
            --output-len 1000 \
            --output-window-size 1000 \
            --batch-size 32 \
            --reverse-complement-average \
            --threads 2 \
            --generate-predicted-profile-bigWigs

## Compute importance score

    SHAP_DIR=$BASE_DIR/shap
    
    # 이것도 최초 한 번만
    # mkdir $SHAP_DIR
    
    bpnet-shap \
            --reference-genome $REFERENCE_GENOME \
            --model $MODEL_DIR/GR_model_split000  \
            --bed-file $DATA_DIR/peaks_inliers.bed \
            --output-dir $SHAP_DIR \
            --input-seq-len 2114 \
            --control-len 1000 \
            --task-id 0 \
            --input-data $INPUT_DATA \
            --chrom-sizes $CHROM_SIZES \
            --generate-shap-bigWigs

# Run tf-modisco to find motif

    conda create --name tfmodisco python=3.10
    conda activate tfmodisco
    pip install modisco-lite
    
    MODISCO_DIR=$BASE_DIR/modisco
    MODISCO_COUNTS_DIR=$BASE_DIR/modisco/counts
    MODISCO_PROFILE_DIR=$BASE_DIR/modisco/profile
    
    mkdir -p $MODISCO_DIR
    mkdir -p $MODISCO_COUNTS_DIR
    mkdir -p $MODISCO_PROFILE_DIR
    
    modisco motifs\
        --max_seqlets 50000 \
        --h5py $SHAP_DIR/profile_scores.h5 \
        -o $MODISCO_PROFILE_DIR/profile_modisco_scores.h5 \
        --trim_size 20 \
        --initial_flank_to_add 5 \
        --final_flank_to_add 10
    
    modisco motifs\
        --max_seqlets 50000 \
        --h5py $SHAP_DIR/counts_scores.h5 \
        -o $MODISCO_COUNTS_DIR/counts_modisco_scores.h5 \
        --trim_size 20 \
        --initial_flank_to_add 5 \
        --final_flank_to_add 10

# Reference

This project is based on [BPNet-refactor](https://github.com/kundajelab/bpnet-refactor.git) by the Kundaje Lab.

If you use this repository, please cite the original BPNet paper:

> Avsec, Ž., Weilert, M., Shrikumar, A. et al.  
> *Base-resolution models of transcription-factor binding reveal soft motif syntax.*  
> Nature Genetics 53, 354–366 (2021). https://doi.org/10.1038/s41588-021-00782-6

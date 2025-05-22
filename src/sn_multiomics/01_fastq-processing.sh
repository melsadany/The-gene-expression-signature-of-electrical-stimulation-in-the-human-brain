#! /bin/bash

#### set-up
P_DIR=/Dedicated/jmichaelson-wdata/msmuhammad/projects/neurosurgery

OUT_DIR=${P_DIR}/data/raw/ME-processing/01
mkdir -p ${OUT_DIR}

cd ${OUT_DIR}


##### tools/refs
# https://timoast.github.io/sinto/basic_usage.html
CA=/Dedicated/jmichaelson-wdata/msmuhammad/workbench/cellranger-arc-2.0.2/bin/cellranger-arc
REF=/Dedicated/jmichaelson-wdata/msmuhammad/workbench/cellranger-arc-2.0.2/refs/refdata-cellranger-arc-GRCh38-2020-A-2.0.0
REF_mm=/Dedicated/jmichaelson-wdata/msmuhammad/workbench/cellranger-arc-2.0.2/refs/refdata-cellranger-arc-mm10-2020-A-2.0.0


#### PROCESS ALL SAMPLES
samples=("E1A" "E1B" "E2A" "E2B" "E3A" "E3B" "E4A" "E4B" "G1A" "G1B" "G2A" "G2B")
for sample in "${samples[@]}"; do
    echo $sample
    
    ${CA} count \
        --id=$sample \
        --reference=${REF} \
        --libraries=${P_DIR}/data/${sample}_sample-libraries.csv \
        --localcores=100 \
        --localmem=512
done


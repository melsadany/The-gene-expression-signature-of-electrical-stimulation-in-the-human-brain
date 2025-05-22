#$ bash


#  create the environment for DCA first
# python3.8 -m venv dca-env
# source ~/dca-env/bin/activate
# pip install --upgrade pip
# pip install tensorflow==2.4.0 keras==2.4.3 dca==0.3.3


source ~/dca-env/bin/activate

P_DIR=/Dedicated/jmichaelson-wdata/msmuhammad/projects/neurosurgery
cd $P_DIR
OUT=data/derivatives/rna-imputation/dca/out

## run for stimulated ones
mkdir -p ${OUT}_stimulated
dca --threads 35 \
    data/derivatives/rna-imputation/dca/dca-input-integrated-stimulated-rna.csv \
    ${OUT}_stimulated
deactivate

## run for unstimulated ones
source ~/dca-env/bin/activate
mkdir -p ${OUT}_unstimulated
dca --threads 35 \
    data/derivatives/rna-imputation/dca/dca-input-integrated-unstimulated-rna.csv \
    ${OUT}_unstimulated

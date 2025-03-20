![alt text](https://github.com/interactivereport/SpliceHarmonization/blob/main/figures/SpliceHarmonization%20LOGO.svg)

## About
SpliceHarmonization is an integrated approach that leverages the strengths of different alternative splicing detection methods, thereby facilitating robust and reliable analysis of splicing events with annotations.
### Required dependencies
- rMATS: https://github.com/Xinglab/rmats-turbo/blob/master/README.md
- LeafCutter: https://github.com/davidaknowles/leafcutter/tree/master/conda_recipe
- MAJIQ: https://biociphers.bitbucket.io/majiq-docs-academic/getting-started-guide/installing.html
- StringTie: https://ccb.jhu.edu/software/stringtie
- SAMTools
- Python
- R
  
### Installation 
- pull SpliceHarmonization from github \
        `git clone https://github.com/interactivereport/SpliceHarmonization.git`
- enter the folder and install SpliceHarmonization \
        `cd install` \
        `conda create -n env -f install.yml`

### SpliceHarmonization
- prepare the config.yml and input data \
- start SpliceHarmonization \
  `python ./spliceharmonization.py -cfig ./config.yml`
  
### Splicing events simulation 
- enter the folder and intall env \
        `cd simulation/install` \
        `conda env create -f simulation.yml`
- prepare the config.yml and input data \
        Details in simulation_data.tar.gz at https://zenodo.org/uploads/15032402 
- start simulation \
        `python ./SpliceSimulator.py -cfig ./config.yml`

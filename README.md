![alt text](https://github.com/interactivereport/SpliceHarmonization/blob/main/figures/SpliceHarmonization%20LOGO.svg)

## About
SpliceHarmonization is an integrated approach that leverages the strengths of different alternative splicing detection methods, thereby facilitating robust and reliable analysis of splicing events with annotations.

### Installation 
- pull SpliceHarmonization from github \
        `git clone https://github.com/interactivereport/SpliceHarmonization.git`
- enter the folder and install SpliceHarmonization \
        `cd install` \
        `conda create -n env -f install.yml`

### SpliceHarmonization
- enter the folder and install env \
      `cd simulation/install` \
      `conda env create -f simulation.yml`
- prepare the config.yml and input data \
- start SpliceHarmonization
  `python ./spliceharmonization.py -cfig ./config.yml`
  
### Splicing events simulation 
- enter the folder and intall env \
        `cd simulation/install` \
        `conda env create -f simulation.yml`
- prepare the config.yml and input data \
        Details in simulation_data.tar.gz at https://zenodo.org/uploads/15032402 
- start simulation \
        `python ./SpliceSimulator.py -cfig ./config.yml`

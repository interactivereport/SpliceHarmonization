![alt text](https://github.com/interactivereport/SpliceHarmonization/blob/main/figures/SpliceHarmonization%20LOGO.svg)

## About
SpliceHarmonization is an integrated approach that leverages the strengths of different alternative splicing detection methods, thereby facilitating robust and reliable analysis of splicing events with annotations.
### Required dependencies
SpliceHarmonization requires execution on a high-performance computing (HPC) platform.
- rMATS: https://github.com/Xinglab/rmats-turbo/blob/master/README.md
- LeafCutter: https://github.com/davidaknowles/leafcutter/tree/master/conda_recipe
- MAJIQ: https://biociphers.bitbucket.io/majiq-docs-academic/getting-started-guide/installing.html
- StringTie: https://ccb.jhu.edu/software/stringtie
- SAMTools
- Anaconda
- Python
- R
  
### Installation 
- pull SpliceHarmonization from github \
        `git clone https://github.com/interactivereport/SpliceHarmonization.git`
- enter the folder and install SpliceHarmonization \
        `cd install` \
        `conda create -n env -f install.yml`

### SpliceHarmonization
- prepare the config.yml and input data 
    - config.yml: please see an exmample in `test/config.yml`
    - samplesheet.csv: please see an example in `test/samplesheet.csv`
    - comparison.csv: please see an example in `test/comparison.csv`
- start SpliceHarmonization \
  `python ./spliceharmonization.py -cfig ./config.yml`

- Output tree view
  ```bash
{output_path} (#output directory)
├── {timestamp}_{userName} (# ouptput folder)
│   ├── {Comparison}
│   │   ├── harm
│   │   ├── junction_prep
│   │   ├── out
│   │   │   ├── all_gene
│   │   │   │   ├── {gene_names}.csv (# events output for each gene)
│   │   │   ├── all_gene_junction (# empty folder if junction_filter = False)
│   │   │   └── mxe
│   │   ├── run1.e.log
│   │   ├── run1.o.log
│   │   ├── run1.sh
│   │   ├── run2.e.log
│   │   ├── run2.o.log
│   │   ├── run2.sh
│   │   ├── run3.e.log
│   │   ├── run3.o.log
│   │   ├── run3.sh
│   │   ├── run4.e.log
│   │   ├── run4.o.log
│   │   └── run4.sh
│   ├── _cmdline_20250321.103736
│   ├── comparison.csv
│   ├── config.yml
│   ├── leafcutter (# orginial leafcutter output)
│   ├── majiq (# orginial majiq output)
│   ├── rmats (# orginial rmats output)
│   ├── samplesheet.csv
│   ├── status.log
│   └── stringtie (# orginial stringtie output)
├── BAM_index_sbatch # tmp folder, can delete later
├── leafcutter_output # tmp folder, can delete later
├── majiq_output # tmp folder, can delete later
├── rmats_output # tmp folder, can delete later
└── stringtie_output # tmp folder, can delete later
```
  
### Splicing events simulation 
- enter the folder and intall env \
        `cd simulation/install` \
        `conda env create -f simulation.yml`
- prepare the config.yml and input data \
        Details in simulation_data.tar.gz at https://zenodo.org/uploads/15032402 
- start simulation \
        `python ./SpliceSimulator.py -cfig ./config.yml`

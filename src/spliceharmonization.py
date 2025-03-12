import sys
import SplicePrep
import SpliceEvent
import os,sys
import shutil
import argparse

def get_arg(base_dir):
    strPipePath =  os.path.join(base_dir, "SpliceEvent")
    parser = argparse.ArgumentParser(description="How to use SpliceEvent pipeline")
    parser.add_argument('-cfig', require=True, dest = 'config_file', help ='file directionary')
    # parser.add_argument('-indir',required=True, dest="indir", help ='data directory')
    parser.add_argument('-wdir',required=True, dest="wdir", help='working directory')
    parser.add_argument('-ref', require=True, dest='ref', help="reference .gtf",default=None)
    parser.add_argument('-bamindir', dest='bamindir', help='bam data directory', default=None)
    parser.add_argument('-specific_comparison', dest= 'specific_comparison', help = "identify your specific comparison substring separated by ,", default=None)
    parser.add_argument('-junction_filter', dest='junction_filter', type=bool, help="open or close this filter", default=False)
    parser.add_argument('-majiq_cutoff_val', dest='majiq_cutoff_val', type=float, default=0.95)
    parser.add_argument('-junctionFC', dest='junctionFC', type=float, default=1.2)
    parser.add_argument('-junctionMAX', dest='junctionMAX', type=float, default=25)
    parser.add_argument('-samplesheet', dest='samplesheet', help= 'sample sheet file', default=None)
    parser.add_argument('-addexons', dest='addexons', help= 'added exons sheet file', default=os.path.join(strPipePath,'data','example_add_exon.csv'))
    parser.add_argument('-cmdonly',dest='cmdonly',help='only save all slurm bash files, not submit',default=False)
    args = parser.parse_args()
    return args

def powerMsg():
    print("\nSpliceHarmonization Pipeline Powered by the Biogen Computational Analytical Group "
          "[yirui.chen@biogen.com]\n")
    
def main():
    args = get_arg(base_dir)
    base_dir = os.path.dirname(os.path.realpath(__file__))
    src_dir = os.path.join(base_dir, "src")
    config_src = os.path.join(src_dir, "SplicePrep/config.yml")

    if len(sys.argv) < 2:
        shutil.copy(config_src, os.getcwd())
        print("An empty config file copied to the current working directory")
        sys.exit(0)

    sysarg = sys.argv[1]
    if os.path.isdir(sysarg):
        shutil.copy(config_src, sysarg)
        print("An empty config file copied to the specified folder: {}".format(sysarg))
        sys.exit(0)

    print("Starting SplicePrep")
    spliceprep_path = SplicePrep.main(args.config_file)
    print("The project folder is located at:", spliceprep_path)
    args.indir = spliceprep_path

    print("Starting SpliceEvent")
    SpliceEvent.main(args)

    powerMsg()

if __name__ == '__main__':
    main()
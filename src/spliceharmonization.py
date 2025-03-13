import sys
import SplicePrep
import SpliceEvent
import os,sys
import shutil
import argparse
from datetime import datetime
import yaml
from SpliceEventSumbit import SpliceHarmJobSubmitter

def get_arg(base_dir):
    strPipePath =  os.path.join(base_dir, "SpliceEvent")
    parser = argparse.ArgumentParser(description="How to use SpliceEvent pipeline")
    parser.add_argument('-cfig', required=True, dest = 'cfig', help ='file directionary')
    parser.add_argument('-wdir',required=True, dest="wdir", help='working directory')
    # parser.add_argument('-ref', required=True, dest='ref', help="reference .gtf",default=None)
    args = parser.parse_args()
    return args

def powerMsg():
    print("\nSpliceHarmonization Pipeline Powered by the Biogen Computational Analytical Group "
          "[yirui.chen@biogen.com]\n")
    
def main():
    
    print("Script started", flush=True)
    timestamp = datetime.now().strftime('%Y%m%d%H%M%S')
    user = os.getenv('USER')
    base_dir = os.path.dirname(os.path.realpath(__file__))
    args = get_arg(base_dir)

    

    strPrj = os.path.join(args.wdir, f"{timestamp}_{user}")
    # strPrj = os.path.join(args.wdir, '20250312233909_ychen12')
    
    with open(args.cfig, 'r') as file:
        g_config = yaml.safe_load(file)
    

    g_config['wdir'] = strPrj
    g_config['indir'] = strPrj

    with open('config.yml', 'w') as file:
        yaml.dump(g_config, file)

    print("Arguments parsed:", args, flush=True)
    
    src_dir = os.path.join(base_dir, "src")
    config_src = os.path.join(src_dir, "SplicePrep/config.yml")
    
    if len(sys.argv) < 2:
        shutil.copy(config_src, os.getcwd())
        print("An empty config file copied to the current working directory", flush=True)
        sys.exit(0)
    
    sysarg = sys.argv[1]
    if os.path.isdir(sysarg):
        shutil.copy(config_src, sysarg)
        print("An empty config file copied to the specified folder: {}".format(sysarg), flush=True)
        sys.exit(0)
    
    print("Before calling SplicePrep.main()", flush=True)
    SplicePrep.main(args.cfig, args.wdir)

    # print(g_config)

    g_config = yaml.safe_load(args.cfig)

    # print("After SplicePrep.main(), project folder is located at:", spliceprep_path, flush=True)
    # delattr(args, 'cfig')
    print("Before calling SpliceEvent.main()", flush=True)

    SpliceHarmJobSubmitter.main(args.cfig)
    
    print(args)
    print("After calling SpliceEvent.main()", flush=True)
    
    powerMsg()

if __name__ == '__main__':
    main()
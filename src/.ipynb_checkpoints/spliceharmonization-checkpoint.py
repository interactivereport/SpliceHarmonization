import sys
import SplicePrep
import SpliceEvent
import os,sys
import shutil
import argparse
from datetime import datetime
import yaml

def get_arg():
    parser = argparse.ArgumentParser(description="How to use SpliceEvent pipeline")
    parser.add_argument('-cfig', required=True, dest = 'cfig', help ='file directionary')
    # parser.add_argument('-wdir',required=True, dest="wdir", help='working directory')
    args = parser.parse_args()
    return args

def powerMsg():
    print("\nSpliceHarmonization Pipeline Powered by the Biogen Computational Analytical Group "
          "[yirui.chen@biogen.com]\n")
    
def yaml_dump(args, config):
    prefix = args.cfig.split('.')[0]
    insert = True
    with open(f'{prefix}.tmp.yaml', 'w') as f_out:
        with open(f'{args.cfig}', 'r') as f_in:
            for line in f_in.readlines():
                if line.startswith('#'):
                    f_out.write(line)
                    if insert:
                        f_out.write(f"wdir: {config['wdir']}\n")
                        insert = False
                else:
                    f_out.write(line)
    shutil.move(f'{prefix}.tmp.yaml', args.cfig)

def main():
    print("Script started", flush=True)
    args = get_arg()
    # load the config file
    with open(args.cfig, 'r') as file:
        g_config = yaml.safe_load(file)

    assert 'output_path' in g_config, "output_path is not defined in the config file"

    if 'wdir' not in g_config:
        # create a subfolder using username and timestamp if there is no wdir specified.
        timestamp = datetime.now().strftime('%Y%m%d%H%M%S')
        user = os.getenv('USER')
        strPrj = os.path.join(g_config['output_path'], f"{timestamp}_{user}")
        os.makedirs(strPrj,exist_ok=True)
        g_config['wdir'] = strPrj
        # write the updated config file
        print('Rewrite config file')
        yaml_dump(args, g_config)

    print("Arguments parsed:", args, flush=True)
  

    print("Before calling SplicePrep.main()", flush=True)
    SplicePrep.main(g_config)

    print("Before calling SpliceEvent.main()", flush=True)
    SpliceEvent.main(g_config)
    print("After calling SpliceEvent.main()", flush=True)
    
    powerMsg()

if __name__ == '__main__':
    main()
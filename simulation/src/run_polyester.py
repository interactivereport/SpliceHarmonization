import subprocess
import argparse
def main(args, num_reads, cntl_alt_ratio):
    

    # program='/home/ychen12/splicing_test/edge_splicing/splicing_harmonization/DIRS/simulators'
    # outdir="/mnt/depts/dept04/compbio/projects/TST12320/simulation_demo/fasta"
    # simulated_odir="/mnt/depts/dept04/compbio/projects/TST12320/simulation_demo/simulation"
    r_script_path = f"{args.progpath}/polyester.R"
    ref_path = f"{args.fastadir}/combined_ref.fasta"
    alt_path = f"{args.fastadir}/combined_alt.fasta"

    for num_read in num_reads:
        for  cntl_p, alt_p in cntl_alt_ratio:
            print(cntl_p)
            print(alt_p)
            command = [
                "Rscript", r_script_path,
                "--ref_input", ref_path,
                "--alt_input", alt_path,
                "--output_dir", args.outdir,
                "--num_reads", str(num_read), 
                "--cntl_propo", str(cntl_p),
                "--alt_propo", str(alt_p)
            ]

            result = subprocess.run(command, capture_output=True, text=True)
            # print(command)


            if result.returncode == 0:
                print("R script executed successfully.")
                print("Output:", result.stdout)
            else:
                print("Error executing R script.")
                print("Error output:", result.stderr)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="run polyester.")
    parser.add_argument('-progpath', dest="progpath", help ='prog path', required=True)
    parser.add_argument('-fastadir', dest='fastadir', help="fasta file path", required=True)
    parser.add_argument('-outdir', dest='outdir', help='output files location', required=True)
    parser.add_argument('-lowCounts', dest='lowCounts', type=int, help='lowCounts', required=True)
    parser.add_argument('-highCounts', dest='highCounts', type=int, help='highCounts', required=True)
    args = parser.parse_args()
    # num_reads = [20000, 500000]
    num_reads = [args.lowCounts, args.highCounts]
    cntl_alt_ratio = cntl_alt_ratio = [(0,1), (0.2, 0.8), (0.5, 0.5), (0.8, 0.2)]
    main(args, num_reads, cntl_alt_ratio)

    
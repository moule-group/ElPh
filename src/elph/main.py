from elph.workflow import run_j0, run_lambda, run_disp_j, run_matrix, run_tlt_mobility, run_svd_projection, submit_slurm_script
import elph.utils as ut
import argparse
import sys

def validate_args(args):
    if args.nmol < 3 or args.nmol > 5:
        raise ValueError(f"{args.nmol} is invalid 'nmol' count. 'nmol' must be between 3 and 5")
    
    if args.workflow < 0 or args.workflow > 3:
        raise ValueError(f"{args.workflow} is not a valid workflow.")

def main():
    """ Main function to run el-ph coupling calculation
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-w", "--workflow", type=int, required=True, help="Workflow number to run.")
    parser.add_argument("--overwrite", action='store_true', default=False, help="Overwrite existing files and directories.")
    parser.add_argument("-q", "--mesh", type=int, default=[8,8,8], nargs=3, help="Defining a mesh grid. Defaults to [8,8,8].")
    parser.add_argument("-b", "--basis", type=str, default=['6-31G**','6-311G**'], nargs=2, help="Gaussian basis sets, first:local; second:non-local")
    parser.add_argument("-f", "--functional", type=str, default=['b3lyp','b3lyp'], nargs=2, help="Gaussian functional, first:local; second:non-local")
    parser.add_argument("-n", "--nmol", type=int, default=3, help="The number of molecules to extract.")
    parser.add_argument("-s", "--supercell", type=int, nargs=3, default=(2,2,2), help="Supercell matrix. Defaults to (2,2,2).")
    parser.add_argument("-o", "--output", type=str, default="tlt_mobility", help="Mobility calculation output name")
    # parser.add_argument("-homo", "--homo", type=str, default=True, help="P-type: HOMO; N-type: LUMO")
    parser.add_argument("-l", "--local", action='store_true', default=False, help="Run ElPh on local machine instead of submitting a slurm script.")
    parser.add_argument("-a", "--account", type=str, help="Account name for slurm job.")
    parser.add_argument("-time", "--time", type=str, default="02:00:00", help="Time limit for slurm job (hh:mm:ss). Defaults to 1:00:00.")
    parser.add_argument("-g", "--gpu", action="store_true", default=True, help="Run slurm job on GPU.")
    parser.add_argument("-H", "--hpc", type=int, nargs=3, default=[8,32,8], help="Slurm job flag values for srun: -n, -c, and -G, respectively.")

    args = parser.parse_args()
    try:
        validate_args(args)
    except Exception as e:
        ut.throw_error(e)
    
    ut.print_start()

    try:
        # submit script
        if not args.local:
            submit_slurm_script(args)
            return
            
        # run locally 
        if args.workflow == 1: 
            run_j0(args.basis, args.functional, args.supercell, args.nmol, args.overwrite) # Run Gaussian with optimization 

        elif args.workflow == 2: # Run workflow 2 (need to finish workflow 1 first)
            run_lambda(args.basis, args.functional)
            run_disp_j(args.basis, args.functional, args.nmol) # Create displaced dimers and calculate J_ij of dimers.
            run_matrix(args.mesh) # Calculate electron phonon coupling matrix (including local and non-local part)
        
        elif args.workflow == 3:  # Calculate the mobility
            run_tlt_mobility(output=args.output)                

        ut.print_end()

    except KeyboardInterrupt:
        ut.throw_error("Interrupted by user!")

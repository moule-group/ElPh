from elph.workflow import run_j0, run_lambda, run_disp_j, run_matrix, run_tlt_mobility, run_svd_projection
import elph.utils as ut
import argparse
import sys

def main():
    """ Main function to run el-ph coupling calculation
    """
    parser = argparse.ArgumentParser() # Create the parser
    parser.add_argument("-q", "--mesh", type=int, default=[8,8,8], nargs=3, help='Defining a mesh grid. (Defaults to [8,8,8]') # Add an argument: mesh
    parser.add_argument("-b", "--basis", type=str, default=['6-311G*','6-311G**'], nargs=2, help='Gaussian basis sets, first:local; second:non-local') # Add an argument: basis
    parser.add_argument("-f", "--functional", type=str, default=['b3lyp','b3lyp'], nargs=2, help='Gaussian functional, first:local; second:non-local') # Add an argument: functional
    parser.add_argument("-n", "--nmol", type=int, default=3, help='The number of molecules will be extracted') # Add an argument: nmol
    parser.add_argument("-s", "--supercell", type=int, nargs=3, default=(2,2,2), help='The supercell matrix (Defaults to (2,2,2) )') # Add an argument: supercell
    parser.add_argument("-o", "--output", type=str, default=None, help='Mobility calculation output name') # Add an argument: filename
    parser.add_argument("-svd", "--svdqpts", type=int, default=1, help='Number of qpoints that SVD projection will apply') # Add an argument: svd
    parser.add_argument("-w", "--workflow", type=int, default=1, help='Type in the workflow number to run corresponding simulation') # Add an argument: workflow
    parser.add_argument("-homo", "--homo", type=str, default=True, help='P-type: HOMO; N-type: LUMO') # Add an argument: homo
    args = parser.parse_args() # Parse the argument
    
    ut.print_start()

    try:
        if args.workflow == 1: 
            run_j0(args.basis, args.functional, args.supercell, args.nmol) # Run Gaussian with optimization 
            ut.print_end()

        if args.workflow == 2: # Run the workflow 2 (need to finish workflow 1 first)
            run_lambda(args.basis, args.functional)
            run_disp_j(args.basis, args.functional, args.nmol) # Create displaced dimers and calculate J_ij of dimers.
            run_matrix(args.mesh) # Calculate electron phonon coupling matrix (including local and non-local part)
            ut.print_end()
        
        if args.workflow == 3:  # Calculate the mobility
            if args.output is None:
                run_tlt_mobility() # Calculate the mobility
                ut.print_end()
            else:
                run_tlt_mobility(output=args.output)
                ut.print_end()

    except KeyboardInterrupt:
        ut.print_error("Interrupted by user!")
        sys.exit(0)
from elph.workflow import getGeometry, run_j0, run_disp_j, run_matrix, run_tlt_mobility
import elph.utils as ut
import argparse
import sys

def main():
    """ Main function to run el-ph coupling calculation
    """
    parser = argparse.ArgumentParser() # Create the parser
    parser.add_argument("-q", "--mesh", type=int, default=[8,8,8], nargs=3, help='Defining a mesh grid. (Defaults to [8,8,8]') # Add an argument: mesh
    parser.add_argument("-b", "--basis", type=str, default='3-21G*', help='Gaussian basis sets') # Add an argument: basis
    parser.add_argument("-m", "--mol", type=int, nargs=3, help='The numbering of molecule 1 2 and 3') # Add an argument: mol
    parser.add_argument("-s", "--supercell", type=int, nargs=3, default=[2,2,2], help='The supercell matrix (Defaults to [2,2,2])') # Add an argument: supercell
    parser.add_argument("-mu", "--mobility", action='store_true', help='Calculate the mobility') # Add an argument: mobility
    parser.add_argument("-o", "--output", type=str, default=None, help='Mobility calculation output name') # Add an argument: filename
    args = parser.parse_args() # Parse the argument
    
    ut.print_start()

    try:
        if not args.mobility:
            run_j0(args.mol, args.basis) # Run Gaussian with optimization 
            run_disp_j(args.basis) # Create displaced dimers and calculate J_ij of dimers.
            run_matrix(args.mesh,args.supercell) # Calculate electron phonon coupling matrix
            ut.print_end()
        
        else:  # Calculate the mobility
            if args.output is None:
                run_tlt_mobility() # Calculate the mobility
                ut.print_end()
            else:
                run_tlt_mobility(output=args.output)
                ut.print_end()

    except KeyboardInterrupt:
        ut.print_error("Interrupted by user!")
        sys.exit(0)
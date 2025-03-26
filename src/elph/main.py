from elph.workflow import getGeometry, run_j0, run_disp_E, run_disp_j, run_matrix, run_tlt_mobility, run_svd_projection
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
    parser.add_argument("-o", "--output", type=str, default=None, help='Mobility calculation output name') # Add an argument: filename
    parser.add_argument("-svd", "--svdqpts", type=int, default=1, help='Number of qpoints that SVD projection will apply') # Add an argument: svd
    parser.add_argument("-w", "--workflow", type=int, default=1, help='Type in the workflow number to run corresponding simulation') # Add an argument: workflow
    args = parser.parse_args() # Parse the argument
    
    ut.print_start()

    try:
        if args.workflow == 1: 
            run_j0(args.mol, args.basis) # Run Gaussian with optimization 
            run_disp_j(args.basis) # Create displaced dimers and calculate J_ij of dimers.
            run_disp_E() # Collect energy for displaced molecules in the materials
            run_matrix(args.mesh,args.supercell) # Calculate electron phonon coupling matrix (including local and non-local part)
            ut.print_end()

        if args.workflow == 2: # Run the workflow 2 (need to finish workflow 1 first)
            run_svd_projection(args.svdqpts) # Run the SVD projection
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
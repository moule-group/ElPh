from elph.workflow import getGeometry, run_j0, run_disp_j, run_matrix
import elph.utils as ut
import sys

def main():
    """ Main function to run el-ph coupling calculation
    """
    parser = argparse.ArgumentParser() # Create the parser
    parser.add_argument("-q", "--mesh", type=int, default=[8,8,8], nargs=3, help='Defining a mesh grid. (Defaults to [8,8,8]') # Add an argument: mesh
    parser.add_argument("-m", "--mol", type=int, nargs=3, help='The numbering of molecule 1 2 and 3') # Add an argument: mol1
    args = parser.parse_args() # Parse the argument
    
    ut.print_start()

    try:
        run_j0(args.mol) # Create monomers and dimers and calculate J_0 of dimers.
        run_disp_j() # Create displaced dimers and calculate J_ij of dimers.
        run_matrix(args.mesh) # Calculate electron phonon coupling matrix
        ut.print_end()

    except KeyboardInterrupt:
        ut.print_error("Interrupted by user!")
        sys.exit(0)

    


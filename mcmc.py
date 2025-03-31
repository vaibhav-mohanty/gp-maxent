import numpy as np
import math
import argparse
import dill
import tqdm

# function to parse arguments
def parse_args():
    parser = argparse.ArgumentParser(description='EvoTrap script.')
    parser.add_argument('--K',type=int,help='Alphabet size.')
    parser.add_argument('--L',type=int,help='Sequence length.')
    parser.add_argument('--nav',type=int,help='(1) if navigability calculations enabled, (0) if not.')
    parser.add_argument('--theta',type=int,help='(1) if theta calculations enabled, (0) if not.')
    args = parser.parse_args()
    return args 

if __name__ == '__main__':
    # parse arguments
    args = parse_args()

    K = args.K
    L = args.L
    nav = args.nav
    theta = args.theta

    # setup directories
    nav_dir = 'navigability_data'
    energy_dir = 'energy_data'
    theta_dir = 'theta_data'

    # import adjacency table
    with open('adjacency_tables/K' + str(K) + 'L' + str(L) + '.pkl','rb') as file:
        adjtab = dill.load(file)
    
    # setup constants
    J = 1
    Tmin = 0.01
    Tmax = 5.01
    Tsteps = 10
    dT = (Tmax - Tmin)/(Tsteps - 1)

    simtemps = np.arange(Tmin,Tmax,dT)
    
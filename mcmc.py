import numpy as np
import math
import argparse
import dill
from tqdm import tqdm

# function to parse arguments
def parse_args():
    parser = argparse.ArgumentParser(description='EvoTrap script.')
    parser.add_argument('--K',type=int,help='Alphabet size.')
    parser.add_argument('--L',type=int,help='Sequence length.')
    parser.add_argument('--nav',type=int,help='(1) if navigability calculations enabled, (0) if not.')
    parser.add_argument('--theta',type=int,help='(1) if theta calculations enabled, (0) if not.')
    args = parser.parse_args()
    return args

def init_gp_map(K,L,Vn):
    gp_map = np.array([])
    for phenotype in range(0,Q):
        pheno_size = Vn[phenotype]
        new_pheno_vec = phenotype * np.ones((pheno_size))
        gp_map = np.hstack((gp_map,new_pheno_vec))
    return np.array(np.flip(gp_map),dtype=int)

def test_swap(ind1,ind2,Tsim,gp_map,adjtab,Q,J):
    e1_neighbors = adjtab[ind1,:]
    e2_neighbors = adjtab[ind2,:]
    num_neighbors = len(e1_neighbors)

    ind1_vec = ind1*np.ones((num_neighbors),dtype=int)
    ind2_vec = ind2*np.ones((num_neighbors),dtype=int)

    e1_pre = np.sum(gp_map[e1_neighbors] == gp_map[ind1_vec])
    e2_post = np.sum(gp_map[e1_neighbors] == gp_map[ind2_vec])
    e2_pre = np.sum(gp_map[e2_neighbors] == gp_map[ind2_vec])
    e1_post = np.sum(gp_map[e2_neighbors] == gp_map[ind1_vec])

    de = ((e1_post + e2_post) - (e1_pre + e2_pre))
    if ind1 in e2_neighbors and gp_map[ind1] != gp_map[ind2]:  # in case 1 & 2 neighbors
        de -= 2
    de = -J*de

    if de < 0:
        bool_out = True
    elif np.random.rand() < np.exp(-de / Tsim):
        bool_out = True
    else:
        bool_out = False
        de = 0
    
    return bool_out, de

def calc_theta(K,L,Vn,adjtab,gp_map):
    Q = len(Vn)
    theta_mat = np.zeros((Q,Q))

    for i in range(Q):
        temp = np.where(gp_map == i)[0]
        F = len(temp)
        left = temp.T

    return None

    

def mcmc_simulation(K,L,Vn,Tsim,nav_yes=0,theta_yes=0):
    # set up sampled H list and other averages
    Hsampled = []
    if nav_yes:
        navigability_avg = 0
    else:
        navigability_avg = None
    
    if theta_yes:
        theta_avg = np.zeros((Q,Q))
    else:
        theta_avg = None

    transient = 50
    mcs = 50
    time_per_mcs = K**L

    # initialize GP map
    gp_map = init_gp_map(K,L,Vn)
    print(gp_map)

    print('Beginning transient MCMC...')
    for a in tqdm(range(transient)):
        # print('Transient step',a,'of',transient)
        for b in range(time_per_mcs):
            rand_ind1 = np.random.randint(low=0,high=K**L)
            rand_ind2 = np.random.randint(low=0,high=K**L)
            while gp_map[rand_ind1] == gp_map[rand_ind2]:
                rand_ind1 = np.random.randint(low=0,high=K**L)
                rand_ind2 = np.random.randint(low=0,high=K**L)

            bool_out, de = test_swap(rand_ind1,rand_ind2,Tsim,gp_map,adjtab,Q,J)
            if bool_out:
                temp = gp_map[rand_ind2]
                gp_map[rand_ind1] = gp_map[rand_ind2]
                gp_map[rand_ind2] = temp

    theta0 = calc_theta(K,L,Vn,adjtab,gp_map)
    E0 = -J*sum(np.diag(theta0))/2

    print('Beginning main MCMC...')
    for a in tqdm(range(transient)):
        # print('Transient step',a,'of',transient)
        for b in range(time_per_mcs):
            rand_ind1 = np.random.randint(low=0,high=K**L)
            rand_ind2 = np.random.randint(low=0,high=K**L)
            while gp_map[rand_ind1] == gp_map[rand_ind2]:
                rand_ind1 = np.random.randint(low=0,high=K**L)
                rand_ind2 = np.random.randint(low=0,high=K**L)

            bool_out, de = test_swap(rand_ind1,rand_ind2,Tsim,gp_map,adjtab,Q,J)
            if bool_out:
                temp = gp_map[rand_ind2]
                gp_map[rand_ind1] = gp_map[rand_ind2]
                gp_map[rand_ind2] = temp



    # return relevant variables
    return Hsampled, navigability_avg, theta_avg

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

    # setup frequency distribution
    Vn = np.sort(np.hstack(([1],np.outer(np.ones((K-1),dtype=int),np.logspace(0,L-1,num=L,base=K,dtype=int)).flatten())))
    Q = len(Vn)

    for Tsim in simtemps:
        print('Tsim =',Tsim)
        Hsampled, nav, theta = mcmc_simulation(K,L,Vn,Tsim,nav_yes=1,theta_yes=1)
import numpy as np
from utils.indsub import *
import argparse
import dill
from tqdm import tqdm

def parse_args():
    parser = argparse.ArgumentParser(description='EvoTrap script.')
    parser.add_argument('--K',type=int,help='Alphabet size.')
    parser.add_argument('--L',type=int,help='Sequence length.')
    args = parser.parse_args()
    return args 

if __name__ == '__main__':
    # parse arguments
    args = parse_args()

    K = args.K
    L = args.L

    inputdims = K*np.ones((L))

    adjlist = np.zeros((K**L,L*(K-1)),dtype=int)

    for i in tqdm(range(K**L)):
        temp = ind2sub(i,K,L)
        neighbor_vec = -1 * np.ones((L*(K-1)),dtype=int)
        for j in range(L):
            for k in range(K):
                neighbor = temp.copy()
                neighbor[j] = k
                if (temp != neighbor).any():
                    my_next = np.where(neighbor_vec == -1)[0][0]
                    neighbor_vec[my_next] = sub2ind(neighbor,K,L)
        adjlist[i,:] = neighbor_vec.copy()
    
    print(adjlist)
    print('Adjacency list generated.')

    with open('adjacency_tables/K' + str(K) + 'L' + str(L) + '.pkl','wb') as file:
        dill.dump(adjlist,file)

import numpy as np

def viterbi(X, Z, E, T):
    '''Compute the optimal state path
    @param Z all possible state annotation
    @param X observed sequence
    @param E emmision probability
    @param T transition matrix
    '''
    optimal_state_path = []

    # number of hidden states
    k = len(Z)

    # number of observed symbols
    n = len(X)

    # n * k matrix storing the probability of ending at state at certain position
    prob_for_optimum = np.zeros((n+1, k+1))

    # Initialization
    prob_for_optimum[0, 0] = 1
    prob_for_optimum[0, 1:k+1] = 0


    for i in range(1, n+1):
        for j in range(1, k+1):
            v = np.zeros(k+1)


if __name__ == '__main__':

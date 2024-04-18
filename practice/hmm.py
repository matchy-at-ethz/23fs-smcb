import numpy as np
import pandas as pd

from typing import List, Tuple, Iterable

state3 = ["H", "E", "C"]
state8 = ["H", "G", "I", "E", "B", "T", "S", "C"]
map8to3 = {
    "H": "H",
    "G": "H",
    "I": "H",
    "E": "E",
    "B": "C",
    "T": "C",
    "S": "C",
    "C": "C",
}


def MLE(data, aa_states, ss_states):
    # Initialize the parameters
    I = pd.Series(0.0, index=ss_states)

    E = pd.DataFrame(0.0, index=ss_states, columns=aa_states)

    Tr = pd.DataFrame(0.0, index=ss_states, columns=ss_states)

    N = len(data)  # num of data points

    # Iterate over each row of the data
    for i in range(N):
        seq = data.loc[i, "sequence"]
        struct = data.loc[i, "structure"]

        ss_1st = struct[0]
        I[ss_1st] += 1.0

        for j in range(len(seq)):
            aa = seq[j]
            ss = struct[j] if struct[j] in ss_states else map8to3[struct[j]]

            E.loc[ss, aa] += 1.0

            if j < len(seq) - 1:
                ss_next = struct[j + 1] if struct[j + 1] in ss_states else map8to3[struct[j + 1]]

                Tr.loc[ss, ss_next] += 1.0

    # Convert the counts to probs
    I /= I.sum()
    E = E.div(E.sum(axis=1), axis=0)
    Tr = Tr.div(Tr.sum(axis=1), axis=0)

    return {"I": I, "E": E, "Tr": Tr}

def MLE_3state(data, aa_states):
    return MLE(data, aa_states=aa_states, ss_states=state3)

def viterbi(
    E: pd.DataFrame,
    Tr: pd.DataFrame,
    I: np.array,
    data: str | Iterable[str],
    is_log: bool = False,
):
    """
    dp framework for hidden Markov models
    :param E: emission matrix, where E[i, j] is the probability of observing the j-th symbol in the i-th state. index: state, columns: symbols
    :param Tr: transition matrix, where Tr[i, j] is the probability of transitioning from the i-th state to the j-th state. index: from state, columns: to state
    :param I: initial state probabilities, where I[i] is the probability of the i-th state being the initial state. must have I.dtype.names != None
    :param data: observed data, either a string or an iterable of strings (for multiple sequences)
    :param is_log: whether the input is log probabilities
    """

    states = E.index

    if isinstance(data, str):
        data = [data]
    _I = I
    _E = E
    _Tr = Tr
    if not is_log:
        _I = np.log(I)
        _E = np.log(E)
        _Tr = np.log(Tr)

    for seq in data:
        # initialization
        T = len(seq)
        N = len(states)

        # likelihood matrix
        delta = np.zeros((T, N))

        # backpointer matrix
        psi = np.zeros((T, N), dtype=int)

        # base case
        for i, state in enumerate(states):
            delta[0, i] = _I[i] + _E.loc[state, seq[0]]
            psi[0, i] = 0

        # recursive case
        for t in range(1, T):
            for i, state in enumerate(states):
                delta[t, i] = np.max(
                    delta[t - 1, :] + _Tr.loc[:, state] + _E.loc[state, seq[t]]
                )
                psi[t, i] = np.argmax(delta[t - 1, :] + _Tr.loc[:, state])

        # termination
        P = max(delta[T - 1, :])

        # path backtracking
        path = [np.argmax(delta[T - 1, :])]
        for t in range(T - 1, 0, -1):
            path.append(psi[t, path[-1]])
        path.reverse()

        path = "".join([states[i] for i in path])

        yield path, P

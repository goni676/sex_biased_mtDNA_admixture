import numpy as np
from scipy.stats import binom

def theoretical_probability(X, N, p, threshold=0.8):

    # Initial state
    k0 = int(round(p * N))

    # dist[k] = P(K_t = k)
    dist = np.zeros(N + 1)
    dist[k0] = 1.0

    for _ in range(X):

        new_dist = np.zeros(N + 1)

        for k in range(N + 1):

            if dist[k] == 0:
                continue

            probs = binom.pmf(
                np.arange(N + 1),
                N,
                k / N,
            )

            new_dist += dist[k] * probs

        dist = new_dist

    cutoff = int(np.ceil(threshold * N))

    return dist[cutoff:].sum()


def theoretical_probability_by_time(times, N, p_A, threshold=0.8):
    results = {}
    for t in times:
        prob = theoretical_probability(
            t,
            N,
            p_A,
            threshold,
        )


        results[t] = round(prob, 3)

    return results

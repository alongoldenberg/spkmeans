import sys
import numpy as np
import pandas as pd
import myspkmeans
import traceback


def initial_centroids(datapoints, k):
    """
    Parameters
    __________
    datapoints: np.array(n,d)
    k: int

    return the random chosen centroids and their indices based on kmeans++ algorithem.
    """
    index = []
    n, d = datapoints.shape
    centroids = np.zeros((k, d))
    rnd = np.random.choice(n, size=1)
    cur_chosen_m = datapoints[rnd[-1]]
    index.append(rnd[0])
    centroids[0] = cur_chosen_m
    i = 1
    while i < k:
        D = calc_D(centroids, datapoints, i)
        probs = D / D.sum()
        i += 1
        centroids[i - 1], rnd = rnd_select_m(probs, datapoints, n)
        index.append(rnd)
    return centroids, index


def rnd_select_m(probs, datapoints, n):
    idx = np.random.choice(np.arange(n), p=probs)
    return datapoints[idx], idx


def calc_D(Mu, X, i):
    return np.array([min([np.dot(x-m, x-m) for m in Mu[:i+1]]) for x in X])


def print_results(centroids):
    for centroid in centroids:
        print(",".join('%.4f' % x for x in centroid))


def print_eigenvalues(eigenvalues):
    diagonal = np.diag(np.array(eigenvalues))
    diagonal = [round(num, 4) for num in diagonal]
    for i in range(len(diagonal)):
        if diagonal[i] == 0: # find all 0 including -0.0
            diagonal[i] = 0.0
    print(",".join('%.4f' % x for x in diagonal))


def main():
    np.random.seed(0)
    try:
        input_args = sys.argv
        if len(input_args) != 4:
            raise Exception
        k = int(input_args[1])

        if(k == 1):
            print("Invalid Input")
            return
        goal = input_args[2]
        file = input_args[3]
        datapoints = pd.read_csv(file, header=None)
        datapoints = datapoints.to_numpy().tolist()
        n = len(datapoints)
        d = len(datapoints[0])
        if k >= len(datapoints) or k < 0 or k == 1:
            raise Exception
    except:
        print("Invalid Input")
        return
    try:
        if goal == "spk":
            T = myspkmeans.get_goal(n, d, k, "spk", datapoints)
            heuristic_k = len(T[0])
            T = pd.DataFrame(T)
            if k == 0:
                k = heuristic_k
                
            centroids, centroids_index = initial_centroids(T.to_numpy(), k)
            kmeans_new_centroids = myspkmeans.kmeans(n, heuristic_k, k, T.values.tolist(), centroids.tolist())
            print(",".join([str(x) for x in centroids_index]))
            print_results(np.array(kmeans_new_centroids))
        elif goal == "jacobi":
            values_and_vectors = myspkmeans.get_goal(n, d, k, goal, datapoints)
            values = values_and_vectors[0]  # diagonal matrix form
            vectors = values_and_vectors[1]
            print_eigenvalues(values)
            print_results(np.array(vectors))
        else:
            print_results(np.array(myspkmeans.get_goal(n, d, k, goal, datapoints)))
    except:
        print(traceback.format_exc())
        print("An Error Has Occurred")


if __name__ == '__main__':
    main()

import sys
import numpy as np
import pandas as pd
import myspkmeans


def initial_centroids(datapoints, k):
    """
    Parameters
    __________
    datapoints: np.array(n,d)
    k: int
    return the random chosen centroids and their indices based on kmeans++ algorithem.
    """
    index = []
    n, d = datapoints.shape[0], datapoints.shape[1]
    centroids = np.zeros((k, d))
    rnd = np.random.choice(n, size=1)
    cur_chosen_m = datapoints[rnd[-1]]
    index.append(rnd[0])
    centroids[0] = cur_chosen_m
    i = 1
    while i < k:
        distance, distance_sum = calc_distances(centroids, datapoints, i, n)
        probs = distance/distance_sum
        i += 1
        centroids[i - 1], rnd = rnd_select_m(probs.tolist(), datapoints)
        index.append(rnd)
    return centroids, index


def rnd_select_m(probs, datapoints):
    n = len(datapoints)
    idx = np.random.choice(np.arange(n), p=probs)
    return datapoints[idx], idx


def calc_distances(centroids, datapoints, i, n):
    """
    Calculate distances of each datapoint from the centroids.
    """
    distance_lst = []
    distance_sum = 0
    for l in range(n):
        cur_min_distance = sqr_distance(centroids[0], datapoints[l])
        for j in range(i):
            cur_min_distance = min(sqr_distance(centroids[j], datapoints[l]),
                                   cur_min_distance)
        distance_lst.append(cur_min_distance)
        distance_sum += cur_min_distance
    return np.array(distance_lst), distance_sum


def sqr_distance(m, x):
    """
    Calcualte Square distnce between two d-dimensional points.
    """
    return np.dot((m-x), (m-x))


def print_results(centroids):
    for centroid in centroids:
        print(",".join('%.4f' % x for x in centroid))


def print_eigenvalues(eigenvalues):
    diagonal = np.diag(np.array(eigenvalues))
    diagonal = [round(num, 4) for num in diagonal]
    for i in range(len(diagonal)):
        if diagonal[i] == 0:  # find all 0 including -0.0
            diagonal[i] = 0.0
    print(",".join('%.4f' % x for x in diagonal))


def main():
    np.random.seed(0)
    valid_goals = {"wam", "ddg", "lnorm", "jacobi", "spk"}
    try:
        input_args = sys.argv
        if len(input_args) != 4:
            raise Exception
        k = int(input_args[1])
        goal = input_args[2]
        if goal not in valid_goals:
            raise Exception
        file = input_args[3]
        datapoints = pd.read_csv(file, header=None)
        datapoints = datapoints.to_numpy().tolist()
        n = len(datapoints)
        d = len(datapoints[0])
        if (goal == "spk") and (k >= len(datapoints) or k < 0 or k==1):
            raise Exception
    except:
        print("Invalid Input!")
        return
    try:
        if goal == "spk":
            T = myspkmeans.get_goal(n, d, k, "spk", datapoints)
            heuristic_k = len(T[0])
            T = pd.DataFrame(T)
            if k == 0:
                k = heuristic_k
            centroids, centroids_index = initial_centroids(T.to_numpy(), k)
            kmeans_new_centroids = myspkmeans.kmeans(n, heuristic_k, k,
                                                     T.values.tolist(),
                                                     centroids.tolist())
            print(",".join([str(x) for x in centroids_index]))
            print_results(np.array(kmeans_new_centroids))
        elif goal == "jacobi":
            values_and_vectors = myspkmeans.get_goal(n, d, k, goal,
                                                     datapoints)
            values = values_and_vectors[0]  # diagonal matrix form
            vectors = values_and_vectors[1]
            print_eigenvalues(values)
            print_results(np.array(vectors))
        else:
            print_results(np.array(myspkmeans.get_goal(n, d, k, goal,
                                                       datapoints)))
    except:
        print("An Error Has Occurred")


if __name__ == '__main__':
    main()
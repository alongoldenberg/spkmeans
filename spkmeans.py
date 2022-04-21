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
    n, d = datapoints.shape[0], datapoints.shape[1]
    centroids = np.zeros((k, d))
    cur_chosen_m = datapoints[np.random.choice(n, size=1)]
    centroids[0] = cur_chosen_m
    index.append(np.where(np.all(datapoints == cur_chosen_m, axis=1))[0][0])
    i = 1
    while i < k:
        distance_lst, distance_sum = calc_distances(centroids, datapoints, i, n)
        probs_lst = calc_probs(distance_lst, distance_sum)
        i += 1
        centroids[i - 1] = rnd_select_m(probs_lst, datapoints)
        index.append(np.where(np.all(datapoints == centroids[i - 1], axis=1))[0][0])
    return centroids, index


def rnd_select_m(probs, datapoints):
    n = len(datapoints)
    idx = np.random.choice(np.arange(n), p=probs)
    return datapoints[idx]


def calc_cumulative_array(arr):
    res = [0 for i in range(len(arr))]
    res[0] = arr[0]
    j = 1
    for num in arr[1:]:
        res[j] = res[j - 1] + num
    return res


def calc_distances(centroids, datapoints, i, n):
    """
    Calculate distances of each datapoint from the centroids.
    """
    distance_lst = []
    distance_sum = 0
    for l in range(n):
        cur_min_distance = sqr_distance(centroids[0], datapoints[l])
        for j in range(i):
            cur_min_distance = min(sqr_distance(centroids[j], datapoints[l]), cur_min_distance)
        distance_lst.append(cur_min_distance)
        distance_sum += cur_min_distance
    return distance_lst, distance_sum


def calc_probs(distance_lst, distance_sum):
    """
    :param distance_lst: list[int]
    :param distance_sum: int
    :return: The distributivity function of datapoints selection
    """
    res = []
    n = len(distance_lst)
    for l in range(n):
        res.append(distance_lst[l] / distance_sum)
    return res


def sqr_distance(m, x):
    """
    Calcualte Square distnce between two d-dimensional points.
    """
    res = 0
    for i in range(len(m)):
        res += ((m[i] - x[i]) ** 2)
    return res


def print_results(centroids):
    for centroid in centroids:
        print(",".join('%.4f' % x for x in centroid))


def print_eigenvalues(eigenvalues):
    diagonal = np.diag(np.array(eigenvalues))
    for i in range(len(diagonal)):
        if diagonal[i] == -0.0000:
            diagonal[i] = 0.0000
    print(",".join('%.4f' % x for x in diagonal))


def main():
    np.random.seed(0)
    try:
        input_args = sys.argv
        if len(input_args) != 4:
            raise Exception
        k = int(input_args[1])

        if(k == 1):
            print("Invalid Input!")
            return
        goal = input_args[2]
        file = input_args[3]
        datapoints = pd.read_csv(file, header=None)
        datapoints = datapoints.to_numpy().tolist()
        n = len(datapoints)
        d = len(datapoints[0])
        if k >= len(datapoints) or k < 0:
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
            real_index = [int(T.iloc[i].name) for i in centroids_index]
            kmeans_new_centroids = myspkmeans.kmeans(n, heuristic_k, k, T.values.tolist(), centroids.tolist())
            print(','.join(map(str, real_index)))
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
        print("An Error Has Occurred!")


if __name__ == '__main__':
    main()

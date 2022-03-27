import sys
import numpy as np
import pandas as pd
import mykmeanssp
from enum import Enum


def main():
    try:
        input_args = sys.argv
        if len(input_args) != 4:
            raise Exception
            k = int(input_args[1])
            goal = int(input_args[2])
            file = input_args[4]
        datapoints = pd.read_csv(file1, index_col=0, header=None).round(4)
        if k >= len(datapoints) or k < 0: raise Exception
    except:
        print("Invalid Input!")
        return

    try:
        init_centroids, centroids_index = initial_centroids(datapoints.to_numpy(), k)
        real_index = [int(datapoints.iloc[i].name) for i in centroids_index]
        centroids = mykmeanssp.fit(datapoints.values.tolist(),
                                   init_centroids.tolist(),
                                   k,
                                   max_iter,
                                   epsilon,
                                   len(datapoints),
                                   len(datapoints.columns))
        centroids = np.round(centroids, 4)
        write_results(real_index, centroids)
    except:
        print(traceback.format_exc())
        print("An Error Has Occurred!")


if __name__ == '__main__':
    main()
# Importing necessary libraries
import argparse
from typing import List
from typing import Optional
import numpy as np
import laspy

def main():
    
    # Parsing the input arguments with ArgumentParser
    parser = argparse.ArgumentParser(
        "LAZ cropper", description="Outputs a LAZ file that is a subset of the input file"
    )
    parser.add_argument("ifile", help="input file in LAZ")
    parser.add_argument("ofile", help="output file in LAZ")
    parser.add_argument("x_min", type=float, help="the longitude of the southwest corner of the output area (EPSG:28992)")
    parser.add_argument("y_min", type=float, help="the latitude of the southwest corner of the output area (EPSG:28992)")
    parser.add_argument("x_size", type=float, help="the west-east size of the output area, in meters")
    parser.add_argument("y_size", type=float, help="the north-south size of the output area, in meters")
    parser.add_argument("--points-per-iter", default=10**6, type=int, help="OPTIONAL: the number of points processed per iteration of the chunk_iterator, by default it is 10^6")

    # Saving parsed arguments to variables
    args = parser.parse_args()
    input_file = args.ifile
    output_file = args.ofile
    x_min = args.x_min
    y_min = args.y_min
    x_size = args.x_size
    y_size = args.y_size
    points_per_iter = args.points_per_iter

    # Reading the .LAZ file
    with laspy.open(input_file) as file:

        # Defining the desired bounds
        x_max = x_min + x_size
        y_max = y_min + y_size
        sub_bounds = [(x_min,y_min,x_max,y_max)]

        # Using the method given in laspy documentation to create a new .LAZ file containing just the desired data points
        writers: List[Optional[laspy.LasWriter]] = [None] * len(sub_bounds)
        try:
            count = 0
            for points in file.chunk_iterator(points_per_iter):
                print(f"{count / file.header.point_count * 100}%")

                # For performance we need to use copy so that the underlying arrays are contiguous
                x, y = points.x.copy(), points.y.copy()

                point_piped = 0
                # For each sub bound
                for i, (x_min, y_min, x_max, y_max) in enumerate(sub_bounds):
                    mask = (x >= x_min) & (x <= x_max) & (y >= y_min) & (y <= y_max)
                    # If there are points in the sub bound
                    if np.any(mask):
                        if writers[i] is None:
                            writers[i] = laspy.open(
                                output_file, mode="w", header=file.header
                            )
                        sub_points = points[mask]
                        writers[i].write_points(sub_points)

                    # If there are no points in the sub bound
                    point_piped += np.sum(mask)
                    if point_piped == len(points):
                        break
                count += len(points)
            print(f"{count / file.header.point_count * 100}%")

        # If there is an error, close all the writers
        finally:
            for writer in writers:
                if writer is not None:
                    writer.close()

if __name__ == "__main__":
    main()
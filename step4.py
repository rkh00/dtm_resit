# Importing necessary libraries
import argparse
from typing import List
from typing import Optional
import numpy as np
import laspy
import startinpy

def main():
    
    # Parsing the input arguments with ArgumentParser
    parser = argparse.ArgumentParser(
        "LAZ TIN simplifier", description="Performs TIN decimation on input LAZ file"
    )
    parser.add_argument("ifile", help="input file in LAZ")
    parser.add_argument("ofile", help="output file in LAZ")
    parser.add_argument("threshold", type=float, help="the vertical distance threshold used for the decimation procedure, in meters")
    parser.add_argument("--points-per-iter", default=10**6, type=int, help="OPTIONAL: the number of points processed per iteration of the chunk_iterator, by default it is 10^6")

    # Saving parsed arguments to variables
    args = parser.parse_args()
    input_file = args.ifile
    output_file = args.ofile
    threshold = args.threshold
    points_per_iter = args.points_per_iter
    
    print("Performing TIN decimation with threshold = " + str(threshold) + "m")
    print("Reading .LAZ file...")
    
    # Initializing TIN containing the input points
    dt = startinpy.DT()
    dt.read_las(f"{input_file}")
    
    continue_simplification = True
    iteration = 0
    importances = np.zeros(len(dt.points))
    first_iter = True
    neighbors = []
    
    print("Performing algorithm...")
    
    # The main loop which simplifies the TIN
    while continue_simplification == True:
        
        # If it's the first iteration...
        if first_iter == True:
            # ...we want to find the importance of all the points
            for i, point in enumerate(dt.points):
                
                # Skip the vertex at infinity
                if i == 0:
                    continue
                
                # Find the coordinates of the given point
                x = point[0]
                y = point[1]
                z_actual = point[2]
                
                # Remove the point from the TIN temporarily
                dt.remove(i)
                
                try:
                    # Find the height of the projected point
                    z_temp = dt.interpolate_tin_linear(x, y)
                    
                except:
                    # If this is not possible because the point is part of the
                    # convex hull, add the point back into the TIN and move
                    # on to the next point
                    dt.insert_one_pt(x, y, z_actual)
                    continue
                
                # Find the importance and add it to the importances array
                d = np.absolute(z_actual - z_temp)
                importances[i] = d
                
                # Finally, add the point back into the TIN
                dt.insert_one_pt(x, y, z_actual)
                
            # We don't want to have to iterate over all the points again
            first_iter = False
            
        # If it's not the first iteration...
        else:
            # ...then only find the importance for the neighboring points to
            # the previous deleted point
            for k in range(len(neighbors)):
                
                # Find the coordinates of the given point
                x = dt.get_point(neighbors[k])[0]
                y = dt.get_point(neighbors[k])[1]
                z_actual = dt.get_point(neighbors[k])[2]
                
                # Remove the point from the TIN temporarily
                dt.remove(neighbors[k])
                
                try:
                    # Find the height of the projected point
                    z_temp = dt.interpolate_tin_linear(x, y)
                    
                except:
                    # If this is not possible because the point is part of the
                    # convex hull, add the point back into the TIN and move
                    # on to the next point
                    dt.insert_one_pt(x, y, z_actual)
                    continue
                
                # Find the importance and add it to the importances array
                d = np.absolute(z_actual - z_temp)
                importances[neighbors[k]] = d
                
                # Finally, add the point back into the TIN
                dt.insert_one_pt(x, y, z_actual)
                
        # Find the index of the least important point
        to_remove = np.argmax(importances)
        
        # If this importance is less than the user-defined threshold...
        if importances[to_remove] < threshold:
            # ...terminate the algorithm
            continue_simplification = False
            
        else:
            # Find the neighbors of the point we are about to remove
            neighbors = dt.adjacent_vertices_to_vertex(to_remove)
                    
            # Finally, remove the point
            dt.remove(to_remove)
            
            # Set the importance of the removed point to zero so we don't
            # try to remove it again
            importances[to_remove] = 0
            
        # Just for user confirmation
        iteration += 1
        if iteration % 1000 == 0:
            print("Removed " + str(iteration) + " points")

    print("done")
    print("Saving simplified TIN to .LAZ...")

    # Reading the .LAZ file
    with laspy.open(input_file) as ifile:
        
        sub_bounds = [(ifile.header.x_min,ifile.header.y_min,ifile.header.x_max,ifile.header.y_max)]
        
        dt_points = dt.points
        dt_points_set = {arr.tobytes() for arr in dt_points}

        # Using the method given in laspy documentation to create a new .LAZ file containing just the desired data points
        writers: List[Optional[laspy.LasWriter]] = [None] * (int(np.ceil(ifile.header.point_count/points_per_iter)))
        try:
            count = 0
            # i = 0
            for points in ifile.chunk_iterator(points_per_iter):
                print(f"{count / ifile.header.point_count * 100}%")

                # For performance we need to use copy so that the underlying arrays are contiguous
                x, y, z = points.x.copy(), points.y.copy(), points.z.copy()

                cloud_pts = np.stack((x,y,z),axis=1)
                mask = [arr.tobytes() in dt_points_set for arr in cloud_pts]

                point_piped = 0
                # For each sub bound
                for i, (x_min, y_min, x_max, y_max) in enumerate(sub_bounds):
                    # If there are points in the sub bound
                    if np.any(mask):
                        if writers[i] is None:
                            writers[i] = laspy.open(
                                output_file, mode="w", header=ifile.header
                            )
                        sub_points = points[mask]
                        writers[i].write_points(sub_points)

                count += len(points)
                # i += 1

                # If there are no points in the sub bound
                point_piped += np.sum(mask)
                if point_piped == len(points):
                    break
            print(f"{count / ifile.header.point_count * 100}%")

        # If there is an error, close all the writers
        finally:
            for writer in writers:
                if writer is not None:
                    writer.close()

if __name__ == "__main__":
    main()
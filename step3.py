# Importing necessary libraries
import argparse
from typing import List
from typing import Optional
import numpy as np
import laspy
import startinpy

def equation_plane(x1, y1, z1, x2, y2, z2, x3, y3, z3):
    """Calculates the equation of the plane spanned by three points in 3D space.
    Parameters:
        The coordinates of three points (x1,y1,z1), (x2,y2,z2), (z1,z2,z3).
    Returns:
        The coefficients a, b, c, and d corresponding to the planar equation 
        ax + by + cz + d = 0."""
    
    a1 = x2 - x1
    b1 = y2 - y1
    c1 = z2 - z1
    a2 = x3 - x1
    b2 = y3 - y1
    c2 = z3 - z1
    a = b1 * c2 - b2 * c1
    b = a2 * c1 - a1 * c2
    c = a1 * b2 - b1 * a2
    d = (- a * x1 - b * y1 - c * z1)
    
    return a, b, c, d

def shortest_distance(x1, y1, z1, a, b, c, d):
    """Calculates the shortest distance between a point and a plane in 3D space.
    Parameters:
        The coordinates of a point (x1,y1,z1) and the coefficients a, b, c, and 
        d corresponding to the planar equation ax + by + cz + d = 0.
    Returns:
        The shortest distance between the input point and the input plane."""
    
    d = abs((a * x1 + b * y1 + c * z1 + d))
    e = (np.sqrt(a * a + b * b + c * c))
    
    return d/e

def main():
    
    # Parsing the input arguments with ArgumentParser
    parser = argparse.ArgumentParser(
        "LAZ TIN refinement processor", description="Performs TIN refinement on the input point cloud"
    )
    parser.add_argument("ifile", help="input file in LAZ")
    parser.add_argument("ofile", help="output file in LAZ")
    parser.add_argument("cell_size", type=float, help="the size of each cell used in the rudimentary TIN construction, in meters")
    parser.add_argument("max_distance", type=float, help="perpendicular distance between p and the plane spanned by the triangle, in meters")
    parser.add_argument("max_angle", type=float, help="the largest angle, in degrees, of the angles between the triangle and the three vectors that connect each vertex with p")
    parser.add_argument("--points-per-iter", default=10**6, type=int, help="OPTIONAL: the number of points processed per iteration of the chunk_iterator, by default it is 10^6")

    # Saving parsed arguments to variables
    args = parser.parse_args()
    input_file = args.ifile
    output_file = args.ofile
    cell_size = args.cell_size
    max_distance = args.max_distance
    max_angle = args.max_angle
    points_per_iter = args.points_per_iter
    
    # Confirmation for the user
    print("Performing TIN refinement with cell_size = " + str(cell_size) + "m, d_max = " + str(max_distance) + "m, and alpha_max = " + str(max_angle) + "°")
    
    # Extracting points for rudimentary TIN
    with laspy.open(input_file) as ifile:
        
        sub_bounds = []
    
        # Calculating the number of rows and columns, given cell_size
        rows = int(np.ceil((float(ifile.header.x_max) - float(ifile.header.x_min))/cell_size))
        cols = int(np.ceil((float(ifile.header.y_max) - float(ifile.header.y_min))/cell_size))
        
        # Defining arrays for construction rudimentary TIN
        lowest_points = [None]*(rows*cols)
        lowest_heights = np.ones(rows*cols)*999999
        
        # Iterating over the grid and adding the bounds to sub_bounds
        for row in range(rows):
        	for col in range(cols):
        		x_min = ifile.header.x_min + cell_size*col
        		x_max = x_min + cell_size
        		y_min = ifile.header.y_min + cell_size*row
        		y_max = y_min + cell_size
        		bounds = (x_min,y_min,x_max,y_max)
        		sub_bounds.append(bounds)
        
        count = 0
        print("Extracting ground points for rudimentary TIN...")
        
        # Iterating over the input point cloud, one chunk at a time
        for points in ifile.chunk_iterator(points_per_iter):
            print(f"{count / ifile.header.point_count * 100}%")

            # For performance we need to use copy so that the underlying arrays are contiguous
            x, y, z = points.x.copy(), points.y.copy(), points.z.copy()

            # For each sub bound
            for i, (x_min, y_min, x_max, y_max) in enumerate(sub_bounds):
                mask = (x >= x_min) & (x <= x_max) & (y >= y_min) & (y <= y_max)
                # If there are points in the sub bound
                if np.any(mask):
                    sub_points = points[mask]
                    x_array = np.array(sub_points.x)
                    y_array = np.array(sub_points.y)
                    z_array = np.array(sub_points.z)
                    
                    # Finding the index of the lowest point in this set of
                    # points, as well as that point's height value
                    lowest_point = np.argmin(z_array)
                    lowest_height = z_array[lowest_point]
                    
                    # If the lowest point is lower than the current lowest
                    # point for that cell, replace it
                    if lowest_height < lowest_heights[i]:
                        lowest_points[i] = [x_array[lowest_point],y_array[lowest_point],z_array[lowest_point]]
                        lowest_heights[i] = z_array[lowest_point]
                        
            count += len(points)
            
        print("done")
    
    print("Creating rudimentary TIN...")
    
    # Initializing the TIN
    dt = startinpy.DT()
    
    # Adding the lowest point for each cell to the TIN
    for j in range(len(lowest_points)):
        x_to_add = lowest_points[j][0]
        y_to_add = lowest_points[j][1]
        z_to_add = lowest_points[j][2]
        dt.insert_one_pt(x_to_add, y_to_add, z_to_add)
    print("done")
    
    print("Performing greedy algorithm. This will take some time! Best to leave it running overnight")
    
    continue_refinement = True     
    iteration = 0

    # The main loop in which all the non-ground points are filtered out
    while continue_refinement == True:
        
        len_before = len(dt.points)
        
        with laspy.open(input_file) as file:
            
            for points in file.chunk_iterator(points_per_iter):
    
                # For performance we need to use copy so that the underlying arrays are contiguous
                x, y, z = points.x.copy(), points.y.copy(), points.z.copy()
    
                # For each sub bound
                for i, (x_min, y_min, x_max, y_max) in enumerate(sub_bounds):
                    mask = (x >= x_min) & (x <= x_max) & (y >= y_min) & (y <= y_max)
                    
                    # If there are points in the sub bound
                    if np.any(mask):
                        sub_x = x[mask]
                        sub_y = y[mask]
                        sub_z = z[mask]
                        
                        # Finding the points from the current TIN which are
                        # inside the current cell
                        dt_points = dt.points
                        dt_local = dt_points[(dt_points[:,0] >= x_min) & (dt_points[:,0] <= x_max) & (dt_points[:,1] >= y_min) & (dt_points[:,1] <= y_max)]
                        
                        if np.any(dt_local):
                        
                            # Finding the highest point in the cell
                            zmax_local = np.max(dt_local[:,-1])
                            
                            # Filtering out all points which are at least d_max
                            # from the current TIN
                            mask2 = [sub_z[i] <= zmax_local + max_distance for i in range(len(sub_z))]
                            sub_x = sub_x[mask2]
                            sub_y = sub_y[mask2]
                            sub_z = sub_z[mask2]
                        
                            # Filtering out all points that are already in the
                            # TIN
                            dt_points_set = {arr.tobytes() for arr in dt_points}    
                            cloud_pts = np.stack((sub_x,sub_y,sub_z),axis=1)
                            mask3 = [arr.tobytes() not in dt_points_set for arr in cloud_pts]
                            sub_x = sub_x[mask3]
                            sub_y = sub_y[mask3]
                            sub_z = sub_z[mask3]
                            
                            for l in range(len(sub_x)):
                                try:
                                    # Finding the triangle in the TIN that the point's projection
                                    # falls within
                                    triangle = dt.locate(sub_x[l],sub_y[l])
                                    
                                except:
                                    # Error handling: If no such triangle exists, do nothing and
                                    # process the next point in the loop
                                    continue
                                
                                # Finding the three vertices that define the triangle
                                v1 = dt.get_point(triangle[0])
                                v2 = dt.get_point(triangle[1])
                                v3 = dt.get_point(triangle[2])
                                
                                # Finding the plane spanned by the triangle
                                a, b, c, d = equation_plane(v1[0],v1[1],v1[2],v2[0],v2[1],v2[2],v3[0],v3[1],v3[2])
                                
                                # Finding the shortest distance between the input point and the 
                                # plane spanned by the triangle
                                distance = shortest_distance(sub_x[l],sub_y[l],sub_z[l],a,b,c,d)
                                
                                # Calculating the angles between the triangle and the three vectors 
                                # that connect each vertex with the input point
                                alpha1 = np.rad2deg(np.arcsin(distance/np.sqrt((sub_x[l]-v1[0])**2 + (sub_y[l]-v1[1])**2 + (sub_z[l]-v1[2])**2)))
                                alpha2 = np.rad2deg(np.arcsin(distance/np.sqrt((sub_x[l]-v2[0])**2 + (sub_y[l]-v2[1])**2 + (sub_z[l]-v2[2])**2)))
                                alpha3 = np.rad2deg(np.arcsin(distance/np.sqrt((sub_x[l]-v3[0])**2 + (sub_y[l]-v3[1])**2 + (sub_z[l]-v3[2])**2)))
                                
                                # Finding the largest of these angles
                                alphas = np.array([alpha1,alpha2,alpha3])
                                max_alpha = np.max(alphas)
                                
                                # Only if the distance and the largest angle fall within the
                                # desired tolerances to we add the point to the TIN
                                if distance < max_distance and max_alpha < max_angle:
                                    dt.insert_one_pt(sub_x[l], sub_y[l], sub_z[l])
            
        # If this iteration hasn't added any points, terminate the algorithm
        len_after = len(dt.points)
        added_points = len_after - len_before
        if added_points == 0:
            continue_refinement = False
        
        iteration += 1
        print("Iteration #" + str(iteration) + " done. " + str(added_points) + " point(s) added to TIN")
        print(f"{dt.number_of_vertices() / file.header.point_count * 100}% of original points classified as ground points")
    
    print("Saving refined TIN to output .LAZ file...")
    
    dt_points = dt.points
    dt_points_set = {arr.tobytes() for arr in dt_points}
    
    with laspy.open(input_file) as ifile:
        # Writing the points to a new point cloud
        writers: List[Optional[laspy.LasWriter]] = [None] * len(sub_bounds)
        try:
            for points in ifile.chunk_iterator(points_per_iter):
    
                # For performance we need to use copy so that the underlying arrays are contiguous
                x, y, z = points.x.copy(), points.y.copy(), points.z.copy()
    
                # Adding only the points that are in the TIN
                cloud_pts = np.stack((x,y,z),axis=1)
                mask = [arr.tobytes() in dt_points_set for arr in cloud_pts]
                
                for i, (x_min, y_min, x_max, y_max) in enumerate(sub_bounds):
                    if np.any(mask):
                        if writers[i] is None:
                            writers[i] = laspy.open(
                                output_file, mode="w", header=ifile.header
                            )
                        sub_points = points[mask]
                        writers[i].write_points(sub_points)
    
        # If there is an error, close all the writers
        finally:
            for writer in writers:
                if writer is not None:
                    writer.close()
    print("done")
    
if __name__ == "__main__":
    main()
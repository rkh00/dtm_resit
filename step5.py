# Importing necessary libraries
import argparse
import numpy as np
import startinpy
import rasterio

def distance_2d(x1, y1, x2, y2):
    """Calculates the distance between two points projected to a 2D plane.
    Parameters:
        The coordinates of two points (x1,y1), (x2,y2).
    Returns:
        The 2D distance d between the two points."""
    d = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    return d

def laplace_xy(dt, x, y):
    """Finds the interpolated height of an input point using Laplace interpolation.
    Parameters:
        The startinpy TIN dt and the coordinates of the point (x,y=.
    Returns:
        The interpolated height of the point (x,y), according to Laplace
        interpolation."""
    
    # If the point is outside the convex hull then it cannot be interpolated
    if not dt.is_inside_convex_hull(x, y):
        raise Exception("Outside convex hull")
    
    # Temporarily add the point to be interpolated to the TIN
    new_point = dt.insert_one_pt(x,y,0)
    
    # Find the neighboring triangles and vertices to the point
    neighboring_triangles = dt.incident_triangles_to_vertex(new_point)
    neighboring_vertices = dt.adjacent_vertices_to_vertex(new_point)
    
    # Initializing heights and weights
    weights = np.zeros(len(neighboring_vertices))
    heights = np.zeros(len(neighboring_vertices))
    
    # The main loop used to interpolate the height of the point. We iterate
    # over every neighboring vertex
    for i, vertex in enumerate(neighboring_vertices):
        
        # Insert the height value of this neighboring vertex into the heights
        # array
        heights[i] = dt.get_point(vertex)[2]
        
        # Find the "overlapping" neighboring triangles between the point to
        # be interpolated and this neighboring vertex
        neighboring_triangles_of_neighboring_vertex = dt.incident_triangles_to_vertex(vertex)
        relevant_mask = [relevant in neighboring_triangles_of_neighboring_vertex for relevant in neighboring_triangles]
        relevant_triangles = neighboring_triangles[relevant_mask]
        
        centers = []
        
        # Find the center of the circumcircle of both triangles in the
        # "overlap"
        for triangle in relevant_triangles:
            x1,y1 = dt.get_point(triangle[0])[0],dt.get_point(triangle[0])[1]
            x2,y2 = dt.get_point(triangle[1])[0],dt.get_point(triangle[1])[1]
            x3,y3 = dt.get_point(triangle[2])[0],dt.get_point(triangle[2])[1]
            a = np.array([[2*(x1-x2),2*(y1-y2)],[2*(x1-x3),2*(y1-y3)]])
            b = np.array([x1**2+y1**2-x2**2-y2**2,x1**2+y1**2-x3**2-y3**2])
            center = np.linalg.solve(a,b)
            centers.append(center)
            
        # Find the distance between these centers and add it to the weights
        # array
        weight = distance_2d(centers[0][0],centers[0][1],centers[1][0],centers[1][1])/distance_2d(x, y, dt.get_point(vertex)[0], dt.get_point(vertex)[1])
        weights[i] = weight
    
    # Finally, calculate the interpolated height of the point, given by the
    # definition of a weighted-average interpolation method
    interpolated_height = np.sum(np.multiply(weights, heights))/np.sum(weights)
    
    # Remove the point to be interpolated from the TIN
    dt.remove(new_point)
    
    return interpolated_height

def main():
    
    # Parsing the input arguments with ArgumentParser
    parser = argparse.ArgumentParser(
        "LAZ Laplace interpolator", description="Performs Laplace interpolation on an input .LAZ point cloud file to output a GeoTIFF raster with interpolated points"
    )
    parser.add_argument("ifile", help="input file in LAZ")
    parser.add_argument("ofile", help="output file in GeoTIFF")
    parser.add_argument("pixel_size", type=float, help="the size of each cell in the output raster, in meters")

    # Saving parsed arguments to variables
    args = parser.parse_args()
    input_file = args.ifile
    output_file = args.ofile
    pixel_size = args.pixel_size
    
    # Initializing TIN containing the input points and removing the infinite
    # vertex
    print("Reading .LAZ file...")
    dt = startinpy.DT()
    dt.read_las(f"{input_file}")
    dt_points_no_inf = np.array(dt.points)[1:]
    
    # Finding the bounding box of the dataset
    x_min = np.min(dt_points_no_inf[:,0])
    x_max = np.max(dt_points_no_inf[:,0])
    y_min = np.min(dt_points_no_inf[:,1])
    y_max = np.max(dt_points_no_inf[:,1])
    
    # Initializing a 2D numpy array to store the raster
    x = np.arange(x_min, x_max, pixel_size)
    y = np.arange(y_min, y_max, pixel_size)
    z = np.zeros([len(y),len(x)])
    
    # Looping over the raster array
    for j in range(len(z)):
        print("row " + str(j))
        
        for k in range(len(z[0])):
            
            # Try to interpolate the height of this point
            try:
                z[j][k] = laplace_xy(dt,x[k],y[j])
                
            # If this does not work, simply set the value to NaN and move on
            except:
                z[j][k] = np.nan
                continue
    
    # Some rasterio Affine transform stuff
    xres = (x_max - x_min) / len(x)
    yres = (y_max - y_min) / len(y)
    transform = rasterio.Affine.translation(x[0] - xres / 2, y[0] - yres / 2) * rasterio.Affine.scale(xres, yres)
    
    # Saving the output raster
    with rasterio.open(
        f'{output_file}',
        'w',
        driver='GTiff',
        height=z.shape[0],
        width=z.shape[1],
        count=1,
        dtype=z.dtype,
        crs='+proj=latlong',
        transform=transform,
        ) as dst:
        dst.write(z,1)

if __name__ == "__main__":
    main()
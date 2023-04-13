# Importing necessary libraries
import argparse
import numpy as np
import rasterio

def main():
    
    # Parsing the input arguments with ArgumentParser
    parser = argparse.ArgumentParser(
        "GeoTIFF hillshade visualizer", description="Performs hillshading on an input GeoTIFF raster to produce an output GeoTIFF raster where the terrain has been visualized"
    )
    parser.add_argument("ifile", help="input file in GeoTIFF")
    parser.add_argument("ofile", help="output file in GeoTIFF")
    parser.add_argument("--sun_height", default=45.0, type=float, help="OPTIONAL: the height of the sun above the horizon, in degrees. By default it is 45")
    parser.add_argument("--sun_azimuth", default=315.0, type=float, help="the angle of the sun clockwise from the north, in degrees. By default it is 315 (northwest)")
    
    # Saving parsed arguments to variables
    args = parser.parse_args()
    input_file = args.ifile
    output_file = args.ofile
    sun_height = np.deg2rad(args.sun_height)
    sun_azimuth = np.deg2rad(args.sun_azimuth)
    
    # Opening the input raster and reading the data
    dataset = rasterio.open(input_file)
    n1 = dataset.read(1)
    
    # Creating a blank array with the same shape as the input raster
    hillshade = np.zeros_like(n1)
    
    # The index of the last row and column
    last_i = len(hillshade) - 1
    last_j = len(hillshade[0]) - 1
    
    # The size of a pixel
    r = (dataset.transform*(1,0))[0] - (dataset.transform*(0,0))[0]
    
    # The main loop, iterating over the entire 2D hillshade array
    for i in range(len(hillshade)):
        print("row " + str(i))
        
        for j in range(len(hillshade[0])):
            
            # If there is no data in this pixel, set it to no_data and move on
            if n1[i][j] == dataset.nodatavals:
                hillshade[i][j] = np.nan
                continue
            
            # Handling exceptions and calculating dy
            if i == 0:
                dy = n1[i+1][j]/r
            elif i == last_i:
                dy = n1[i-1][j]/r
            else:
                dy = (n1[i-1][j] - n1[i+1][j])/(2*r)
            
            # Handling exceptions and calculating dx
            if j == 0:
                dx = n1[i][j+1]/r
            elif j == last_j:
                dx = n1[i][j-1]/r
            else:
                dx = (n1[i][j-1] - n1[i][j+1])/(2*r)

            # Calculating the gradient and aspect for this pixel
            gradient = np.arctan(np.sqrt(dx**2 + dy**2))
            aspect = np.arctan(dy/dx)
            
            # Calculating the hillshade for this pixel
            hillshade[i][j] = 250*((np.cos(sun_height)*np.cos(aspect)) + (np.sin(sun_height)*np.sin(aspect)*np.cos(sun_azimuth - gradient)))
            
    # Saving the output raster
    with rasterio.open(
        f'{output_file}',
        'w',
        driver='GTiff',
        height=hillshade.shape[0],
        width=hillshade.shape[1],
        count=1,
        dtype=hillshade.dtype,
        crs='+proj=latlong',
        transform=dataset.transform,
        ) as dst:
        dst.write(hillshade,1)

if __name__ == "__main__":
    main()
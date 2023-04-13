# DTM re-sit assignment

The assignment in question can be found [here](https://3d.bk.tudelft.nl/courses/geo1015/hw/resit/).  

## step2.py

Requires the following libraries:  
- numpy  
- laspy  

The program takes in a .LAZ point cloud file and outputs a .LAZ file which contains a subset of this point cloud, given by the boundaries specified by the user.  

The code is run with the following command:  
```
python step2.py input.LAZ output.LAZ longitude latitude x_size y_size [points_per_iter]
```

- input.LAZ: The input point cloud. This must be in the same folder as the program. Replace "input" with the name of the .LAZ file, for example C_68GZ1.LAZ  
- output.LAZ: The output point cloud. This will appear in the same folder as the program. Replace "output" with the name of the .LAZ file, for example step2.LAZ  
- longitude: The longitude of the southwest corner of the output area (given in EPSG:28992). Must be an int. Example: 191135  
- latitude: The latitude of the southwest corner of the output area (given in EPSG:28992). Must be an int. Example: 325200  
- x_size: Number of meters east from the given longitude to include in the output area. Must be an int. Example: 500  
- y_size: Number of meters north from the given latitude to include in the output area. Must be an int. Example: 500  
- *points\_per\_iter: Optional. The number of points processed per iteration of the chunk_iterator. Must be an int. By default set to 10^6. Increase this if you are feeling adventuruous*  

## step3.py

Requires the following libraries:  
- numpy  
- laspy  
- startinpy  

The program takes in a .LAZ point cloud file and outputs a .LAZ file which contains only the points classified as ground points, given ground classifiers specified by the user.  

The code is run with the following command:  
```
python step3.py input.LAZ output.LAZ cell_size max_distance max_angle [points_per_iter]
```

- input.LAZ: The input point cloud. This must be in the same folder as the program. Replace "input" with the name of the .LAZ file, for example step2.LAZ  
- output.LAZ: The output point cloud. This will appear in the same folder as the program. Replace "output" with the name of the .LAZ file, for example step3.LAZ  
- cell\_size: The program first constructs a rudimentary TIN consisting of the lowest points in each of a set of smaller cells that make up the dataset. This parameter specifies the size of those cells, in meters. Example: 35  
- max\_distance: For each iteration the greedy algorithm, a calculation is performed for each point in the input point cloud that is not already in the TIN. This parameter determines the distance d_max in meters from the TIN which the point must fall within in order to be added to the TIN. Example: 1.5  
- max\_angle: For each iteration the greedy algorithm, a calculation is performed for each point in the input point cloud that is not already in the TIN. This parameter determines the angle alpha_max in degrees between the relevant triangle in the TIN and the point, seen from each of the triangle's vertices. All angles must fall within this value in order to be added to the TIN. Example: 5    
- *points\_per\_iter: Optional. The number of points processed per iteration of the chunk_iterator. Must be an int. By default set to 10^6. Increase this if you are feeling adventuruous*  

**NB: This program is slow!** Depending on the parameters given, it could take several hours to finish.  

## step4.py

Requires the following libraries:  
- numpy  
- laspy  
- startinpy  

The program takes in a .LAZ point cloud file and outputs a .LAZ file which is a simplified version of the point cloud, given a simplification parameter specified by the user.  

The code is run with the following command:  
```
python step4.py input.LAZ output.LAZ threshold [points_per_iter]
```

- input.LAZ: The input point cloud. This must be in the same folder as the program. Replace "input" with the name of the .LAZ file, for example step3.LAZ  
- output.LAZ: The output point cloud. This will appear in the same folder as the program. Replace "output" with the name of the .LAZ file, for example step4.LAZ  
- threshold: The vertical distance threshold used for the decimation procedure, in meters. Example: 0.1  
- *points\_per\_iter: Optional. The number of points processed per iteration of the chunk_iterator. Must be an int. By default set to 10^6. Increase this if you are feeling adventuruous*  

## step5.py

Requires the following libraries:  
- numpy  
- startinpy  
- rasterio  

The program takes in a .LAZ point cloud file and outputs a .tif raster file containing a rasterized version of the point cloud file, where the values of the raster's pixels are interpolated using Laplace interpolation. The size of the pixels are decided by the user.  

The code is run with the following command:  
```
python step5.py input.LAZ output.tif pixel_size
```

- input.LAZ: The input point cloud. This must be in the same folder as the program. Replace "input" with the name of the .LAZ file, for example step4.LAZ  
- output.tif: The output raster. This will appear in the same folder as the program. Replace "output" with the name of the .tif file, for example step5.tif  
- pixel\_size: The size of each pixel in the output raster, in meters. Example: 0.5  

## step6.py

Requires the following libraries:  
- numpy  
- startinpy  
- pyinterpolate  
- rasterio  

The program takes in a .LAZ point cloud file and outputs a .tif raster file containing a rasterized version of the point cloud, where the values of the raster's pixels are interpolated using ordinary kriging. The size of the pixels are decided by the user.  

The code is run with the following command:  
```
python step6.py input.LAZ output.tif pixel_size
```

- input.LAZ: The input point cloud. This must be in the same folder as the program. Replace "input" with the name of the .LAZ file, for example step4.LAZ  
- output.tif: The output raster. This will appear in the same folder as the program. Replace "output" with the name of the .tif file, for example step6.tif  
- pixel\_size: The size of each pixel in the output raster, in meters. Example: 0.5  

## step7.py

Reqeuires the following libraries:  
- numpy  
- rasterio  

The program takes in a .tif raster file and outputs a .tif raster file with hillshade visualization on the input file.  

The code is run with the following command:  
```
python step7.py input.tif output.tif [sun_height] [sun_azimuth]
```

- input.tif: The input raster. This must be in the same folder as the program. Replace "input" with the name of the .tif file, for example step5.tif  
- output.tif: The output raster. This will appear in the same folder as the program. Replace "output" with the name of the .tif file, for example step7.tif  
- *sun\_height: Optional. The height of the sun above the horizon, in degrees. By default it is 45.*  
- *sun\_azimuth: Optional. The angle of the sun clockwise from the north, in degrees. By default it is 315 (northwest).*  

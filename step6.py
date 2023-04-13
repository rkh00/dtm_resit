# Importing necessary libraries
import argparse
import numpy as np
import startinpy
import pyinterpolate
import rasterio

def test_kriging(train_data, variogram_model, ktype, test_values, number_of_neighbors, sk_mean_value=None):
    predictions = pyinterpolate.kriging(observations=train_data,
                          theoretical_model=variogram_model,
                          points=test_values[:, :-1],
                          how=ktype,
                          no_neighbors=number_of_neighbors,
                          number_of_workers=1,
                          sk_mean=sk_mean_value)
    mse = np.mean((predictions[:, 0] - test_values[:, -1])**2)
    rmse = np.sqrt(mse)
    return rmse

def create_train_test(dataset: np.ndarray, training_set_ratio=0.3):
    """
    Function divides base dataset into a training and a test set.

    Parameters
    ----------
    dataset : np.ndarray

    training_set_ratio : float, default = 0.3

    Returns
    -------
    training_set, test_set : List[np.ndarray]
    """

    np.random.seed(101)  # To ensure that we will get the same results every time

    indexes_of_training_set = np.random.choice(range(len(dataset) - 1), int(training_set_ratio * len(dataset)), replace=False)
    training_set = dataset[indexes_of_training_set]
    validation_set = np.delete(dataset, indexes_of_training_set, 0)
    return training_set, validation_set

def main():
    
    # Parsing the input arguments with ArgumentParser
    parser = argparse.ArgumentParser(
        "LAZ ordinary kriging interpolator", description="Performs ordinary kriging on an input .LAZ point cloud file to output a GeoTIFF raster with interpolated points"
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
    dt.snap_tolerance = 0.1
    dt.read_las(f"{input_file}")
    dt_points_no_inf = np.array(dt.points)[1:]
    
    # Selecting 10000 random points from the input
    mask = np.zeros(len(dt_points_no_inf),dtype=bool)
    indices = np.random.choice(len(dt_points_no_inf), size=10000, replace=False)
    mask[indices] = True
    random_points = dt_points_no_inf[mask]

    # Splitting the random points into a training dataset and testing dataset
    train_set, test_set = create_train_test(random_points)
    
    # Prepare experimental semivariogram
    step_radius = 25  # meters
    max_range = 250  # meters
    
    # Create experimental semivariogram
    exp_semivar = pyinterpolate.build_experimental_variogram(input_array=train_set, step_size=step_radius, max_range=max_range)
    print("Experimental semivariogram data:")
    print(exp_semivar)
    exp_semivar.plot(plot_semivariance=True, plot_covariance=True, plot_variance=True)
    
    # Fit data into a theoretical model
    semivar = pyinterpolate.TheoreticalVariogram()
    fitted = semivar.autofit(experimental_variogram=exp_semivar)
    print("Theoretical semivariogram data:")
    print(fitted)
    print(semivar)
    semivar.plot()
    
    # Select one known value
    random_index = np.random.randint(0,len(train_set))
    known_value = train_set[random_index]
    print(known_value)
    
    # Predict with Ordinary Kriging
    ok_interpolation = pyinterpolate.kriging(train_set, semivar, [known_value[:-1]])
    print(ok_interpolation)
    
    # Number of neighbors
    no_of_n = [4, 8, 16, 32, 64, 128, 256]
    errors = np.zeros(len(no_of_n))
    print('Ordinary Kriging: tests')
    print('')
    
    # Testing the model using different numbers of neighbors and finding their
    # root mean squared error (RMSE)
    for i, nn in enumerate(no_of_n):
        print('Number of neighbors:', nn)
        rmse_pred = test_kriging(train_data=train_set, variogram_model=semivar, ktype='ok', test_values=test_set, number_of_neighbors=nn)
        print('RMSE:', rmse_pred)
        print('')
        errors[i] = rmse_pred
        
    # Finding the number of neighbors parameter with the smallest error
    least_rmse = np.argmin(errors)
    best_nn = no_of_n[least_rmse]

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
                point = np.array(x[k],y[j])
                z[j][k] = pyinterpolate.ordinary_kriging(theoretical_model=semivar, known_locations=train_set, unknown_location=point, no_neighbors=best_nn)
                
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
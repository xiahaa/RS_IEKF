# RS_IEKF
Parameter estimation using IEKF for RS camera. The following figure shows the illustration.

![fig](https://github.com/xiahaa/RS_IEKF/blob/master/figs/filter.png)

## Src
In `src` folder, you will find all source code.

## Dependency
To run with the `main.m`, you need to download the [dataset](http://users.ece.utexas.edu/~bevans/projects/dsc/software/calibration/) provided in this paper.
> Jia, C., & Evans, B. L. (2013, December). Online calibration and synchronization of cellphone camera and gyroscope. In 2013 IEEE Global Conference on Signal and Information Processing (pp. 731-734). IEEE.

1. download the file;
2. create a folder called `thirdparty`;
3. copy the downloaded file to `thirdparty` and uncompress it.
4. make sure the directory is like `thirdparty/gyro_camera_calibration/`

## Usage
1. open `main.m` and run.
2. you can config the filter type in the config.m.

## Note
if you use the mentioned dataset, please cite the original work
> Jia, C., & Evans, B. L. (2013, December). Online calibration and synchronization of cellphone camera and gyroscope. In 2013 IEEE Global Conference on Signal and Information Processing (pp. 731-734). IEEE.

if you use the re-implemented version `ekf_epipolar_analytic.m`, please cite the original work by the author
> Jia, C., & Evans, B. L. (2014). Online camera-gyroscope autocalibration for cell phones. IEEE Transactions on Image Processing, 23(12), 5070-5081.

if you use other files, please cite this work (format release later). 

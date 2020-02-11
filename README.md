# RS_IEKF
Parameter estimation using IEKF for RS camera. The following figure shows the illustration.

![fig](https://github.com/xiahaa/RS_IEKF/blob/master/figs/filter2.png)

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

## More results

Experiments on public datasets provided in this paper.
> Hannes Ovrén, Per-Erik Forssén, "Gyroscope-based video stabilisation with auto-calibration", 2015 IEEE INTERNATIONAL CONFERENCE ON ROBOTICS AND AUTOMATION (ICRA), IEEE International Conference on Robotics and Automation ICRA, 2090-2097, 2015.

**rccar dataset**
Estimations are:
tr: 0.0294 
td: 5.2006;
wd: [-0.0043 -6.4969e-04 -0.0013];
Ric: [0.0293    0.9996    0.0011; -0.9994    0.0293   -0.0185;-0.0185   -0.0006    0.9998]

The references are:
tr: 0.0316734
td: 5.20867366868152E+00
wd: [-5.09433180476134E-03	-7.9750612785219E-03	7.11096622185005E-03];
Ric: [  -0.0712    0.9974   -0.0130; -0.9954   -0.0719   -0.0640; -0.0648    0.0084    0.9979]

**rotation dataset**
Estimations are:
tr: 0.0315;
td: 2.8005;
wd: [2.0156e-04 1.8922e-04 -2.1981e-04];
Ric: [ 0.0307    0.9995   -0.0001; -0.9995    0.0307    0.0000; 0.0000    0.0001    1.0000];

The references are:
tr: 0.0316734
td: 2.83818847389866E+00
wd: [-4.89355970494164E-03	-7.69152770619337E-03	9.78420813209227E-03];
Ric: [ 0.0288    0.9996   -0.0011;-0.9996    0.0288   -0.0030;-0.0029    0.0012    1.0000];

**walk dataset**
Estimations are:
tr: 0.0309;
td: 3.3904;
wd: [0.0013 -0.0021 3.7599e-04];
Ric: [ 0.0310    0.9995   -0.0002; -0.9995    0.0310    0.0027; 0.0027    0.0001    1.0000];

The references are:
tr: 0.0316734
td: 3.38581144830956E+00
wd: [-3.76554195856584E-03	-1.13762523611197E-02	1.05866716566729E-02];
Ric: [ 0.0213    0.9998    0.0024; -0.9996    0.0212    0.0164; 0.0164   -0.0028    0.9999];


## Note
if you use the mentioned dataset, please cite the original work.

if you use the re-implemented version `ekf_epipolar_analytic.m`, please cite the original work by the author
> Jia, C., & Evans, B. L. (2014). Online camera-gyroscope autocalibration for cell phones. IEEE Transactions on Image Processing, 23(12), 5070-5081.

if you use the dataset in the `data` folder, please cite the original work.

if you use other files, please cite this work (format release later). 

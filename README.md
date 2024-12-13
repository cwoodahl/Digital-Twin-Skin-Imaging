Our multiple_layers.m is an edited version of the file in github repository: https://github.com/PranilJ/skin_reflectance_mcmatlab. 
This uses MCMatlab versionn 3.6.2 to simulate skin reflectance. 
Mu absorption and scattering, and geometries can be adjusted to test different skin conditions.

sensorCalibration.m uses the simulated reflectances from multiple_layers.m file to process the reflected light under 
various illuminants, and at various sensors.

plotSensorPerformance.m pulls raw pixel data and processed image data, and plots averages. 
This code is specifed for like the scenes we create in sensorCalibration.m, but can be adjusted for other formats.

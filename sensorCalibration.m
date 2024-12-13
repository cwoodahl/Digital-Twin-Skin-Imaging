%This code is designed to run an scene and illuminant through a sensor and output
%sensor values and a processed image from the sensor.

%ieInit

patchSize = 24;
wavelength = 300:4:800;
reflectances1 = load('Ep1_v4.mat');
%reflectances1 = load('Epidermis_1_Results.mat');
data = transpose(reflectances1.output(:,:,1));
%dataFileNames = which('Ep1_v3.mat');
reflectances1 = load('Ep6_v3.mat');
data = [data, transpose(reflectances1.output(:,:,1))];
data = [data,data];
%dataFileNames = which('Ep6_v3.mat');
save('Ep6.mat','data','wavelength');

save('Ep1.mat','data','wavelength');

sFiles = cell(1,1);
sFiles{1} = 'Ep1.mat';

% The number of samples from each of the data sets, respectively
sSamples = 16;     %number of columns to consider

% How many row/col spatial samples in each patch (they are square)
pSize    = patchSize;           % Patch size
wave     = wavelength;           % Whatever is in the reflectance data file
grayFlag = 0;            % No gray strip
sampling = 'no replacement';
scene1 = sceneCreate('reflectance chart',pSize,sSamples,sFiles,wave,grayFlag,sampling);

scene = scene1;

wave = sceneGet(scene,'wave');  d65 = ieReadSpectra('D65',wave);
sceneD65 = sceneAdjustIlluminant(scene,d65);
sceneD65 = sceneSet(sceneD65,'name','Reflectance Chart D65');

custom = normpdf(wave,618,20); custom = custom/max(custom)*200;
sceneCustom1 = sceneAdjustIlluminant(scene,custom);
sceneCustom1 = sceneSet(sceneCustom1,'name','Custom');

custom2 = normpdf(wave,520,20); custom2 = custom2/max(custom2)*100;
sceneCustom2 = sceneAdjustIlluminant(scene,custom2);
sceneCustom2 = sceneSet(sceneCustom2,'name','Custom2');

d50 = ieReadSpectra('d75',wave);
sceneD50 = sceneAdjustIlluminant(scene,d50);
sceneD50 = sceneSet(sceneD50,'name','D50');


scene = sceneCombine(scene,sceneD65,'direction','horizontal');
scene = sceneCombine(scene,sceneCustom1,'direction','horizontal');
scene = sceneCombine(scene,sceneCustom2,'direction','horizontal');


hfov = 20;
scene = sceneSet(scene,'fov',hfov);
vfov  = sceneGet(scene,'v fov');
sceneWindow(scene);

%%
oi = oiCreate;
oi = oiSet(oi,'optics fnumber',1.2);
oi = oiCompute(oi,scene);
oiWindow(oi);

%% Now run through some sensors

%sensorList = {'bayer-rggb','imx363','rgbw','mt9v024','imec44','cyym','monochrome'};
%pick a sensor from list, or other iSetCam/personal source
sensorList = {'cyym'};


for ii=1:numel(sensorList)
    if isequal(sensorList{ii},'mt9v024') 
        sensor = sensorCreate(sensorList{ii},[],'rccc');
    elseif isequal(sensorList{ii},'ar0132at')
        sensor = sensorCreate(sensorList{ii},[],'rccc');
    else
        sensor = sensorCreate(sensorList{ii});
    end

    sensor = sensorSet(sensor,'pixel size constant fill factor',1.4e-6);
    sensor = sensorSet(sensor,'hfov',hfov,oi);
    sensor = sensorSet(sensor,'vfov',vfov);
    sensor = sensorSet(sensor,'auto exposure',true);
    sensor = sensorCompute(sensor,oi);
    sensorWindow(sensor);

    switch sensorList{ii}
        case 'imx363'
            ip = ipCreate('imx363 RGB',sensor);
            ip = ipCompute(ip,sensor);
            ipWindow(ip);
        case 'mt9v024'
            ip = ipCreate('mt9v024 RCCC', sensor);
            % NOTE: ipCreate doesn't seem to take its cue from the 
            %       sensor that it is rccc, so we do it manually
            ip = ipSet(ip,'demosaic method','analog rccc');
            ip = ipCompute(ip,sensor);
            ipWindow(ip);
        case 'bayer-rggb'
            ip = ipCreate('bayer-rggb',sensor);
            ip = ipCompute(ip,sensor);
            ipWindow(ip);
        case 'rgbw'
            ip = ipCreate('rgbw',sensor);
            ip = ipCompute(ip,sensor);
            ipWindow(ip);
        case 'imec44'
            ip = ipCreate('imec44',sensor);
            ip = ipCompute(ip,sensor);
            ipWindow(ip);
        case 'cyym'
            ip = ipCreate('cyym',sensor);
            ip = ipCompute(ip,sensor);
            ipWindow(ip);
        case 'monochrome'
            ip = ipCreate('monochrome',sensor);
            ip = ipCompute(ip,sensor);
            ipWindow(ip);
        case 'ar0132at'
            ip = ipCreate('ar0132at RCCC', sensor);
            % NOTE: ipCreate doesn't seem to take its cue from the 
            %       sensor that it is rccc, so we do it manually
            ip = ipSet(ip,'demosaic method','analog rccc');
            ip = ipCompute(ip,sensor);
            ipWindow(ip);
    end


end

%%

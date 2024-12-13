%% Description
% In this example, we show two things: The use of cyclic boundary
% conditions and how to calculate the reflectance of a sample, including
% how to provide an estimate of the statistical error of the value.
%
% The geometry is similar to that of example 1, although we will run it
% both with and without matched interfaces here. We place a thin layer of
% air, 1 voxel thick, at the top of the cuboid and fill the rest with
% "standard tissue". The thin layer of air is necessary to enforce specular
% reflectance of the incident beam and to include the reflection and
% refraction effects at the surface.
%
% We want to calculate the total and diffuse reflectance of the tissue, so
% we don't want photons to escape at the side walls and disappear from the
% simulation. One possibility for how to avoid this is to use a cuboid with
% very large Lx and Ly. Another is to use cyclic boundary conditions, which
% we emply here by setting model.MC.boundaryType = 3. With cyclic boundary
% conditions, all photons that hit a side wall will immediately enter the
% cuboid again on the opposite wall. Because our geometry is supposed to
% represent a slab of tissue with infinite horizontal extent, this is a
% valid way for us to avoid losing photons to side wall effects.
%
% We will calculate the reflectance by integrating the
% model.MC.normalizedIrradiance_zneg 2D array, which contains all the power
% that has hit the top cuboid boundary from the inside. In principle, the x
% and y resolution can be set as low as the minimum value, 2x2 and still
% get the same result (try it). For visualization purposes, however, we
% keep nx and ny at reasonable values of 101 in this example. Also Lx and
% Ly could be set arbitrarily low due to the use of the cyclic boundaryType.
%
% The reflectance that we calculate is the total refletance, including the
% specular reflection that we get when simulating without matched
% interfaces. We can obtain the diffuse reflectance by subtracting the
% Fresnel reflectivity, ((1 - 1.4)/(1 + 1.4))^2, from the total
% reflectance. Another way of avoiding the specular reflectance would
% have been to set depositionCriteria.minInterfaceTransitions = 1, but we
% do not use that method here.
%
% To get some statistics, we launch 1e6 photons 5 times and collect the
% diffuse reflectance value for each run. Then we finally write out the
% mean and the standard error of the mean for both the calculation with
% matched interfaces and without.

%load spectral_library_1.mat;

%% Global mu_s values
global mus_air;
mus_air_ = 100; % [cm^-1]
global mus_st;
mus_st_ = 1.8; % [cm^-1]
global mus_ep1;
mus_ep1_ = 127.47; % [cm^-1]
global mus_ep6;
mus_ep6_ = 127.47; % [cm^-1]
global mus_ud;
mus_ud_ = 127.54; % [cm^-1]
global mus_blood;
mus_blood_ = 31.87; % [cm^-1]
global mus_ld;
mus_ld_ = 124.54; % [cm^-1]
global mus_fat;
mus_fat_ = 109.78; % [cm^-1]
global mus_black;
mus_black_ = 1; % [cm^-1]

global mua_ep1;
mua_ep1_ = 1.865; % [cm^-1]
global mua_blood;
mua_blood_ = 1; % [cm^-1]

load('C:\Users\achip\Documents\GitHub\isetcam\data\surfaces\reflectances\skin\absorbances\SkinComponentAbsorbances.mat')
MU = data;
wavelengthsData = wavelength;

%% Ouptput Matrix
wavelengths = (300:12:800);
bo_arr = (0.85:0.05:1);

output = zeros(size(bo_arr, 2), size(wavelengths, 2), 2);
%% 



i = 0;
j = 0;
for bo_val = bo_arr
    j = j + 1;
    i = 0;
    for wave = wavelengths
    wave
        %% Geometry definition
        model = MCmatlab.model;
        model.G.nx                = 101; % Number of bins in the x direction
        model.G.ny                = 101; % Number of bins in the y direction
        model.G.nz                = 500; % Number of bins in the z direction
        % Least counts: X = 50um, Y = 50um, Z = 10um 
        model.G.Lx                = .5; % [cm] x size of simulation cuboid
        model.G.Ly                = .5; % [cm] y size of simulation cuboid
        model.G.Lz                = .5; % [cm] z size of simulation cuboid
        
        model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
        model.G.geomFunc          = @geometryDefinition; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file
        model = plot(model,'G');

        i = i + 1;
        %% Monte Carlo simulation
        model.MC.nPhotonsRequested        = 1e5; % Number of photons requested
        model.MC.boundaryType             = 3; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping
        model.MC.wavelength               = wave; % [nm] Excitation wavelength, used for determination of optical properties for excitation ligh
        model.MC.lightSource.sourceType          = 2; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
        model.MC.lightSource.xFocus              = 0; % [cm] x position of focus
        model.MC.lightSource.yFocus              = 0; % [cm] y position of focus
        model.MC.lightSource.zFocus              = 0; % [cm] z position of focus

        %% Update Mus based on wavelength 
        mus_air = mus_air_ * ((523 / model.MC.wavelength)^0.6);
        mus_st = mus_st_ * ((523 / model.MC.wavelength)^0.6); % [cm^-1]
        mus_ep1 = mus_ep1_ * ((523 / model.MC.wavelength)^0.6); % [cm^-1]
        mus_ep6 = mus_ep6_ * ((523 / model.MC.wavelength)^0.6);
        mus_ud = mus_ud_ * ((523 / model.MC.wavelength)^0.6);
        mus_blood = mus_blood_ * ((523 / model.MC.wavelength)^0.6); % [cm^-1]
        mus_ld = mus_ld_ * ((523 / model.MC.wavelength)^0.6); % [cm^-1]
        mus_fat = mus_fat_ * ((523 / model.MC.wavelength)^0.6);
        mus_black = mus_black_ * ((523 / model.MC.wavelength)^0.6);

        mua_ep1 = get_mua_epidermis(wave, MU,wavelengthsData);

        %bo_level = 0.9;
        % model.MC.wavelength
        [mua_h, mua_d] = get_mua_from_blood(model.MC.wavelength, MU,wavelengthsData);
        mua_blood = bo_val * mua_h + (1 - bo_val) * mua_d;

        % Set up arrays to collect statistics in:
        n = 5; % Times to run
        Rd_matched = NaN(1,n);
        Rd_mismatched = NaN(1,n);

        % First run with matched interfaces:
        model.MC.matchedInterfaces = true; % Assumes all refractive indices are the same
        tic
        for iRun = 1:n
            model = runMonteCarlo(model);
            Rd_matched(iRun) = sum(model.MC.normalizedIrradiance_zneg(:))*model.G.dx*model.G.dy; % Diffuse reflectance
            %Rd_matched_(iRun) = sum(model.MC.normalizedIrradiance_zpos(:))*model.G.dx*model.G.dy; % Diffuse reflectance
        end
        t_Matched = toc;

        model = plot(model,'MC');

        %% Print outputs:
        output(j, i, 1) = mean(Rd_matched);
        output(j, i, 2) = std(Rd_matched)/sqrt(n);    
        fprintf('\nReflectance    = %.6f +- %.6f (n = %d, total time elapsed = %d s)\n'   ,mean(Rd_matched   ),std(Rd_matched   )/sqrt(n),n,round(t_Matched   ));
        % fprintf('Rd mismatched = %.6f +- %.6f (n = %d, total time elapsed = %d s)\n\n',mean(Rd_mismatched),std(Rd_mismatched)/sqrt(n),n,round(t_Mismatched));

        %fprintf('\nRd matched    = %.6f +- %.6f (n = %d, total time elapsed = %d s)\n'   ,mean(Rd_matched_   ),std(Rd_matched_   )/sqrt(n),n,round(t_Matched   ));
        %fprintf('Rd mismatched = %.6f +- %.6f (n = %d, total time elapsed = %d s)\n\n',mean(Rd_mismatched_),std(Rd_mismatched_)/sqrt(n),n,round(t_Mismatched));
        %% Geometry function(s)
        % A geometry function takes as input X,Y,Z matrices as returned by the
        % "ndgrid" MATLAB function as well as any parameters the user may have
        % provided in the definition of Ginput. It returns the media matrix M,
        % containing numerical values indicating the media type (as defined in
        % mediaPropertiesFunc) at each voxel location.
    end
end
output

function M = geometryDefinition(X,Y,Z,parameters)
    M = 1*ones(size(X)); % "Standard" tissue
    M(:,:,:) = 9; % Make everything black first
    M(:,:,1) = 1; % air
    M(:,:,2:11) = 3; % epidermis
    M(:,:,12:15) = 4; % UD1
    M(:,:,16:25) = 5; % blood 1
    M(:,:,26:114) = 4; % UD2
    M(:,:,115:190) = 6; % LD
    for x = 1:5:101     % Blood 2, at the center of LD layer, each vessel is 2x2, and spaced 3 boxes apart
        M(x,:,152) = 5;
        M(x,:,153) = 5;
        sizeM = size(M);
        if x < sizeM(1)
            M(x+1,:,152) = 5;
            M(x+1,:,153) = 5;
        end
    end
end

%% Media Properties function
% The media properties function defines all the optical and thermal
% properties of the media involved by constructing and returning a
% "mediaProperties" struct with various fields. As its input, the function
% takes the wavelength as well as any other parameters you might specify
% above in the model file, for example parameters that you might loop over
% in a for loop. Dependence on excitation fluence rate FR, temperature T or
% fractional heat damage FD can be specified as in examples 12-15.
function mediaProperties = mediaPropertiesFunc(wavelength, parameters)
%   mediaProperties = MCmatlab.mediumProperties;
  global mus_air;
  global mus_st;
  global mus_ep1;
  global mus_ep6;
  global mus_blood;
  global mus_ud;
  global mus_ld;
  global mus_black;
  global mus_fat;

  global mua_blood;
  global mua_ep1;
  
  j=1;
  mediaProperties(j).name  = 'air';
  mediaProperties(j).mua   = 1e-8; % [cm^-1]
  mediaProperties(j).mus   = 1e-8; %mus_air; % [cm^-1]
  mediaProperties(j).g     = 1;
  mediaProperties(j).n     = 1;

  j=2;
  mediaProperties(j).name  = 'standard tissue';
  mediaProperties(j).mua   = 0.2; % [cm^-1]
  mediaProperties(j).mus   = mus_st; % [cm^-1]
  mediaProperties(j).g     = 0.9; % 0.99; % 1; % 0.90;
  mediaProperties(j).n     = 1.4;

  j=3;
  mediaProperties(j).name  = 'Epidermis1';
  mediaProperties(j).mua   = mua_ep1; % [cm^-1]
  mediaProperties(j).mus   = mus_ep1; % [cm^-1]
  mediaProperties(j).g     = 0.7; % 0.99; % 1; % 0.90;
  mediaProperties(j).n     = 1.47;

  j=8;
  mediaProperties(j).name  = 'Epidermis6';
  mediaProperties(j).mua   = 24.668; % [cm^-1]
  mediaProperties(j).mus   = mus_ep6; % [cm^-1]
  mediaProperties(j).g     = 0.9; % 0.99; % 1; % 0.90;
  mediaProperties(j).n     = 1.47;

  j=4;
  mediaProperties(j).name  = 'Upperdermis';
  mediaProperties(j).mua   = 1.949; % [cm^-1]
  mediaProperties(j).mus   = mus_ud; % [cm^-1]
  mediaProperties(j).g     = 0.7; % 0.99; % 1; % 0.90;
  mediaProperties(j).n     = 1.47;

  j=5;
  mediaProperties(j).name  = 'Blood';
  mediaProperties(j).mua   = mua_blood*0.8; % [cm^-1]
  mediaProperties(j).mus   = mus_blood; % [cm^-1]
  mediaProperties(j).g     = 0.9; % 0.99; % 1; % 0.90;
  mediaProperties(j).n     = 1.47;

  j=6;
  mediaProperties(j).name  = 'Lowerdermis';
  mediaProperties(j).mua   = 17.881; % [cm^-1]
  mediaProperties(j).mus   = mus_ld; % [cm^-1]
  mediaProperties(j).g     = 0.7; % 0.99; % 1; % 0.90;
  mediaProperties(j).n     = 1.47;

  j=7;
  mediaProperties(j).name  = 'Fat';
  mediaProperties(j).mua   = 0.0001; % [cm^-1]
  mediaProperties(j).mus   = mus_fat; % [cm^-1]
  mediaProperties(j).g     = 0.7; % 0.99; % 1; % 0.90;
  mediaProperties(j).n     = 1.47;

  j=9;
  mediaProperties(j).name  = 'Black';
  mediaProperties(j).mua   = 1000000; % [cm^-1]
  mediaProperties(j).mus   = mus_black; % [cm^-1]
  mediaProperties(j).g     = 0.7; % 0.99; % 1; % 0.90;
  mediaProperties(j).n     = 1.47;
end

function [mua_h, mua_d] = get_mua_from_blood(wavelength, MU,wavelengthsData)
    
    [minVal, index] = min(abs(wavelengthsData - wavelength));
    %index = wavelength - 300;

    hemoglobin_absorption_spectrum = MU(:,1);
    deoxyhemoglobin_absoprtion_spectrum = MU(:,2);

    mua_h = hemoglobin_absorption_spectrum(index);
    mua_d = deoxyhemoglobin_absoprtion_spectrum(index);

end

%{
function [mu_a] = get_mua_epidermis(wavelength, MU)
    arr = [2.305555556, 1, 0.7847222222, 0.6527777778, 0.5, 0.3402777778, 0.2708333333, 0.1180555556, 0.09027777778];
    wave_arr = [405, 532, 595, 632, 694, 755, 800, 980, 1064];
    
    mu_a = -1;
    for i = 1:length(arr)-1
        if wave_arr(i) == wavelength
            mu_a = MU * (arr(i));
            break
        elseif wave_arr(i) < wavelength && wave_arr(i+1) > wavelength
            mu_a = MU * (arr(i) + ((arr(i+1) - arr(i))*(wavelength - wave_arr(i)) / (wave_arr(i+1) - wave_arr(i))));
            break
        end
    end
    assert(mu_a ~= -1, "Hagaa");
end
%}

function [mu_a] = get_mua_epidermis(wavelength,MU,wavelengthsData)
    [minVal, index] = min(abs(wavelengthsData - wavelength));
    
    wavelengths = wavelengthsData;

    %skinBaseLine = (7.84*10^8).*(wavelengths.^(-3.255));
    skinBaseLine = 0.244 + 85.3*exp(-(wavelengths-154)/64);
    mel = (6.6*10^11)*(wavelengths.^(-3.33));
    
    freqmel = 0.02;
    epi1 = freqmel*mel+ (1-freqmel)*skinBaseLine;
    
    freqmel = 0.27;
    epi6 = freqmel*mel+ (1-freqmel)*skinBaseLine;
    
    mu_a = epi1(index)*0.5;
end
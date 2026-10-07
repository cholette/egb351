clearvars, close all
clc
addpath('../../../matlab_library/sun_position')

%% inputs
fileName = 'brisbane_tmy.csv';
lat = -27.470; 
lon = 360-153.030; % our convention for lon
lon_tz = 360 - 150.0;

%% Import data
opts = detectImportOptions(fileName);
opts.DataLines = [4 Inf];
opts.VariableNamesLine  = 3;
opts.SelectedVariableNames = opts.VariableNames(1:20);
T = readtable(fileName,opts);

year = 2026;
T.datetime = datetime(year,T.Month,T.Day,T.Hour,T.Minute,0,'TimeZone','+10:00');
T.doy = day(T.datetime,'dayofyear');

%% Compute zenith and azimuth angle of sun
tod = timeofday(T.datetime);
T.SolarTime = tod + minutes(4*(lon_tz-lon) + equation_of_time(T.doy));
T.HourAngle = 15*(hours(T.SolarTime)-12);

[z,az]= solar_angles(T.doy,T.HourAngle,lat);
T.SolarZenithAngle = z;
T.SolarAzimuthAngle = az;

figure(name="Solar angles")
tiledlayout(2,1)
ax1 = nexttile;
plot(T.datetime,z)
title('Solar zenith (deg)')


ax2 = nexttile;
plot(T.datetime,az)
title('Solar azimuth (deg, south zero')

% [fig,ax] = sun_path_diagram(T)

%% b) Compute the DHI & compare with TMY values
T.DHI2 = T.GHI - T.DNI.*cosd(T.SolarZenithAngle);


figure
tiledlayout(2,1)
nexttile
plot(T.datetime,T.DHI,T.datetime,T.DHI2,'r--')

nexttile
scatter(T.DHI,T.DHI2)
xlabel('TMY DHI')
ylabel('My DHI')
hold on
plot([min(T.DHI),max(T.DHI)],[min(T.DHI),max(T.DHI)],'g')

%% part c & d POA irradiance & annual insolation (flat panel)
beta = -lat;
gamma = 180;
ground_albedo = 0;
[GPOA,~] = plane_of_array_irradiance(beta,gamma,ground_albedo,T.GHI,...
    T.DNI,T.SolarZenithAngle,T.SolarAzimuthAngle);


figure
plot(T.datetime,GPOA)
title(sprintf('Plane of array Irradiance (beta=%.1f, gamma=%.2f)',...
    beta,gamma))

dt = 1; % time step is one hour
insolation = sum(GPOA*dt,'omitnan')/1000

%% Change the tilt
ground_albedo = 0;
obj = @(x) -sum(plane_of_array_irradiance(x,gamma,ground_albedo,...
    T.GHI,T.DNI,T.SolarZenithAngle,T.SolarAzimuthAngle));
[xs,fval] = fminbnd(obj,0,90);

beta = xs
POAI = -fval;

%% Question 2
beam_shading = readtable('shading_2ai.csv','NumHeaderLines',4);
fb = 100-beam_shading(:,end).Variables;
fd = 100-15.47;

[~,AOI] = plane_of_array_irradiance(0,138,0,T.GHI,T.DNI,...
    T.SolarZenithAngle,T.SolarAzimuthAngle);

GPOA_2ai = fb/100 .* T.DNI .* cosd(AOI) ... 
    + fd/100 .* T.DHI .* 0.5*(1+cosd(0));

insolation_2ai = sum(GPOA_2ai)/1000


beam_shading = readtable('shading_2aii.csv','NumHeaderLines',4);
fb = 100-beam_shading(:,end).Variables;
fd = 100-15.47;

[~,AOI] = plane_of_array_irradiance(beta,0,0,T.GHI,T.DNI,...
    T.SolarZenithAngle,T.SolarAzimuthAngle);

GPOA_2aii = fb/100 .* T.DNI .* cosd(AOI) ... 
    + fd/100 .* T.DHI .* 0.5*(1+cosd(beta));

insolation_2aii = sum(GPOA_2aii)/1000

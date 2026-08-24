clearvars, close all

% Add path to the MATLAB library for EGB351 from
% https://github.com/cholette/egb351. This will
% depend on where you save it.
addpath('../../../code/matlab_library/') 

% inputs
offset = 10;        % UTC hour offset for Brisbane (thankfully constant!)
lat = -27.470;
lon = 360-153.030;  % Convention is degrees west of 0
lon_tz = 360-150;   % Brisbane local timzone meridian (careful here!!)

%% Add solar angles to the MATLAB Table (part a)

% Read in data from the TMY File
opts = detectImportOptions("brisbane_tmy.csv");
opts.DataLines = [4 Inf]; % for TMY, data starts on line 4
opts.VariableNamesLine = 3;
opts.SelectedVariableNames = opts.VariableNames(1:20); 
T = readtable("brisbane_tmy.csv",opts);

% Add datetime, DOY to table
year = 2025; % don't use year in TMY (it's the year that was "typical")
T.datetime = datetime(year,T.Month,T.Day,T.Hour,T.Minute,0,'TimeZone','+10:00');
T.doy = day(T.datetime,'dayofyear');

% Compute solar time (be careful with datetimes!!!!!)
tod = timeofday(T.datetime);
T.SolarTime = tod + minutes(4*(lon_tz-lon) + ...
    equation_of_time(T.doy)); % hour of the day, note that adjustments have units of minutes
T.HourAngle = 15*(hours(T.SolarTime)-12.0);

% Zenith and azimuth angle
[z,az] = solar_angles(T.doy,T.HourAngle,lat);
T.SolarZenithAngle = z;
T.SolarAzimuthAngle = az;

figure(name='Solar angles')
tiledlayout(2,1)
ax1 = nexttile;
plot(T.datetime,z);
title('Solar zenith (deg)')

ax2 = nexttile;
plot(T.datetime,az);
title('Solar azimuth (deg, South zero)')

linkaxes([ax1,ax2],'x')

% sun paths
[fig,ax] = sun_path_diagrams(T); % no need to do this in class (show result?)

%% Compute DHI & Compare with TMY Values (Part b)
T.DHI2 = T.GHI - T.DNI.*cosd(T.SolarZenithAngle); % On equation sheet
figure
tiledlayout(2,1)
nexttile
plot(T.datetime,T.DHI,T.datetime,T.DHI2,'r--')
ylabel('DHI')
legend("TMY File","Computations")
nexttile
scatter(T.DHI,T.DHI2)
hold on
plot([min(T.DHI),max(T.DHI)],[min(T.DHI),max(T.DHI)],'k',linewidth=3)
xlabel('DHI TMY')
ylabel('DHI Computed')
legend('Data','Perfect match')
hold off

%% Part c & d: POA Irradiance & annual insolation (flat panel)
beta = 0;
gamma = 180; % South is zero
ground_albedo = 0.2;
[GPOA,~] = plane_of_array_irradiance(beta,gamma,ground_albedo,...
    T.GHI,T.DNI,z,az);

figure
plot(T.datetime,GPOA)
title(sprintf('Plane of array irradiance (beta=%.1f,gamma=%.2f)',beta,gamma))

insolation = sum(GPOA,'omitnan')/1e3; % in kWh/m2 because time step is 1 hour
fprintf('Annual Insolation (no shade): %.2f kWh/m2\n',insolation)

%% Maximize insolation (Part e & f)
ground_albedo = 0;
ng = 73;
panel_gamma = linspace(0,360,ng); %azimuth angle
beta = nan(ng,1);
POAI = nan(ng,1);
for g = panel_gamma
    obj = @(x) -sum(plane_of_array_irradiance(x,g,ground_albedo,...
        T.GHI,T.DNI,z,az))/1e3; % this is annual insolation in kWh since timestep is one hour
    [xs,fval] = fminbnd(obj,0,90);
    beta(g==panel_gamma) = xs;
    POAI(g==panel_gamma) = -fval;
end

figure
yyaxis left
plot(panel_gamma,POAI)
xlabel('Panel azimuth (south zero)')
ylabel('Annual insolation kWh/m^2')
yyaxis right
plot(panel_gamma,beta)
ylabel('Optimal panel tilt')

best_idx = find(POAI == max(POAI));
best_gamma = panel_gamma(best_idx);
best_beta = beta(best_idx);

GPOA = plane_of_array_irradiance(best_beta,best_gamma,ground_albedo,...
    T.GHI,T.DNI,z,az);
fprintf('Optimal insolation without shading(az=%.1f,tilt=%.1f): %.1f kWh/m2\n',best_gamma,best_beta,sum(GPOA)/1e3)

% re-run with ground_albedo == 0

%% Question 2a, case i (flat panel)
beam_shading = readtable("shading_2ai.csv",'NumHeaderLines',4);
fb = 100-beam_shading(:,end).Variables; % [%]
fd = 100-15.47; % from SAM [%]
[~,AOI] = plane_of_array_irradiance(0,138,0,T.GHI,T.DNI,z,az);

GPOA_2ai = fb/100 .* T.DNI .*cosd(AOI) + fd/100 .* T.DHI .* 0.5*(1+cosd(0));
fprintf('Insolation with shading(az=%.1f,tilt=%.1f): %.1f kWh/m2\n',138,0,sum(GPOA_2ai)/1e3)

%% Question 2a, case ii (optimal orientation)
beam_shading = readtable("shading_2aii.csv",'NumHeaderLines',4);
fb = 100-beam_shading(:,end).Variables; % [%]
fd = 100-15.85; % from SAM [%]
[~,AOI] = plane_of_array_irradiance(best_beta,best_gamma,0,T.GHI,T.DNI,z,az);

GPOA_2aii = fb/100 .* T.DNI .*cosd(AOI) + fd/100 .* T.DHI .* 0.5*(1+cosd(0));
fprintf('Insolation with shading(az=%.1f,tilt=%.1f): %.1f kWh/m2\n',best_gamma,best_beta,sum(GPOA_2aii)/1e3)
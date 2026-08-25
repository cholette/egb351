clearvars, close all

%%%%%%%%%%%%%%%%% Exercise 1 %%%%%%%%%%%%%%%%%%%%%%

%%(b)
wind_data = import_nasa_power_wind(...
    'POWER_Point_Hourly_20150101_20241231_026d73S_151d47E_UTC.csv');
wind_data = localize_time(wind_data,10)

%% (c)
pax = wind_rose(wind_data.WD50M);

%%(d)--(e)
a=1/7.0;
wind_data.WS100M = wind_data.WS50M * (100/50)^a;
wind_data.P100M = wind_data.PS 

%% (e)
v_bar_100 = mean(wind_data.WS100M);
fprintf('Mean: %.2f, Std dev: %.2f\n',v_bar_100,std(wind_data.WS100M))

%% (f)--(g)
ws_grid_100 = linspace(min(wind_data.WS100M),max(wind_data.WS100M),1000);
dist = fitdist(wind_data.WS100M,"Weibull");

figure
histogram(wind_data.WS100M,"DisplayName","Data",'Normalization','pdf')
hold on
plot(ws_grid_100,dist.pdf(ws_grid_100),"DisplayName","Fit")
legend()
xlabel('Wind speed @ 100m (m/s)')
ylabel('Density')
hold off

%% (h) WPD
rho_air = 1.20; % kg/m3
c = dist.a; % scale parameter
Ke = gamma(1+3/dist.b)/(gamma(1+1/dist.b)^3);
WPD = 0.5*rho_air*Ke*v_bar_100^3;
WPD_analytical = 0.5*rho_air*gamma(1+3/dist.b)*c^3;
fprintf('WPD (using average): %.2f kW/m^2 \n',WPD)
fprintf('WPD (analytical): %.2f kW/m^2 \n',WPD_analytical)

%% (i) velocity-duration curve
sorted_ws = sort(wind_data.WS100M);
hours_above = (length(sorted_ws)-1):-1:0; 
figure()
plot(hours_above*8760/length(sorted_ws),sorted_ws,'LineWidth',2)
ax = gca;
xlabel('Number of hours per year')
ylabel('Wind speed')

%%%%%%%%%%%%%%%%%%%%%%%%%% Exercise 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% (a) -- (c)
% locate SAM turbine definition file and paste strings in here
parse_sam =@(x) double(regexp(x,'\|','split'));

ge.speeds = parse_sam("1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|23|24|25|26|27|28|29|30");
ge.power_kw = parse_sam("0|0|0|38|154|405|693|986|1323|1506|1578|1600|1600|1600|1600|1600|1600|1600|1600|1600|1600|1600|1600|1600|1600|0|0|0|0|0");
ge.hub_height = 80;
wind_data.WS80M = wind_data.WS50M * (ge.hub_height/50)^a;
ws_grid_80 = linspace(min(wind_data.WS80M),max(wind_data.WS80M),1000);
dist_ge = fitdist(wind_data.WS80M,"Weibull");

siemens.speeds = parse_sam("1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|23|24|25|26|27|28|29|30");
siemens.power_kw = parse_sam("0|0|0|127|285|512|795|1140|1587|2120|2770|3287|3500|3578|3600|3600|3600|3600|3600|3600|3600|3600|3600|3600|3600|0|0|0|0|0");
siemens.hub_height = 90;
wind_data.WS90M = wind_data.WS50M * (siemens.hub_height/50)^a;
ws_grid_90 = linspace(min(wind_data.WS90M),max(wind_data.WS90M),1000);
dist_siemens = fitdist(wind_data.WS90M,"Weibull");

figure('Name','GE Little One')
yyaxis left
plot(ge.speeds,ge.power_kw,'LineWidth',2,DisplayName="Power Curve"); 
ylim([0,4000])
ylabel('Power (kW)')
yyaxis right
histogram(wind_data.WS100M,Normalization="pdf",DisplayName="Data")
hold on
plot(ws_grid_100,dist_ge.pdf(ws_grid_100),LineWidth=1,DisplayName='Fit')
hold off
xlabel('Wind speed (m/s)')
legend()

figure('Name',"Siemens Big One")
yyaxis left
plot(siemens.speeds,siemens.power_kw,LineWidth=2,DisplayName="Power Curve"); 
ylim([0,4000])
ylabel('Power (kW)')
yyaxis right
histogram(wind_data.WS90M,Normalization="pdf",DisplayName="Data")
hold on
plot(ws_grid_90,dist_siemens.pdf(ws_grid_90),LineWidth=1,DisplayName='Fit')
hold off
ylabel('Density')
xlabel('Wind speed (m/s)')

%% (d) aep and power duration curves
aep_ge = compute_aep(ge,dist_ge);
aep_siemens = compute_aep(siemens,dist_siemens);
fprintf('GE AEP: %.2f GWh \n',aep_ge/1e6)
fprintf('Siemens AEP: %.2f GWh \n',aep_siemens/1e6)


%% (e) Turbine Power duration curves
figure("Name","Power duration curves")
% [power_ge,hours_exceeded_ge] = power_duration(ge,dist_ge); % alternate
powerfun_ge = griddedInterpolant(ge.speeds,ge.power_kw);
power_ge = powerfun_ge(ws_grid_80);
hours_exceeded_ge = 8760*(1-dist_ge.cdf(ws_grid_80));
plot(hours_exceeded_ge,power_ge,'b-',LineWidth=2.0,DisplayName='GE (Weibull)')

% [power_siemens,hours_exceeded_siemens] = power_duration(siemens,dist_siemens); % alternate
powerfun_siemens = griddedInterpolant(siemens.speeds,siemens.power_kw);
power_siemens = powerfun_siemens(ws_grid_90);
hours_exceeded_siemens = 8760*(1-dist_siemens.cdf(ws_grid_90));
hold on
plot(hours_exceeded_siemens,power_siemens,'r-',LineWidth=2.0,DisplayName='Siemens (Weibull)')
xlabel("Hours per year")
ylabel("Power (kW)")

%% 2023 data only
power_ge = griddedInterpolant(ge.speeds,ge.power_kw);
power_siemens = griddedInterpolant(siemens.speeds,siemens.power_kw);
wind_data_2023 = wind_data(year(wind_data.timestamp)==2023,:);

% use same technique for wind speed duration curve
ep_2023_ge = sum(power_ge(wind_data_2023.WS80M));
sorted_power_ge = sort(power_ge(wind_data_2023.WS80M));
hours_ge = (length(sorted_power_ge)-1):-1:0; 

ep_2023_siemens = sum(power_siemens(wind_data_2023.WS90M));
sorted_power_siemens = sort(power_siemens(wind_data_2023.WS90M));
hours_siemens = (length(sorted_power_siemens)-1):-1:0; 

fprintf('GE 2023 Production: %.2f GWh \n',ep_2023_ge/1e6)
fprintf('Siemens 2023 Production: %.2f GWh \n',ep_2023_siemens/1e6)

% Old, slow method for power curve
% hours_above = (length(sorted_ws)-1):-1:0; 
% speeds_ge = linspace(0,max(ge.speeds),1000);
% power_2023_ge = power_ge(wind_data_2023.WS80M);
% powers_ge = unique(power_2023_ge);
% hours_ge = zeros(size(powers_ge));
% for ii = 1:length(powers_ge)
%     hours_ge(ii) = sum(power_2023_ge>powers_ge(ii)); 
% end
% 
% speeds_siemens = linspace(0,max(siemens.speeds),1000);
% power_2023_siemens = power_siemens(wind_data_2023.WS90M);
% powers_siemens = unique(power_2023_siemens);
% hours_siemens = zeros(size(powers_siemens));
% for ii = 1:length(powers_siemens)
%     hours_siemens(ii) = sum(power_2023_siemens>powers_siemens(ii)); 
% end

plot(hours_ge,sorted_power_ge,'b--',LineWidth=2,DisplayName="GE 2023 Only")
hold on
plot(hours_siemens,sorted_power_siemens,'r--',LineWidth=2,DisplayName="Siemens 2023 Only")
legend()
hold off
title("Power Duration Curves")

%% Adding direction columns and exporting 2023 data
wind_data.WD80M = wind_data.WD50M; % assuming direction is the same
wind_data.WD90M = wind_data.WD50M;
sam_export(wind_data,[datetime(2023,1,1,0,0,0),...
    datetime(2023,12,31,23,59,59)],-26.73,151.47,10)

%% OPTIONAL: check density correction in SAM for GE Turbine
% SAM does an air density correction (sec 7.3 of [1]). 
%  [1]  J. Freeman, J. Jorgenson, P. Gilman, and T. Ferguson, 
%       "Reference Manual for the System Advisor Model's Wind Power Performance 
%       Model," National Renewable Energy Laboratory (NREL), Golden, CO., 
%       NREL/TP-6A20-60570, Aug. 2014. doi: 10.2172/1150800.

pressure_pa = wind_data_2023.PS; % Pa
temp = wind_data_2023.T2M + 273.15; % K
R_air = 287.058; % J/kg/K
air_density = pressure_pa./R_air./temp; % in Pa
ep_2023_ge_SAM = sum(power_ge(wind_data_2023.WS80M).*air_density/1.225);
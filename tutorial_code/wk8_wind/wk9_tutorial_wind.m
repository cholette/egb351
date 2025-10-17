clearvars, close all

%% Exercise 1
wind_data = import_nasa_power_wind(...
    'POWER_Point_Hourly_20150101_20241231_026d73S_151d47E_UTC.csv');
wind_data = localize_time(wind_data,10)

a=0.14;
pax = wind_rose(wind_data.WD50M);
wind_data.WS100M = wind_data.WS50M * (100/50)^a;
v_bar_100 = mean(wind_data.WS100M);

ws_grid_100 = linspace(min(wind_data.WS100M),max(wind_data.WS100M),1000);

fprintf('Mean: %.2f, Std dev: %.2f\n',v_bar_100,std(wind_data.WS100M))
dist = fitdist(wind_data.WS100M,"Weibull");

figure
histogram(wind_data.WS100M,"DisplayName","Data",'Normalization','pdf')
hold on
plot(ws_grid_100,dist.pdf(ws_grid_100),"DisplayName","Fit")
legend()
xlabel('Wind speed @ 100m (m/s)')
ylabel('Density')
hold off

rho_air = 1.20; % kg/m3
WPD = 0.5*rho_air*v_bar_100^3;
fprintf('WPD: %.2f kW/m^2 \n',WPD)

% velocity-duration curve
% [fig,ax] = velocity_duration(ws); % brute force way for case with lots of
% duplicates
sorted_ws = sort(wind_data.WS100M);
hours_above = (length(sorted_ws)-1):-1:0; % smarter way from Chiara when speeds are unique
figure()
plot(hours_above*8760/length(sorted_ws),sorted_ws,'LineWidth',2)
ax = gca;
xlabel('Number of hours per year')
ylabel('Wind speed')

%% Exercise 2
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

% aep and power duration curves
aep_ge = compute_aep(ge,dist_ge);
aep_siemens = compute_aep(siemens,dist_siemens);

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

% 2023 data only
power_ge = griddedInterpolant(ge.speeds,ge.power_kw);
power_siemens = griddedInterpolant(siemens.speeds,siemens.power_kw);
wind_data_2023 = wind_data(year(wind_data.timestamp)==2024,:);
ep_2023_ge = sum(power_ge(wind_data_2023.WS80M));
ep_2023_siemens = sum(power_siemens(wind_data_2023.WS90M));

speeds_ge = linspace(0,max(ge.speeds),1000);
power_2023_ge = power_ge(wind_data_2023.WS80M);
powers_ge = unique(ge.power_kw);
hours_ge = zeros(size(powers_ge));
for ii = 1:length(powers_ge)
    hours_ge(ii) = sum(power_2023_ge>powers_ge(ii)); 
end

speeds_siemens = linspace(0,max(siemens.speeds),1000);
power_2023_siemens = power_siemens(wind_data_2023.WS90M);
powers_siemens = unique(siemens.power_kw);
hours_siemens = zeros(size(powers_siemens));
for ii = 1:length(powers_siemens)
    hours_siemens(ii) = sum(power_2023_siemens>powers_siemens(ii)); 
end

plot(hours_ge,powers_ge,'b--',LineWidth=2,DisplayName="GE 2023 Only")
plot(hours_siemens,powers_siemens,'r--',LineWidth=2,DisplayName="Siemens 2023 Only")
legend()
hold off
title("Power Duration Curves")

% Adding direction columns and exporting 2023 data
wind_data.WD80M = wind_data.WD50M; % assuming direction is the same
wind_data.WD90M = wind_data.WD50M;
sam_export(wind_data,[datetime(2023,1,1,0,0,0),datetime(2023,12,31,23,59,59)],-27.73,151.47,10)

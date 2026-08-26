clearvars, close all

% Exercise 1
%% (a+b)
file_name = "POWER_Point_Hourly_20150101_20251231_026d73S_151d47E_UTC.csv";
wind_data = import_nasa_power_wind(file_name);
wind_data = localize_time(wind_data,10);

%% (c) wind rose
pax = wind_rose(wind_data.WD50M)

%% (d) get wind at 100m
alpha = 1/7.0;
wind_data.WS100M = wind_data.WS50M * (100/50)^alpha;

%% (e)
v_bar_100 = mean(wind_data.WS100M)
v_std_100 = std(wind_data.WS100M)

%% (f)--(g)
dist = fitdist(wind_data.WS100M,"Weibull") 

figure
ws_grid = linspace(0,max(wind_data.WS100M),1000);
histogram(wind_data.WS100M,'Normalization','pdf')
hold on
plot(ws_grid,dist.pdf(ws_grid))

%% (h) WPD
rho_air = 1.22; % kg/m3
c = dist.a; 
k = dist.b;
Ke = gamma(1+3/k)/(gamma(1+1/k)^3)

WPD_average = 0.5*rho_air*Ke*dist.mean^3

%% (i) 
sorted_ws = sort(wind_data.WS100M);
N = length(sorted_ws);
samples_above = (N-1):-1:0;
num_years = N/8760;
plot(samples_above/num_years,sorted_ws)


% Exercise 2

%% (a) import turbine power production data
parse_sam = @(x) double(regexp(x,'\|','split'));

ge.speeds = parse_sam("1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|23|24|25|26|27|28|29|30");
ge.power_kw = parse_sam("0|0|0|38|154|405|693|986|1323|1506|1578|1600|1600|1600|1600|1600|1600|1600|1600|1600|1600|1600|1600|1600|1600|0|0|0|0|0");
ge.hub_height = 80;
wind_data.WS80M = wind_data.WS50M * (80/50)^alpha;
wind_data.WD80M = wind_data.WD50M; % assume direction is the same
dist_ge = fitdist(wind_data.WS80M,'Weibull')

figure
ws_grid_ge = linspace(0,max(wind_data.WS80M),1000);
yyaxis left
plot(ge.speeds,ge.power_kw)
hold on
yyaxis right
histogram(wind_data.WS80M,'Normalization','pdf')
plot(ws_grid_ge,dist_ge.pdf(ws_grid_ge),'r-')
title('GE')

siemens.speeds = parse_sam("1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|23|24|25|26|27|28|29|30");
siemens.power_kw = parse_sam("0|0|0|127|285|512|795|1140|1587|2120|2770|3287|3500|3578|3600|3600|3600|3600|3600|3600|3600|3600|3600|3600|3600|0|0|0|0|0");
siemens.hub_height = 90;
wind_data.WS90M = wind_data.WS50M * (90/50)^alpha;
wind_data.WD90M = wind_data.WD50M; % assume direction is the same
dist_siemens = fitdist(wind_data.WS90M,'Weibull')

figure
ws_grid_siemens = linspace(0,max(wind_data.WS90M),1000);
yyaxis left
plot(siemens.speeds,siemens.power_kw)
hold on
yyaxis right
histogram(wind_data.WS90M,'Normalization','pdf')
plot(ws_grid_siemens,dist_siemens.pdf(ws_grid_siemens),'r-')
title('Siemens')

%% (c) power-duration curve
powerfun_ge = griddedInterpolant(ge.speeds,ge.power_kw);
power_ge = powerfun_ge(ws_grid_ge);
hours_exceeded_ge = 8760*(1-dist_ge.cdf(ws_grid_ge));

powerfun_siemens = griddedInterpolant(siemens.speeds,siemens.power_kw);
power_siemens = powerfun_siemens(ws_grid_siemens);
hours_exceeded_siemens = 8760*(1-dist_ge.cdf(ws_grid_siemens));

plot(hours_exceeded_ge,power_ge,'b-','linewidth',2,'DisplayName','GE (Weibull)')
hold on
plot(hours_exceeded_siemens,power_siemens,'r-','linewidth',2,'DisplayName','Siemens (Weibull)')
legend()

%% (d) AEP
aep_ge = compute_aep(ge,dist_ge)/ 1e6 % GWh
aep_siemens = compute_aep(siemens,dist_siemens)/ 1e6 % GWh

%%  Export to SAM
sam_export(wind_data,[datetime(2023,1,1,0,0,0),...
    datetime(2023,12,31,23,59,59)],-26.73,151.47,10);
clearvars, close all
file = "POWER_Point_Hourly_20150101_20241231_026d73S_151d47E_UTC.csv";

wind_data = import_nasa_power_wind(file);

pax = wind_rose(wind_data.WD50M);

% project up to 100m
a = 0.14;
wind_data.WS100M = wind_data.WS50M*(100/50)^a; % typical power law
ave_ws_100 = mean(wind_data.WS100M)
std_ws_100 = std(wind_data.WS100M)
wind_power_density = 0.5*1.2*ave_ws_100^3

dist_100 = fitdist(wind_data.WS100M,'Weibull');
ws_grid_100 = linspace(0,max(wind_data.WS100M),100);

figure
histogram(wind_data.WS100M,'Normalization','pdf','DisplayName','Data')
hold on
plot(ws_grid_100,dist_100.pdf(ws_grid_100),'DisplayName','Fit')
xlabel('Wind speed (m/s)')
ylabel('density')
legend()

% velocity-duration curve
sorted_ws = sort(wind_data.WS100M);
hours_exceeding = (length(sorted_ws)-1):-1:0;

figure
plot(hours_exceeding,sorted_ws)
xlabel('Hours')
ylabel('Wind speed (m/s)')



%% Exercise 2
parse_sam = @(s) double(regexp(s,'\|','split'));
ge.speeds = parse_sam("1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|23|24|25|26|27|28|29|30");
ge.power_kw = parse_sam("0|0|0|38|154|405|693|986|1323|1506|1578|1600|1600|1600|1600|1600|1600|1600|1600|1600|1600|1600|1600|1600|1600|0|0|0|0|0");
ge.hub_height = 80;
siemens.speeds = parse_sam("1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|23|24|25|26|27|28|29|30");
siemens.power_kw = parse_sam("0|0|0|127|285|512|795|1140|1587|2120|2770|3287|3500|3578|3600|3600|3600|3600|3600|3600|3600|3600|3600|3600|3600|0|0|0|0|0");
siemens.hub_height = 90;

wind_data.WS80M = wind_data.WS50M * (80/50)^a;
wind_data.WS90M = wind_data.WS50M * (90/50)^a;

dist_80 = fitdist(wind_data.WS80M,'Weibull');
ws_grid_80 = linspace(0,max(wind_data.WS80M),1000);
dist_90 = fitdist(wind_data.WS90M,'Weibull');
ws_grid_90 = linspace(0,max(wind_data.WS90M),1000);

figure
yyaxis left
plot(ge.speeds,ge.power_kw,'DisplayName','Power Curve')
ylim([0,4000])
ylabel('Power (kW)')
yyaxis right
hold on
histogram(wind_data.WS80M,'Normalization','pdf')
plot(ws_grid_80,dist_80.pdf(ws_grid_80))
xlabel('Wind Speed (w/s)')
hold off

figure
yyaxis left
plot(siemens.speeds,siemens.power_kw,'DisplayName','Power Curve')
ylim([0,4000])
ylabel('Power (kW)')
yyaxis right
hold on
histogram(wind_data.WS90M,'Normalization','pdf')
plot(ws_grid_90,dist_90.pdf(ws_grid_90))
xlabel('Wind Speed (w/s)')
hold off

% figure
% plot(ge.speeds,ge.power_kw,'DisplayName','GE')
% hold on
% plot(siemens.speeds,siemens.power_kw,'DisplayName','Siemens');
% legend()
% xlabel('wind speed (m/s)')
% ylabel('power (kW)')

% Power duration curve
powerfun_ge = griddedInterpolant(ge.speeds,ge.power_kw);
powerfun_siemens = griddedInterpolant(siemens.speeds,siemens.power_kw);

power_ge = powerfun_ge(ws_grid_80);
hours_exceeding_ge = 8760*(1-dist_80.cdf(ws_grid_80));

power_siemens = powerfun_siemens(ws_grid_90);
hours_exceeding_siemens = 8760*(1-dist_90.cdf(ws_grid_90));

figure
plot(hours_exceeding_ge,power_ge)
hold on
plot(hours_exceeding_siemens,power_siemens)
hold off
xlabel('Hours')
ylabel('Power (kW)')
legend("GE",'Siemens')
title('Power-duration curve')
hold on

% AEP
lb = 0;
ub = max(ge.speeds);
integrand =@(x) dist_80.pdf(x).*powerfun_ge(x);
aep_ge = 8760*integral(integrand,lb,ub)/1e6 % GWh

ub = max(siemens.speeds);
integrand =@(x) dist_90.pdf(x).*powerfun_siemens(x);
aep_siemens = 8760*integral(integrand,lb,ub)/1e6 % GWh

% 2023 data
wind_data_2023 = wind_data(year(wind_data.timestamp)==2023,:);

dist_2023_80 = fitdist(wind_data_2023.WS80M,'Weibull');
dist_2023_90 = fitdist(wind_data_2023.WS90M,'Weibull');

power_ge_2023 = powerfun_ge(ws_grid_80);
hours_exceeding_ge_2023 = 8760*(1-dist_2023_80.cdf(ws_grid_80));

power_siemens_2023 = powerfun_siemens(ws_grid_90);
hours_exceeding_siemens_2023 = 8760*(1-dist_2023_90.cdf(ws_grid_90));

figure
plot(hours_exceeding_ge_2023,power_ge_2023)
hold on
plot(hours_exceeding_siemens_2023,power_siemens_2023)
hold off
xlabel('Hours')
ylabel('Power (kW)')
legend("GE",'Siemens')
title('Power-duration curve (2023)')
hold on




%% Export to SAM
wind_data.WD80M = wind_data.WD50M; 
wind_data.WD90M = wind_data.WD50M; 
sam_export(wind_data,...
    [datetime(2023,1,1,0,0,0),datetime(2023,12,31,23,59,59)],...
    -27.73,151.47,10);


























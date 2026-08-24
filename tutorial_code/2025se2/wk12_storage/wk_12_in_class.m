clearvars, close all

%% Load in data and compute surplus
load = readtable('load.csv');
gen = readtable('generation.csv');
surplus_gen = gen.SystemPowerGenerated__kW_/1000 - ...
    load.Load_MW_; % MW

%% Plot to have a look at the data
figure
tiledlayout(3,1)
nexttile
plot(load.HourOfYear,load.Load_MW_,'DisplayName','Load');
hold on
plot(gen.HourOfYear,gen.SystemPowerGenerated__kW_/1000,'DisplayName','Generation')
% legend()
nexttile
plot(load.HourOfYear,surplus_gen,'DisplayName','Surplus');

sum(surplus_gen(surplus_gen>0))

%% Simulate the battery
eta_charge = 1.0;
eta_discharge = 1.0;
battery_capacity = 25; % MWh
start_soc = 1.0;

N = length(surplus_gen);
stored = zeros(N+1,1);
unsatisfied_demand = zeros(N,1);
curtailed = zeros(N,1);
stored(1) = battery_capacity * start_soc;

for ii=1:N
    % battery change in energy
    if surplus_gen(ii) > 0 % charge
        d = eta_charge * surplus_gen(ii);
    else % discharge
        d = 1/eta_discharge * surplus_gen(ii);
    end
    
    stored(ii+1) = stored(ii) + d;
    % battery limits
    if stored(ii+1)>=battery_capacity
        curtailed(ii) = stored(ii+1) - battery_capacity;
        stored(ii+1) = battery_capacity;
    elseif stored(ii+1)<=0
        unsatisfied_demand(ii) = -stored(ii+1);
        stored(ii+1) = 0;
    end
end

nexttile
plot(1:8761,stored,'DisplayName','Stored Energy')
ylabel('Energy (MWh)')
xlabel('Hour of the year')

fprintf('======= %.0f MWh Capacity ===\n',battery_capacity)
fprintf('Curtailment: %.1f\n',sum(curtailed))
fprintf('Unmet demand: %.1f\n',sum(unsatisfied_demand))





















%% plot and outputs 
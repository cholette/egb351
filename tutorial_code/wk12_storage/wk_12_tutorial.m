clearvars, close all

%% Load in data and compute residual load
load = readtable('load.csv');
gen = readtable('generation.csv');
surplus_gen = load;
surplus_gen.Load_MW_ = gen.SystemPowerGenerated__kW_/1000 - load.Load_MW_ ;

figure
tiledlayout(2,1)
nexttile
plot(load.HourOfYear,load.Load_MW_,'DisplayName','Load')
hold on
plot(gen.HourOfYear,gen.SystemPowerGenerated__kW_/1000,'DisplayName','Generation')
plot(surplus_gen.HourOfYear,surplus_gen.Load_MW_,'DisplayName','Surplus Generation')
ylabel('Power (MW)')
legend()

%% Simulate battery
eta_charge = 0.95;
eta_discharge = 0.9;
battery_capacity = 1000; % MWh
start_soc = 0;
surplus = surplus_gen.Load_MW_;

N_times = length(surplus);
stored = zeros(N_times+1,1);
unsatisfied_demand = zeros(N_times,1);
curtailed = zeros(N_times,1);
stored(1,:) = start_soc*battery_capacity; 
for ii = 1:N_times
    if surplus(ii)>0 % charge if possible
        d = eta_charge*surplus(ii);
    else % discharge if possible
        d = surplus(ii)/eta_discharge;
    end

    stored(ii+1) = stored(ii) + d;
    if stored(ii+1) >= battery_capacity
        curtailed(ii) = stored(ii+1) - battery_capacity; % curtail back to battery capacity 
        stored(ii+1) = battery_capacity;
    elseif stored(ii+1,:) <=0
        unsatisfied_demand(ii) = -stored(ii+1);
        stored(ii+1) = 0;
    end
end
max_discharge = max(-surplus);
max_charge = max(surplus);

nexttile
plot(1:8761,stored,'DisplayName','Stored Energy')
ylabel("Energy (MWh)")
xlabel('Hour of the year')

fprintf('Curtailment: %.1f MWh\n',sum(curtailed))
fprintf('Unsatisfied demand: %.1f MWh\n',sum(unsatisfied_demand))

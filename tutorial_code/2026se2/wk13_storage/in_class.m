clearvars

%% Load in data and compute residual load
load = readtable('load.csv');
gen = readtable('generation.csv');
surplus_gen = load;
surplus_gen.Load_MW_ = gen.SystemPowerGenerated__kW_/1000 - load.Load_MW_;

figure
tiledlayout(2,1)
nexttile
plot(load.HourOfYear,load.Load_MW_,'DisplayName','Load')
hold on
plot(gen.HourOfYear,gen.SystemPowerGenerated__kW_/1000.0,'DisplayName','Generation')
nexttile
plot(surplus_gen.HourOfYear,surplus_gen.Load_MW_)
legend()

%% Simulate battery
eta_c = 0.95;
eta_d = 0.9;
batt_cap = 1000; % MWh
start_soc = 0.5; % fraction of capacity
surplus = surplus_gen.Load_MW_;


N = length(surplus);
stored = zeros(N+1,1);
unsatisfied_demand = zeros(N,1);
curtailed = zeros(N,1);
stored(1) = start_soc*batt_cap;

for ii = 1:N
    if surplus(ii)>0
        d = eta_c * surplus(ii);
    else
        d = surplus(ii)/eta_d ;
    end

    stored(ii+1) = stored(ii) + d;

    % check limits
    if stored(ii+1) > batt_cap
        curtailed(ii) = stored(ii+1) - batt_cap;
        stored(ii+1) = batt_cap;
    elseif stored(ii+1) < 0
        unsatisfied_demand(ii) = -stored(ii+1);
        stored(ii+1) = 0;
    end
end
max_discharge = max(-surplus);
max_charge = max(surplus);

% plotting
hold on 
plot(1:8760,unsatisfied_demand,'DisplayName','Unsat. demand')
legend


figure
plot(1:8761,stored,'DisplayName','Stored Energy')
ylabel('Energy (MWh)')
xlabel('Hour of the year')

fprintf('Curtailment: %.1f MWh\n',sum(curtailed))
fprintf('Demand shortfall: %.1f MWh\n',sum(unsatisfied_demand))



















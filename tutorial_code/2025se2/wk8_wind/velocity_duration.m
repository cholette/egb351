function [fig,ax] = velocity_duration(ws)
% Hourly data with no leap years assumed

num_years = length(ws)/(24*365); 
speeds = linspace(min(ws),max(ws),1000);
counts = zeros(size(speeds));
for ii = 1:length(speeds)
    counts(ii) = sum(ws>=speeds(ii));
end

fig = figure();
plot(counts/num_years,speeds,'LineWidth',2)
ax = gca;
xlabel('Number of hours per year')
ylabel('Wind speed')
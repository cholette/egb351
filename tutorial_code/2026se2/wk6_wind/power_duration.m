function [powers,hours] = power_duration(turb,dist)
    
% isolate monotonic/invertable part
P_rated = max(turb.power_kw);
idx = find( (turb.power_kw<P_rated) & (turb.power_kw>0) );
g_inv = griddedInterpolant(turb.power_kw(idx),turb.speeds(idx));

powers = linspace(0,P_rated,1000);
hours = 8760*(1-dist.cdf(g_inv(powers)));
hours(end) = 0; % cannot produce more than rated power



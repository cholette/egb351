clearvars
close all
sigma = 5.67e-8; % [W/m2/K]
T_p = 120 + 273.15; % surface temperature [K]
T_sky = 263.15; % [K]
T_inf = 30+273.15; % [K]
G_s = 750; % W/m2
epsilon_p = 0.1;
epsilon_g =0.9;
alpha_s = 0.95;
k_air = 0; % air conductivity [W/m/K]
gap = 10e-3; % air gap [m]

fun = @(T_g) epsilon_g*sigma*(T_sky.^4 - T_g.^4 ) - ...
    0.22*sign(T_g-T_inf).*abs(T_g-T_inf).^(4/3) + ...
    0.0989*sigma*(T_p.^4-T_g.^4);

T_g = fzero(fun,T_sky);
t = linspace(T_sky,T_p,1000);
plot(t,fun(t))
hold on
plot(T_g,fun(T_g),'r*')
xlabel('T_g (K)')
ylabel('Residual for energy balance')
fprintf("T_g = %.2f K\n",T_g)

qwater = alpha_s*G_s - 0.0989*sigma*(T_p.^4-T_g.^4);

fprintf("q_water = %.2f W/m2 \n",qwater)



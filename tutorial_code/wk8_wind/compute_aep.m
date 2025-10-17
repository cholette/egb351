function aep = compute_aep(turb,dist,type,N,num_annual_timesteps)

    arguments
        turb struct;
        dist prob.WeibullDistribution;
        type string = "analytical";
        N double = 100000;
        num_annual_timesteps double = 8760; %
    end
    
    power = griddedInterpolant(turb.speeds,turb.power_kw);
    if lower(type)=='analytical'
        lb = min(turb.speeds);
        ub = max(turb.speeds);
        fun = @(x) dist.pdf(x).*power(x);
        aep = num_annual_timesteps*integral(fun,lb,ub);   
    elseif lower(type) == 'monte carlo' || lower(type)=='montecarlo'
        samples = power(dist.random(N,1));
        s = std(samples)/sqrt(N);
        m = mean(samples);
        aep = num_annual_timesteps*[m-1.96*s,m+1.96*s];
    end
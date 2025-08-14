function [GPOA,varargout] = plane_of_array_irradiance(beta,gamma,rho,...
    GHI,DNI,solar_zenith,solar_azimuth)
    % beta, gamma are panel tilt, aziumth (south convention), respectively
    % rho is the ground albedo
    % GHI, DNI, solar_zenith, solar_azimuth obvious

    DHI = GHI - DNI.*cosd(solar_zenith);

    % this can be done without another function in class
    AOI = angle_of_incidence([beta,gamma],[solar_zenith,solar_azimuth]);  
    
    direct = DNI.*cosd(AOI);
    direct(AOI>=90) = 0.0; % behind panel
    GPOA = direct + ...
        0.5*DHI*(1+cosd(beta)) + rho * GHI * 0.5*(1-cosd(beta));

    if nargout>1
        varargout{1} = AOI;
    end
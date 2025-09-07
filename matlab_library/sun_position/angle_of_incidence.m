function aoi = angle_of_incidence(panel_angles,solar_angles)
% panel angles are [tilt,azimuth] in degrees
% solar angles are [zenith,azimuth] in degrees

if size(panel_angles,1)>1
    if size(panel_angles,1)~=size(solar_angles,1)
        error('Angles of panel need to either be constant or one angle'+...
            'needs to be supplied for each solar angle.')
    end
    b = panel_angles(:,1);
    g = panel_angles(:,2);
else
    b = panel_angles(1);
    g = panel_angles(2);
end

z = solar_angles(:,1);
az = solar_angles(:,2);

aoi = acosd( cosd(z).*cosd(b) + sind(z).*sind(b).*cosd(az-g) );
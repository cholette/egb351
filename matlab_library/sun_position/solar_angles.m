function [z,az] = solar_angles(doy,hour_angle,lat)

delta = 23.45.*sind(360/365*(doy-81)); % degrees
z = acosd( sind(delta).*sind(lat) + ...
    cosd(delta).*cosd(lat).*cosd(hour_angle)  );

num = cosd(z).*sind(lat) - sind(delta);
den = sind(z).*cosd(lat);
az = sign(hour_angle) .* abs(acosd(num./den));

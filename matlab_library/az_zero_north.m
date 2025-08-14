function azNorth = az_zero_north(azSouth)

azNorth = nan(size(azSouth));
azNorth(azSouth<0) = -180 - azSouth(azSouth<0);  % east of north
azNorth(azSouth>=0) = 180 - azSouth(azSouth>=0);  % west of north

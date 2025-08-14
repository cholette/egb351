function E = equation_of_time(doy)

B = 360/364 *(doy-81); % in degrees
E = 9.87*sind(2*B) - 7.53*cosd(B) - 1.5*sind(B);

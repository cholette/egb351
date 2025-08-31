function row = wind_power_curve_to_sam(wind_speeds,power_kw,rated_power,...
    rotor_diameter,name)
%WIND_POWER_CURVE_TO_SAM Format wind turbine power curve for SAM library.
%
%   ROW = WIND_POWER_CURVE_TO_SAM(WIND_SPEEDS, POWER_KW, RATED_POWER,
%   ROTOR_DIAMETER) creates a CSV-formatted string representing a wind
%   turbine definition for the System Advisor Model (SAM) Wind Turbines
%   library file.
%
%   The output ROW can be appended to the file:
%       <SAM INSTALLATION DIRECTORY>\libraries\Wind Turbines.csv
%
%   Inputs:
%       WIND_SPEEDS   - Vector of wind speeds (m/s) at which the power
%                       curve is defined.
%
%       POWER_KW      - Vector of turbine power output (kW) corresponding
%                       to WIND_SPEEDS. Must be the same length.
%
%       RATED_POWER   - Rated power of the turbine in MW.
%
%       ROTOR_DIAMETER- Rotor diameter in meters.
%
%   Optional Input:
%       NAME          - Turbine name string for identification in SAM.
%                       Default: 'selected_turbine'.
%
%   Output:
%       ROW           - A single CSV-formatted string containing the
%                       turbine definition in SAM’s expected format:
%
%            Name, RatedPower(W), RotorDiameter(m), Manufacturer,
%            WindSpeeds(m/s, pipe-delimited), Power(W, pipe-delimited)
%
%   Notes:
%       • Rated power is converted from MW to kW internally.
%       • Manufacturer is set to "unknown" by default.
%       • The pipe symbol "|" is used to delimit wind speeds and powers
%         as required by SAM.
%       • To use the result, append ROW to the file
%         <SAM INSTALLATION DIRECTORY>\libraries\Wind Turbines.csv
%         and restart SAM.
%
%   Example:
%       ws = 0:2:25;                        % wind speeds (m/s)
%       pcurve = [0 0 50 200 500 900 1500 ... % power curve (kW)
%                 2000 2500 3000 3300 3400 3450 3450];
%       row = wind_power_curve_to_sam(ws, pcurve, 3.4, 90, ...
%                                     "MyTurbine");
%
%   See also: JOIN, FPRINTF

arguments
    wind_speeds
    power_kw
    rated_power
    rotor_diameter
    name = 'selected_turbine';
    
end
speeds = join(string(wind_speeds),"|");
powers = join(string(power_kw),"|");
row = join([name,rated_power*1000,rotor_diameter,"unknown",speeds,powers],',');

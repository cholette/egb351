function wind_and_solar_data = import_nasa_power(wind_and_solar_file)
%IMPORT_NASA_POWER Import NASA POWER single point data into MATLAB
%
%   wind_and_solar_data = IMPORT_NASA_POWER(wind_and_solar_file)
%   reads a NASA POWER | Data Access Viewer (DAV) wind and solar
%   data from a Single Point download file and converts it into a MATLAB
%   table with timestamps and corrected variables.
%
%   Input:
%       wind_and_solar_file - String or character vector specifying the
%                             path to the NASA POWER DAV single point CSV
%                             file.
%
%   Output:
%       wind_and_solar_data - MATLAB table containing the imported
%                             variables. The table includes a new
%                             'timestamp' column constructed from
%                             YEAR, MO, DY, and HR fields. Negative values
%                             of ALLSKY_SFC_SW_DNI are set to NaN, and
%                             surface pressure PS is converted from kPa
%                             to Pa.
%
%   Notes:
%       - The function automatically detects the end of the file header
%         (line containing '-END HEADER-') to configure import options.
%       - Time columns are combined into a MATLAB datetime variable.
%       - Unit conversion: PS [kPa] → PS [Pa].
%
%   Example:
%       data = import_nasa_power('POWER_singlepoint.csv');
%       plot(data.timestamp, data.ALLSKY_SFC_SW_DNI);

opts_solar = detectImportOptions(wind_and_solar_file);

contents = readlines(wind_and_solar_file);
end_header = find(contains(contents,"-END HEADER-"),1);

opts_solar.DataLines = [end_header+2 Inf];
opts_solar.VariableNamesLine = end_header+1;
wind_and_solar_data = readtable(wind_and_solar_file,opts_solar);

% add a datetime column for convenience
wind_and_solar_data.timestamp = datetime( wind_and_solar_data.YEAR, ...
                                wind_and_solar_data.MO, ...
                                wind_and_solar_data.DY, ...
                                wind_and_solar_data.HR, ...
                                zeros(size(wind_and_solar_data,1),1), ...
                                zeros(size(wind_and_solar_data,1),1));

nan_idx = find(wind_and_solar_data.ALLSKY_SFC_SW_DNI<0);
wind_and_solar_data.ALLSKY_SFC_SW_DNI(nan_idx) = nan;
wind_and_solar_data.PS = 1000*wind_and_solar_data.PS; % kPa -> Pa
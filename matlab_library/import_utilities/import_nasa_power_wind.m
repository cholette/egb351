function wind_data = import_nasa_power_wind(csv_file)
%import_nasa_power_wind Import NASA POWER single point wind data into MATLAB
%
%   wind_data = IMPORT_NASA_POWER_WIND(nasa_power_dav_csv)
%   reads a NASA POWER | Data Access Viewer (DAV) wind
%   data from a Single Point download file and converts it into a MATLAB
%   table with timestamps and corrected variables.
%
%   Input:
%       nasa_power_dav_csv  - String or character vector specifying the
%                             path to the NASA POWER | DAV single point CSV
%                             file.
%
%   Output:
%       wind_data           - MATLAB table containing the imported
%                             variables. The table includes a new
%                             'timestamp' column constructed from
%                             YEAR, MO, DY, and HR fields. Surface pressure 
%                             PS is converted from kPa to Pa.
%                            
%
%   Example:
%       data = import_nasa_power('POWER_singlepoint.csv');
%       plot(data.timestamp, data.WS50M);

opts_solar = detectImportOptions(csv_file);

contents = readlines(csv_file);
end_header = find(contains(contents,"-END HEADER-"),1);

opts_solar.DataLines = [end_header+2 Inf];
opts_solar.VariableNamesLine = end_header+1;
wind_data = readtable(csv_file,opts_solar);

% check to ensure that all columns are available
required = {'YEAR','MO','DY','HR','T2M','PS','WS10M','WD10M','WS50M','WD50M'};
req_names = {'Year','Month','Day','Hour','Temperature at 2M','Surface Pressure',...
    'Wind speed at 10M','Wind direction at 10M',...
    'Wind speed at 50M','Wind direction at 50M'};
vars = wind_data.Properties.VariableNames;
for ii = 1:length(required)
    if ~any(contains(vars,required{ii}))
        error("Missing variable "+required{ii} +...
            ". Go back to NASA POWER | DAV and select " + req_names{ii})
    end
end

% add a datetime column for convenience
wind_data.timestamp = datetime( wind_data.YEAR, ...
                                wind_data.MO, ...
                                wind_data.DY, ...
                                wind_data.HR, ...
                                zeros(size(wind_data,1),1), ...
                                zeros(size(wind_data,1),1));

wind_data.PS = 1000*wind_data.PS; % kPa -> Pa
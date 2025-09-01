function [load_data,wind_and_solar_data,costs] = get_project_data( ...
    data_dir,student_number,base_load_file,wind_and_solar_file)
% GET_MY_DATA Load demand, resource, and cost data EGB351 project.
%
%   [LOAD_DATA, WIND_AND_SOLAR_DATA, COSTS] = GET_MY_DATA(DATA_DIR,
%   STUDENT_NUMBER, LOAD_FILE, WIND_AND_SOLAR_FILE) reads and processes
%   demand (load) data, NASA POWER|DAV Single POint wind and solar resource 
%   data,(see:https://power.larc.nasa.gov/data-access-viewer/) and
%   provides default cost assumptions for wind, solar, and battery
%   technologies.
%
%   Inputs:
%       DATA_DIR           - String. Directory path containing the input
%                            data files. Must include trailing file
%                            separator "/".
%
%       STUDENT_NUMBER     - Your student number.
%
%       BASE_LOAD_FILE     - String. File name of the Excel/CSV file
%                            containing base load data.
%
%       WIND_AND_SOLAR_FILE- String. File name of the NASA POWER | DAV
%                            CSV file containing solar and wind
%                            resource data.
%
%   Outputs:
%       LOAD_DATA          - Table with two columns:
%                               timestamp : datetime vector
%                               load_MW   : demand in MW
%
%       WIND_AND_SOLAR_DATA- Table containing NASA POWER solar and wind
%                            variables with an added timestamp column.
%                            • Negative DNI values are replaced with NaN.
%                            • Surface pressure (PS) is converted from kPa
%                              to Pa.
%
%       COSTS              - Struct with default technology costs:
%                               costs.wind
%                               costs.solar
%                               costs.battery
%                               costs.fcr (fixed charge rate, 1/yr)
%
%   Notes:
%       • Load data are scaled by a factor depending on STUDENT_NUMBER.
%         This ensures unique variations for different users while keeping
%         overall magnitudes realistic.
%       • WIND_AND_SOLAR_FILE is assumed to follow NASA DAV formatting,
%         with -END HEADER- denoting the beginning of the data.
%       • All cost data are illustrative and would need updating for real 
%         projects.
%
%   Example:
%       % Load my datasets with deterministic scaling
%       [loadTbl, resourceTbl, costs] = get_my_data("data/", 0, ...
%                       "queensland_load.xlsx", "nasa_power_data.csv");
%
%       % Plot load profile
%       plot(loadTbl.timestamp, loadTbl.load_MW);
%       ylabel('Load (MW)'); xlabel('Time');

arguments
    data_dir(1,1) string;
    student_number(1,1) int64;
    base_load_file(1,1) string;
    wind_and_solar_file string;
end

[load_data,wind_and_solar_data,costs] = get_my_data( ...
    data_dir,student_number,base_load_file,wind_and_solar_file);
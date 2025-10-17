function new_data = localize_time(nasa_data,new_offset,original_offset)
%LOCALIZE_TIME Adjusts timestamps in a dataset to a specified timezone offset.
%
%   new_data = LOCALIZE_TIME(nasa_data, new_offset, original_offset)
%
%   This function converts the timestamps in the input dataset `nasa_data`
%   from the original timezone offset (default is UTC) to a new timezone
%   offset specified by `new_offset`. It also extracts and adds separate
%   columns for year, month, day, and hour based on the localized time.
%
%   Inputs:
%       nasa_data       - MATLAB table containing NASA POWER data with
%                         timestamps and meteorological variables. This table
%                         should be created originally from get_project_data.
%
%       new_offset      - The desired timezone offset (e.g., -5 for EST, +9 for JST).
%       original_offset - (Optional) The original timezone offset of the timestamps.
%                         Default is 0 (UTC).
%
%   Output:
%       new_data        - A updated table with timestamps adjusted to the new timezone.
%
%   Example:
%       localized_data = localize_time(data, -5); % Convert from UTC to EST
%
%   Notes:
%       - Timezone offsets are interpreted as hours relative to UTC.
%       - The timestamp column is temporarily assigned a timezone for conversion,
%         and then cleared to avoid timezone conflicts in further processing.


    arguments
        nasa_data
        new_offset
        original_offset = 0 % UTC
    end

    % localize timestamps to new timezone
    new_data = nasa_data;
    original_offset = offset_str(original_offset);
    new_data.timestamp = datetime(new_data.timestamp,TimeZone=original_offset);
    new_data.timestamp.TimeZone = offset_str(new_offset);

    % put into year,month,day columns
    new_data.YEAR = year(new_data.timestamp);
    new_data.MO = month(new_data.timestamp);
    new_data.DY = day(new_data.timestamp);
    new_data.HR = hour(new_data.timestamp);
    new_data.timestamp.TimeZone = '';
    
    function os = offset_str(v)
       if v>=0
           os = "+"+num2str(v);
       else
           os = "-"+num2str(v);
       end
    end
end

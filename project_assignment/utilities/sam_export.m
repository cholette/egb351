function sam_export(dat,date_range,lat,lon,tz,name,type,elevation)
%SAM_EXPORT Export NASA POWER | DAV data into SAM-compatible weather files.
%
%   SAM_EXPORT(DAT, DATE_RANGE, LAT, LON, TZ) exports NASA POWER | DAV data
%   (https://power.larc.nasa.gov/data-access-viewer/) contained in the 
%   table DAT into a file formatted for the System Advisor Model (SAM).
%
%   Inputs:
%       DAT        - MATLAB table containing NASA POWER data with
%                    timestamps and meteorological variables. This table
%                    should be created originally from get_project_data.
%                    
%
%       DATE_RANGE - Two-element vector [START, END] specifying the 
%                    datetime range of interest. Only rows within this 
%                    range are exported.
%
%       LAT        - Site latitude in decimal degrees.
%
%       LON        - Site longitude in decimal degrees.
%
%       TZ         - Time zone offset from UTC (hours). Both site and 
%                    data timezones are set to this value.
%
%   Optional Inputs (Name-Order Arguments):
%       NAME       - Output file name (default: "sam_export.csv").
%
%       TYPE       - Data type, either:
%                        "wind"  (default) - Export wind resource data
%                        "solar"           - Export solar resource data
%
%       ELEVATION  - Site elevation above sea level in meters.
%                    (Used only for solar exports; default: 0)
%
%   Notes:
%       • For TYPE="wind":
%           - Wind speed/direction variable names are reformatted with
%             heights in meters.
%           - Air pressure is assumed at 10m (warning displayed).
%           - Temperature at 2m is used as a proxy for 10m (warning displayed).
%
%       • For TYPE="solar":
%           - Extracts GHI, DNI, DHI, temperature, wind speed, and pressure.
%           - Pressure is converted from Pa to mbar.
%           - Adds SAM-required latitude, longitude, timezone, and elevation
%             header rows.
%
%   Example:
%       % Export NASA POWER table "dat" between Jan 1 and Dec 31, 2020
%       % for a site at 40°N, 105°W, UTC-7, as solar resource data:
%
%       sam_export(dat, [datetime(2020,1,1), datetime(2020,12,31)], ...
%                  40, -105, -7, "mysite_solar.csv", "solar", 1700);
%
%   The resulting CSV file can be imported directly into SAM.
%
%   See 
%
%   See also: WRITETABLE, WRITECELL

arguments
    dat
    date_range
    lat
    lon 
    tz
    name = "sam_export.csv";
    type = "wind";
    elevation = 0;
end

dat = dat(  dat.timestamp>=date_range(1) &...
            dat.timestamp<=date_range(2),:);
switch lower(type)
    case "wind"
        row1 = ['Longitude, ',num2str(lon),', Latitude, ', num2str(lat), ...
            ', Site Timezone, ', num2str(tz),', Data Timezone, ',num2str(tz)];
        
        varnames = dat.Properties.VariableNames;
        % vars = ["YEAR","MO","DY","HR","PS","T2M"];
        newNames = varnames;
        for ii = 1:length(varnames)
            vn = upper(varnames{ii});
            if contains(vn,"WS")
                e = strfind(vn,"M");
                hh = vn(3:e-1);
                newNames{ii} = "Wind speed at " + hh + "m";
            elseif contains(vn,"WD")
                e = strfind(vn,"M");
                hh = vn(3:e-1);
                newNames{ii} = "Wind direction at " + hh + "m";
            elseif contains(vn,"YEAR")
                newNames{ii} = "Year";
            elseif contains(vn,"MO")
                newNames{ii} = "Month";
            elseif contains(vn,"DY")
                newNames{ii} = "Day";
            elseif contains(vn,"HR")
                newNames{ii} = "Hour";
            elseif contains(vn,"PS")
                newNames{ii} = "Air pressure at 10m (Pa)";
            elseif contains(vn,"T2M")
                newNames{ii} = "Air Temperature at 10m (C)";
            end
        end
        
        newNames = cellfun(@char,newNames,'UniformOutput',false); 
        dat = renamevars(dat,varnames,newNames);
        
        disp("Warning: Arbitrarily setting pressure data to be at 10m.")
        disp("Warning: Temperature at 10m is not available in MERRA-2 data. Using temperature at 2m instead.")
        disp("Writing wind data to "+name)
        
        writetable(dat,name)
        T = readcell(name);
        prepend = cell(1,size(T,2)); prepend{1}=row1;
        writecell([prepend;T],name,"QuoteStrings","none")

    case "solar"
        row1 = 'Longitude, Latitude, Timezone, Elevation';
        row2 = join(string([lon,lat,tz,elevation]),',');
        
        % headers = [ 'Year','Month','Day','Hour','Minute','GHI','DNI',...
                    % 'DHI','Temperature'];
        varnames = dat.Properties.VariableNames;
        newNames = varnames;
        cols = 1:4;
        for ii = 5:length(varnames)
            vn = upper(varnames{ii});
            if contains(vn,'ALLSKY_SFC_SW_DWN')
                newNames{ii} = 'GHI';
                cols(end+1)=ii;
            elseif contains(vn,"ALLSKY_SFC_SW_DNI")
                newNames{ii} = 'DNI';
                cols(end+1)=ii;
            elseif contains(vn,"ALLSKY_SFC_SW_DIFF")
                newNames{ii} = 'DHI';
                cols(end+1)=ii;
            elseif contains(vn,'T2M')
                newNames{ii} = 'Temperature';
                cols(end+1)=ii;
            elseif contains(vn,"WS10M")
                newNames{ii} = 'Wind Speed';
                cols(end+1)=ii;
            elseif contains(vn,"PS")
                newNames{ii} = 'Pressure';
                dat.(vn) = dat.(vn)/100; % Pa -> mbar
                cols(end+1)=ii;
            end
        end
        dat = renamevars(dat,varnames,newNames);
        new = dat(:,cols);
        newNames = cellfun(@char,newNames,'UniformOutput',false); 
        

        disp("Writing solar data to "+name)
        
        writetable(new,name)
        T = readcell(name);
        prepend = cell(2,size(T,2)); 
        prepend{1,1} = row1;
        prepend{2,1} = row2;
        writecell([prepend;T],name,"QuoteStrings","none")
end


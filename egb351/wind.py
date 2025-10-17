import scipy.stats as stats
from scipy.optimize import fsolve,minimize_scalar
from scipy.optimize import minimize,LinearConstraint
from scipy.integrate import quad
import numpy as np
from numpy import log
import matplotlib.pyplot as plt
from windrose import WindroseAxes
import pandas as pd
from turbine_models.parser import Turbines
import csv

# References
# [1] Wind Energy Explained, 1st ed. John Wiley & Sons, Ltd, 2009. doi: 10.1002/9781119994367.

class merra_wind_speed_data:
    def __init__(self):
        self.data_source_url = "https://power.larc.nasa.gov/data-access-viewer/"
        self.latitude = None
        self.longitude = None
        self.date_range = None
        self.average_elevation = None
        self.data = None

    def import_data(self,merra_single_point_csv):
        file = merra_single_point_csv
        
        # read header
        with open(file,'r') as f:
            reader = csv.reader(f)
            header_count = 0
            for row in reader:
                header_count += 1
                if row and row[0] == '-END HEADER-':
                    break
                
                if "location" in row[0].lower():
                    r_split = row[0].split()
                    islat = [('latitude' in r.lower()) for r in r_split]
                    if any(islat):
                        self.latitude = float(r_split[islat.index(True)+1])
                    
                    islon = [('longitude' in r.lower()) for r in r_split]
                    if any(islon):
                        self.longitude = float(r_split[islon.index(True)+1])
                
                if "elevation" in row[0].lower():
                    self.average_elevation = row[0].split('=')[-1]

        # read the rest
        df = pd.read_csv(file,skiprows=header_count)
        df['Timestamp'] = pd.to_datetime(dict(
                                year=df['YEAR'],
                                month=df['MO'],
                                day=df['DY'],
                                hour=df['HR'],
                                minute=df.get('minute', 0),  # use 0 if not present
                                second=df.get('second', 0),  # use 0 if not present
                                ))
        df = df.drop(['YEAR','MO','DY','HR'],axis=1)
        df = df[['Timestamp',*df.columns[:-1]]]
        self.data = df

    def add_speed_at_height(self,
                        heights:list[int],
                        z0=None,
                        alpha=None,
                        model="power_law",
                        plot=False):

        assert model.lower() in ['power_law','logarithmic'], 'Invalid model. Must be "power_law" or "logarithmic"'
        
        # ensure it is a list
        if isinstance(heights,int):
            heights = [heights]

        ws10,ws50= self.data['WS10M'], self.data['WS50M']
        wd10,wd50= self.data['WD10M'], self.data['WD50M']
        if model.lower() == 'logarithmic':
            if z0 is None:
                print('z0 not provided. Will try to fit it based on data available from 10m and 50m')
                p = self.fit_logarithmic()
            else:
                p = z0
                            
            for height in heights:
                new_col = f'WS{height}M'
                scale = log(height/p)/log(50/p)
                self.data[new_col] = scale * ws50
                self.data[f'WD{height}M'] = stats.circmean(np.c_[wd10,wd50],low=0,high=360,axis=1)

        else: # model.lower() == 'power_law':
            if alpha is None:
                print('alpha not provided. Will try to fit it based on data available from 10m and 50m')
                p = self.fit_power_law()
            else:
                p = alpha

            for height in heights:
                new_col = f'WS{height}M'                
                scale = (height/50.0)**p
                self.data[new_col] = scale * ws50
                self.data[f'WD{height}M'] = stats.circmean(np.c_[wd10,wd50],low=0,high=360,axis=1)

        if plot:
            fig,ax = plt.subplots()
            ax.plot(ws50,ws10,'o')
            xl = ax.get_xlim()
            x = np.linspace(xl[0],xl[1],1000)
            if model.lower() == 'logarithmic':
                ax.plot(x,log(10.0/p)/log(50.0/p) * x,linewidth=3,label="Logarithmic fit")
            else:
                ax.plot(x,(10.0/50.0)**p *x,linewidth=3,label='Power law fit')
            ax.set_xlabel('Wind speed at 50m [m/s]')
            ax.set_ylabel('Wind speed at 10m [m/s]')
            ax.legend()
    
    def fit_power_law(self):
        ws10,ws50= self.data['WS10M'], self.data['WS50M']
        rss = lambda log_alpha: np.sum((ws10 - power_law(ws50,50,10,np.exp(log_alpha)) )**2)
        alpha = np.exp(minimize_scalar(rss)['x'])
        print(f'alpha = {alpha}')
        return alpha
    
    def fit_logarithmic(self):
        ws10,ws50= self.data['WS10M'], self.data['WS50M']
        rss = lambda log_z0: np.sum((ws10 - logarithmic(ws50,50,10,np.exp(log_z0)) )**2)
        z0 = np.exp(minimize_scalar(rss)['x'])
        print(f'z0 = {z0:.3e}m')
        return z0
    
    def export_to_sam_csv(self,site_timezone,filename:str):
        # See SAM CSV Format for Wind for more details on format. 

        # header
        row1 = ['Site Timezone',site_timezone,
                'Data Timezone',site_timezone,
                'Latitude',self.latitude,
                'Longitude',self.longitude,
                'Elevation',self.average_elevation]
        row1 = [str(r) for r in row1]
        row1 = ','.join(row1)
        
        # rename wind speed and direction columns
        heights = [name[2:-1] for name in self.data.columns if 'WS' in name]
        df2 = self.data.rename(columns={f'WS{h}M':f'wind speed at {h}m (m/s)' for h in heights})
        df2.rename(columns={f'WD{h}M':f'wind direction at {h}m (degrees)' for h in heights},inplace=True)
        df2['PS']  = df2['PS'] * 1000 # convert from kPa to Pa

        print("Warning: Arbitrarily setting pressure data to be at 10m.")
        df2.rename(columns={'PS':'air pressure at 10m (Pa)'},inplace=True)

        print("Warning: Temperature at 10m is not available in MERRA-2 data. Using temperature at 10m instead.")
        df2.rename(columns={'T2M':'air temperature at 10m (C)'},inplace=True)
               
        # add back columns for year, month, day, hour, minute if not there
        def prepend(c,v): 
            if c not in df2.columns: 
                df2.insert(0,c,v)
        prepend('Minute',df2.Timestamp.dt.minute)
        prepend('Hour',df2.Timestamp.dt.hour)
        prepend('Day',df2.Timestamp.dt.day)
        prepend('Month',df2.Timestamp.dt.month)
        prepend('Year',df2.Timestamp.dt.year)
        
        with open(filename,'w',newline='') as f:
            f.write(row1+'\n')
            df2.to_csv(f,header=True,index=False)        

def power_law(ws,height,new_height,alpha):
    return ws * (new_height/height)**alpha

def logarithmic(ws,height,new_height,z0):
    return ws * log(new_height/z0)/log(height/z0)
    
def speed_fit(data,plot=False,type='Weibull',hist_kwargs={},fit_kwargs={},handles=None):
    
    assert type.lower() is not None, 'Must provide a distribution type'
    assert type.lower() in ['weibull','rayleigh'], 'Invalid distribution type. Must be "weibull" or "rayleigh'
    if type.lower() == 'weibull':
        shape,loc,scale = stats.weibull_min.fit(data,floc=0) # don't fit location parameter
        dist = stats.weibull_min(shape,scale=scale)
    elif type.lower() == 'rayleigh':
        loc,scale = stats.rayleigh.fit(data,floc=0)
        dist = stats.rayleigh(scale=scale)
    else:
        raise ValueError('Invalid distribution type. Must be "weibull" or "rayleigh')
    
    if (not plot) and (handles is not None):
        print('Not adding to plot since plot = False')
        handles = None

    if plot:
        if handles is None:
            fig,ax = plt.subplots()
        else:
            fig,ax = handles
        ax.hist(data,label='Data',density=True,**hist_kwargs)
        XL = ax.get_xlim()
        x = np.linspace(XL[0],XL[1],1000)
        ax.plot(x,dist.pdf(x),label=f'{type} fit',**fit_kwargs)
        ax.set_xlabel('Wind speed')
        ax.set_ylabel('Density')
        ax.legend()
        return dist,fig,ax
    else:
        return dist

def change_height(dist,height_old,height_new,alpha,type='weibull'):
    # Modify scale parameter for different height. See notes on wind for details.
    if type.lower() == 'weibull':
        return stats.weibull_min(dist.args[0],scale=dist.kwds['scale']*(height_new/height_old)**alpha)
    elif type.lower() == 'rayleigh':
        return stats.rayleigh(scale=dist.kwds['scale']*(height_new/height_old)**alpha)
    else:
        raise ValueError('Invalid distribution type. Must be "weibull" or "rayleigh')

def wind_rose(direction,speed,num_bins=16,s_units=None):
    ax = WindroseAxes.from_ax()
    ax.bar(direction,speed,nsector=num_bins,normed=True)
    ax.set_yticklabels([f"{a:.2f}%" for a in ax.get_yticks()])
    ax.legend(title=f'Wind speed {s_units}',loc='best')
    return ax

def power_cdf(dist,power_curve,n_bins=10,power_units='W'):
    wmax = dist.ppf(0.999)
    w = np.linspace(0,wmax,n_bins+1)
    prob_w = dist.cdf(w)
    power_ = power_curve.get_power(w,power_units=power_units)

    power = np.unique(power_)
    prob = np.zeros_like(power)
    for ii,p in enumerate(power):
        prob[ii] = np.max(prob_w[power_== p])

    return np.r_[-1,power],np.r_[0,prob]
class speed_and_direction_dist:
    def __init__(self,type=None,shapes=[],scales=[],n_az_bins=[],probabilities=[]):

        assert type is not None, 'Must provide a distribution type'
        assert type.lower() in ['weibull','rayleigh'], 'Invalid distribution type. Must be "weibull" or "rayleigh'
        if type.lower() == 'weibull':
            self.dists = [stats.weibull_min(shape,scale=scale) for shape,scale in zip(shapes,scales)]
        elif type.lower() == 'rayleigh':
            assert shapes == [], 'Rayleigh distribution does not have shape parameter'
            self.dists = [stats.rayleigh(scale=scale) for scale in scales]
        else:
            raise ValueError('Invalid distribution type. Must be "weibull" or "rayleigh')
        dth = 360/(n_az_bins)
        azbin = np.arange(-dth/2,360-dth/2,dth)
        self.azimuth_bin_edges = azbin # azimuth is [0,360) with east-of-north positive
        self.azimuth_bin_centers = self.azimuth_bin_edges + dth/2
        azbin += dth
        self.probabilities = probabilities
        self.height = None
        self.type = type.lower()

    def set_height(self,height):
        self.height = height

    def change_height(self,new_height,alpha):
        assert self.height is not None, 'Must set height of the original data using set_height()'
        
        p = self.probabilities
        if self.type == 'weibull':
            for ii,d in enumerate(self.dists):
                self.dists[ii] = change_height(d,self.height,new_height,alpha,self.type)

    def rvs(self,N=1):
        rng = np.random.default_rng()
        direction_bin = rng.choice(len(self.probabilities),size=N,p=self.probabilities,replace=True)
        speeds = np.array([self.dists[ii].rvs() for ii in direction_bin])
        directions = self.azimuth_bin_centers[direction_bin]
        return directions,speeds
    
    def get_params(self):

        probs = self.probabilities
        scales,shapes = [],[]
        for d in self.dists:
            shapes.append(d.args[0])
            scales.append(d.kwds['scale'])

        return np.array(probs),np.array(scales),np.array(shapes)

    def pdf(self,x):
        return np.sum([self.probabilities[ii]*d.pdf(x) for ii,d in enumerate(self.dists)],axis=0)

    def cdf(self,x):
        return np.sum([self.probabilities[ii]*d.cdf(x) for ii,d in enumerate(self.dists)],axis=0)
    
    def mean_speed(self):
        return np.sum([self.probabilities[ii]*d.mean() for ii,d in enumerate(self.dists)])
    
    def var(self):
        return np.sum([self.probabilities[ii]*(d.var()+d.mean()**2) for ii,d in enumerate(self.dists)]) - self.mean()**2  # see https://en.wikipedia.org/wiki/Mixture_distribution#Moments 

    def ppf(self,p):
        return fsolve(lambda x: self.cdf(x)-p,self.mean())

    def fit(self,direction,speed,az_edges=None, plot=False,hist_kwargs={},fit_kwargs={},handles=None):
        # data is [direction,speed]
        assert (self.azimuth_bin_edges is not []) or (az_edges is not None), 'Must provide azimuth bin edges either in object or as a keyword argument to the fitting'

        if self.azimuth_bin_edges is []:
            self.azimuth_bin_edges = az_edges
            self.azimuth_bin_centers = az_edges[:-1] + np.diff(az_edges)/2

        # fits a distribution to the wind data for each direction
        dir,sp = direction,speed
        dir_bins = np.digitize(dir,self.azimuth_bin_edges,right=True)
        dir_bins[dir_bins == len(self.azimuth_bin_edges)] = 0 # wrap around to zero
        prob_hat = np.bincount(dir_bins,minlength=len(self.probabilities))
        self.probabilities = prob_hat/np.sum(prob_hat)

        for ii in np.unique(dir_bins):
            d = speed_fit(sp[dir_bins==ii],plot=False)
            self.dists.append(d)

    def wind_speed_curve(self,x,hours_per_year=8760,plot=True):
        R = 1-self.cdf(x)

        if plot:
            fig,ax = plt.subplots()
            ax.plot(R*hours_per_year,x)
            ax.set_xlabel('Hours per year')
            ax.set_ylabel('Wind speed')
        else:
            fig = None
            ax = None

        return R,fig,ax

    def mean_direction(self):
        # Circular mean of the wind direction
        m = np.angle( np.sum(np.exp(np.pi/180 * self.azimuth_bin_centers * 1j) * self.probabilities))
        if m > 0:  # the angle returns the principal value between -np.pi and np.pi, but I want 0 to 360
            return m * 180/np.pi
        else: 
            return 360 + m*180/np.pi
class turbine:
    def __init__(self):
        self.wind_speeds = []       # m/s
        self.power = []             # W
        self.Cp = []                # [-]
        self.Ct = []                # [-]
        self.thrust = []            # N
        self.cut_in_speed = None    # m/s
        self.cut_out_speed = None   # m/s
        self.rotor_diameter = None  # m
        self.hub_height = None      # m
        self.rated_power = None     # W

    def set_performance(self,wind_speeds,power=None,Cp=None,
                        Ct=None,thrust=None):
        self.wind_speeds = wind_speeds
        self.power = power
        self.Cp = Cp
        self.Ct = Ct
        self.thrust = thrust
    
    def import_nrel_power_curve(self,turbine_name):
        # From https://github.com/NREL/turbine-models

        turb = Turbines()
        specs = turb.specs(turbine_name)
        cols = [c.lower() for c in specs['power_curve'].columns]

        def find_col(substr):
            mask = [substr in c for c in cols]
            if True in mask:
                return cols[mask.index(True)]
            else:
                return None

        self.wind_speeds = specs['power_curve']['wind_speed_ms'].values
        self.power = specs['power_curve']['power_kw'].values*1000 # NREL data in kW
        self.Cp = specs['power_curve']['cp'].values
        self.cut_in_speed = specs['cut_in_wind_speed']
        self.cut_out_speed = specs['cut_out_wind_speed']
        self.rotor_diameter = specs['rotor_diameter']
        self.hub_height = specs['hub_height']
        self.rated_power = specs['rated_power']*1000 # NREL data in kW

        # import thrust if available
        thrust_column = find_col('thrust')
        ct_column = find_col('ct')
        if thrust_column is not None:
            self.thrust = specs['power_curve'][thrust_column].values
        else:
            print('No thrust data. Skipping import.')

        if ct_column is not None:
            self.Ct = specs['power_curve'][ct_column].values
        else:
            print('No thrust data. Skipping import.')
                                               
    def plot(self,ax=None,nonzero_only=False,power_units='MW',plt_kwargs={}):
        if ax is None:
            fig,ax = plt.subplots()
        
        if nonzero_only:
            x = self.wind_speeds
            fx = self.power
        else:
            ws = self.wind_speeds
            x = np.linspace(0,ws[-1]+0.1*(max(ws)-min(ws)),1000)
            fx = self.get_power(x)

        if power_units.lower() == 'w':
            ax.plot(x,fx,**plt_kwargs)
            ax.set_ylabel('Power [W]')
        elif power_units.lower() == 'kw':
            ax.plot(x,fx/1e3,**plt_kwargs)
            ax.set_ylabel('Power [kW]')
        elif power_units.lower() == 'mw':
            ax.plot(x,fx/1e6,**plt_kwargs)
            ax.set_ylabel('Power [MW]')
        else:
            raise ValueError('Invalid units. Must be "W", "kW", or "MW"')
        
        ax.set_xlabel('Wind speed [m/s]')
        return ax

    def get_power(self,wind_speed,power_units = 'W'):
        y = np.interp(wind_speed,self.wind_speeds,self.power)

        if power_units.lower() == 'w':
            scale = 1.0
        elif power_units.lower() == 'kw':
            scale = 1e-3
        elif power_units.lower() == 'mw':
            scale = 1e-6
        else:
            raise ValueError('Invalid units. Must be "W", "kW", or "MW"')

        if isinstance(wind_speed, np.ndarray):
            mask = (wind_speed <= self.cut_in_speed) | (wind_speed > self.cut_out_speed)
            y[mask]=0
        else:
            if (wind_speed <= self.cut_in_speed) | (wind_speed > self.cut_out_speed):
                y = 0.0

        return scale*y
        
    def get_thrust(self,wind_speed):
        if (wind_speed >= self.cut_in_speed) and (wind_speed < self.cut_out_speed):
            return np.interp(wind_speed,self.wind_speeds,self.thrust)

    def get_aep(self,dist,ref_height=None,alpha=None,units='MW',quad_kwargs={},dt=1.0,dist_type='weibull',wind_profile='power_law'):
        # dt is in hours. Default is 1.0 for hourly data.
    
        num_dt_per_year = int(8760/dt) # number of time steps per year

        if ref_height is None:
            print('No reference height provided. Assuming wind speed distribution is at the hub height.')
            integrand = lambda x: self.get_power(x,power_units=units)*dist.pdf(x)    
        else:
            if wind_profile.lower() == 'power_law':
                assert alpha is not None, 'Must provide alpha for power law model'
                print(f"Changing height from {ref_height} to {self.hub_height}m using power law with alpha = {alpha}")
                dist_at_hub_height = change_height(dist,ref_height,self.hub_height,alpha,dist_type)
                # integrand = lambda x: self.get_power(np.array([x*(self.hub_height/ref_height)**alpha]),power_units=units) * dist.pdf(x)
                integrand = lambda x: self.get_power(np.array([x]),power_units=units) * dist_at_hub_height.pdf(x)
            else:
                raise ValueError('Invalid height model. Only power_law implemented in this method so far.')

        res = quad(integrand,0,self.cut_out_speed,**quad_kwargs) # integrate from 0 to cut out speed
        return res[0]*dt*num_dt_per_year # AEP in MWh

    def export_to_sam_format(self,name):
        out = ",".join([name,str(self.rated_power/1000),str(self.rotor_diameter),'unknown'])
        
        wind_speeds = [f'{v}' for v in self.wind_speeds]
        wind_speeds = "|".join(wind_speeds)
        out = f"{out},{wind_speeds}"
        
        power = [f'{v/1000}' for v in self.power]
        power = "|".join(power)
        out = f"{out},{power}"

        return out

def initial_layout_rectangle(x_bnds:list,
                             y_bnds:list,
                             n_turbines:int,
                             n_grid_points:int = None,
                             rotation:float = 0.0,
                             placed_turbine_positions:list[np.ndarray]=None):
    
    rotation = np.deg2rad(rotation)
    
    # Farthest first over a bounding box
    if n_grid_points is None:
        n_grid_points = 100*n_turbines

    xg = np.linspace(x_bnds[0], x_bnds[1], int( np.floor(np.sqrt(n_grid_points))) )
    yg = np.linspace(y_bnds[0], y_bnds[1], int(np.ceil(np.sqrt(n_grid_points))))
    xg,yg = np.meshgrid(xg,yg)
    xg = xg.flatten()
    yg = yg.flatten()
    pg = np.c_[xg,yg]

    if placed_turbine_positions is None:
        idx = np.random.random_integers(0,pg.shape[0]-1,1)
        points = [np.r_[xg[idx],yg[idx]]]
    else:
        points = placed_turbine_positions

    for n in range(len(points),n_turbines):
        points,pg = _add_min_squared_distance(points,pg)
        
    # scale the rectangle to the bounding box
    if np.abs(rotation) > 0.0:
        corners = np.array([[x_bnds[0],y_bnds[0]],
                            [x_bnds[1],y_bnds[0]],
                            [x_bnds[0],y_bnds[1]],
                            [x_bnds[1],y_bnds[1]]])
        W = x_bnds[1] - x_bnds[0]
        H = y_bnds[1] - y_bnds[0]

        R = rotation_matrix(rotation)
        center = np.mean(corners,axis=0)
        rot_corners = np.dot(corners-center,R.T) + center
        M = np.max(np.abs(rot_corners-center),axis=0)
        scale = np.min([W/2/M[0],H/2/M[1]])
        
        points = scale*np.dot(points-center,R.T) + center

    return np.array(points)

def _add_min_squared_distance(points,grid):
    M = len(points)
    d2 = np.inf*np.ones((grid.shape[0],M))
    for ii,p in enumerate(points):
        d2[:,ii] = np.sum((grid - p)**2,axis=1)
    
    d2 = np.min(d2,axis=1)
    mask = np.argsort(d2)
    points.append(grid[mask[-1],:])
    grid = np.delete(grid,mask[-1],axis=0)
    
    return points,grid

def rotation_matrix(ϕ,homogeneous=False):
    
    if homogeneous:
        R = np.array([[np.cos(ϕ),-np.sin(ϕ),0],[np.sin(ϕ),np.cos(ϕ),0],
                      [0,0,1]])
    else:
        R = np.array([[np.cos(ϕ),-np.sin(ϕ)],[np.sin(ϕ),np.cos(ϕ)]])
    
    return R

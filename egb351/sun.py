import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
from mpl_toolkits.mplot3d import Axes3D
from datetime import datetime, timezone, timedelta
from .math import cosd, sind, acosd
from IPython.display import display, Markdown
from calendar import isleap
from ipywidgets import interact
import ipywidgets as widgets
from .constants import h,c,kB,σ
tol = np.finfo(np.float64).tiny

def get_day_of_year(d,ignore_leapday=False):
    doy = d.timetuple().tm_yday # day of year
    if ignore_leapday and isleap(d.year):
        if d.month > 2:
            doy -= 1
    return doy

def sun_solid_angle(R_sun,D_sun_earth):
    return 2*np.pi*(1-np.sqrt( (D_sun_earth**2-R_sun**2)/D_sun_earth**2 )) # [sr]

def black_body_spectral_radiance(λ,T,display_eqn=False):
    C1 = 2*h*c**2 * 1e30 * 1e-6 # 1e30 for m -> µm then 1e6 for W/m3 -> W/m2/µm
    C2 = h*c/kB * 1e6 # 1e6 for m -> µm 
    if display_eqn:
        display(Markdown( rf""" $E(\lambda,T) = \frac{{{2*h*c**2 * 1e30 * 1e-6:.4e}}}{{\lambda^5 \left[exp\left(\frac{{{C2:.4e}}}{{\lambda T}}\right) -1\right]}}\qquad [W/m^2 \, sr \, \mu m]$
                    """))
    return C1 / (λ**5 * (np.exp(C2/λ/T)-1) ) # [W/m2/µm/sr]

def black_body_spectral_exitance(λ,T):
    # assumes that the body is convex
    return np.pi()*black_body_spectral_radiance(λ,T) # [W/m2/µm]
    
def declination(doy):
    return 23.45*sind(360/365.0*(doy-81)) # δ = 23.45*sind(360/365.0*(284+doy))

def equation_of_time(doy):
    B = 360.0*(doy-81)/364.0
    E = 9.87*sind(2*B)-7.53*cosd(B)-1.5*sind(B)
    return E # minutes

def solar_time(t,L_tz,L_loc,ignore_leapday=False):
    doy = get_day_of_year(t,ignore_leapday=ignore_leapday)
    E = equation_of_time(doy)
    st = t + timedelta(minutes = 4*(L_tz-L_loc) + E)
    return st

def compute_solar_angles(standard_clock_time: datetime,φ,L,L_tz,force_south_as_zero=False,ignore_leapday=False):
    azimuth_zero = "South"
    doy = get_day_of_year(standard_clock_time,ignore_leapday=ignore_leapday)
    sdt = solar_time(standard_clock_time,L_tz,L,ignore_leapday=ignore_leapday)
    t_solar = sdt.hour + sdt.minute/60 + sdt.second/3600 + sdt.microsecond/1e6/3600
    ω = 15*(t_solar-12.0)
    δ = declination(doy)
    θ = acosd( sind(δ)*sind(φ)+cosd(δ)*cosd(φ)*cosd(ω) )
    γ = np.sign(ω+tol) * np.abs( acosd( (cosd(θ)*sind(φ) - sind(δ))/(sind(θ)*cosd(φ)) ) )

    if (φ < 0) and (not force_south_as_zero): # southern hemisphere and convention is not forced to use south as zero
        azimuth_zero = "North" # set azimuth zero to N
        γ = _set_azimuth_zero_to_north(γ)

    return {"zenith":θ, 
            "azimuth":γ,
            "declination":δ, 
            "solar_time": t_solar,
            "azimuth_zero":azimuth_zero}

def _set_azimuth_zero_to_north(γ):
    if γ<0: # east of north
        γ = -180 - γ 
    else: # west of north
        γ = 180 - γ
    
    return γ

def _all_sun_angles(dt0, hour_grid, doy_grid,φ,L,L_tz,force_south_as_zero=False,ignore_leapday=False,location_name=""):
    zeniths = np.zeros( (len(doy_grid),len(hour_grid)) )
    azimuths = np.zeros( (len(doy_grid),len(hour_grid)) )
    for ii,d in enumerate(doy_grid):
        for jj,h in enumerate(hour_grid):
            dt = dt0 + timedelta(days=d) + timedelta(hours=h)
            sun_angles = compute_solar_angles(dt,φ,L,L_tz,
                                              force_south_as_zero=force_south_as_zero,
                                              ignore_leapday=ignore_leapday)
            zeniths[ii,jj] = sun_angles['zenith']
            azimuths[ii,jj] = sun_angles['azimuth']
    
    return zeniths, azimuths, sun_angles['azimuth_zero']

def sun_path_diagram(dt0, hour_grid, doy_grid,φ,L,L_tz,force_south_as_zero=False,ignore_leapday=False,location_name=""):
    
    zeniths,azimuths,az_zero = _all_sun_angles(dt0, hour_grid, doy_grid,φ,L,L_tz,
                                       force_south_as_zero=force_south_as_zero,
                                       ignore_leapday=ignore_leapday,
                                       location_name=location_name)

    # Create a colormap
    cmap = plt.get_cmap("twilight")

    # Normalize the color values to be between 0 and 1
    norm = plt.Normalize(doy_grid.min(), doy_grid.max())

    fig,ax = plt.subplots()
    for ii,d in enumerate(doy_grid):
        ax.plot(azimuths[ii,:],90-zeniths[ii,:], color=cmap(norm(d)))
    ax.set_ylim((0,90))
    ax.set_xlabel(f"Azimuth Angle (degrees, + is West of {az_zero})")
    ax.set_ylabel("Altitude angle (degrees)")
    ax.set_title(f"Sun path for {location_name} ({φ:.2f},{L:.2f})")

    # Add a color bar to show the color scale
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)


    # Customize the colorbar labels to show date format
    def date_formatter(x, pos):
        return mdates.num2date(x).strftime('%b %d')
    # cbar.set_ticks([doy_grid.min(),121,242,doy_grid.max()])
    cbar.set_ticks(doy_grid)
    cbar.ax.yaxis.set_major_formatter(plt.FuncFormatter(date_formatter))
    
    return {    'figure_handle':fig, 
                'axis_handle':ax,
                'data': {'zenith':zeniths,
                         'azimuth':azimuths,
                         'azimuth_zero': az_zero}}

def sun_path_3d(dt0, hour_grid, doy_grid,φ,L,L_tz,
                     ignore_leapday=False,
                     location_name="",
                     convention='sam',
                     radius=1.0):
    
    θz,γs,az_zero = _all_sun_angles(dt0, hour_grid, doy_grid,φ,L,L_tz,
                                       force_south_as_zero=True,
                                       ignore_leapday=ignore_leapday,
                                       location_name=location_name)
    
    sun_vecs = np.zeros((len(doy_grid),len(hour_grid),3))
    for ii,_ in enumerate(doy_grid):
        for jj,_ in enumerate(hour_grid):
            sun_vecs[ii,jj,:] = solar_vector_from_angles(γs[ii,jj],θz[ii,jj],convention=convention)
    sun_vecs = radius*sun_vecs

    # doys = [f"Arc {i+1}" for i in range(len(doy_grid))]
    doys = [(dt0+timedelta(days=dt)).strftime('%d-%b') for dt in doy_grid]
    hods = [f"{h:.1f}" for h in hour_grid]
    def _plot_sun(doy,hod,angle=0,elev=30):
        di = doys.index(doy)
        hi = hods.index(hod)
        mag = 1.5*np.sum(sun_vecs**2,axis=2).max()

        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

        # sun arcs for each day
        for ii in range(sun_vecs.shape[0]):
            hours_above_zero = np.where(sun_vecs[ii,:,2]>=0)[0]
            ax.plot(sun_vecs[ii,hours_above_zero,0],
                    sun_vecs[ii,hours_above_zero,1],
                    sun_vecs[ii,hours_above_zero,2],
                    color='gray',
                    alpha=0.5,
                    linewidth=0.5)
    
        # highlight selected arc
        hours_above_zero = np.where(sun_vecs[di,:,2]>=0)[0]
        ax.plot(    sun_vecs[di,hours_above_zero,0],
                    sun_vecs[di,hours_above_zero,1],
                    sun_vecs[di,hours_above_zero,2],
                    color='black',
                    linewidth=0.5)
        
        # plot sun
        ax.scatter( sun_vecs[di,hi,0],
                    sun_vecs[di,hi,1],
                    sun_vecs[di,hi,2],
                    color='#FFDF00',
                    s=100)
        
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_zlabel('')

        # compass
        x_start,x_end,y_start,y_end,z_start,z_end = -mag,mag,0,0,0.0*mag,0.0*mag
        # ax.plot([x_start, x_end], [y_start, y_end],[0,0], marker=None)
        ax.text(x_start, y_start, z_start, "W", color='black', fontsize=12, ha='center')
        ax.text(x_end, y_end, z_end, "E", color='black', fontsize=12, ha='center')
        dx = x_end - x_start
        dy = y_end - y_start
        dz = z_end - z_start
        ax.quiver(x_start, y_start, z_start, dx, dy, dz, arrow_length_ratio=0.05, color='black')
        ax.quiver(x_end, y_end, z_end, -dx, -dy, -dz, arrow_length_ratio=0.05, color='black')

        x_start,x_end,y_start,y_end,z_start,z_end = 0,0,-mag,mag,0.0*mag,0.0*mag
        # ax.plot([x_start, x_end], [y_start, y_end],[0,0], marker=None)
        ax.text(x_start, y_start, z_start, "S", color='black', fontsize=12, ha='center')
        ax.text(x_end, y_end, z_end, "N", color='black', fontsize=12, ha='center')
        dx = x_end - x_start
        dy = y_end - y_start
        dz = z_end - z_start
        ax.quiver(x_start, y_start, z_start, dx, dy, dz, arrow_length_ratio=0.05, color='black')
        ax.quiver(x_end, y_end, z_end, -dx, -dy, -dz, arrow_length_ratio=0.05, color='black')

        ax.set_zlim((0,mag))
        ax.set_ylim((-mag,mag))
        ax.set_ylim((-mag,mag))
        ax.set_xticks(())
        ax.set_yticks(())
        ax.set_zticks(())

        ax.view_init(elev=elev, azim=angle)
        plt.draw()
    
    # Create sliders   
    arc_slider = widgets.SelectionSlider(options=doys, value=doys[0], description='Date:')
    point_slider = widgets.SelectionSlider(options=hods, value=hods[0], description='Hour:')
    angle_slider = widgets.FloatSlider(value=0, min=0, max=360, step=1, description='View Azimuth:')
    elev_slider = widgets.FloatSlider(value=30, min=0, max=90, step=1, description='View Elevation:')
    
    # interactive plot
    intp = interact(_plot_sun, 
                    doy=arc_slider, 
                    hod=point_slider,
                    angle=angle_slider,
                    elev=elev_slider)
    return intp

def solar_vector_from_angles(γ_s,θ_zenith,convention='sam'):
    g,z = γ_s,θ_zenith
    assert convention.lower() in ['soltrace','sam'], "Convention must be either soltrace or sam"

    if convention.lower() == 'soltrace':
        # x-west, y-zenith, z-north
        return np.array([sind(z)*sind(g),cosd(z),-sind(z)*cosd(g)])
    elif convention.lower() == 'sam':
        # x-east, y-north, z-zenith. See SAM manual 3D Shade Calculator for more on this
        return np.array([sind(z)*sind(-g),-sind(z)*cosd(-g),cosd(z)])

def angle_of_incidence(θ_z,γ_s,β,γ):
    return np.arccos(cosd(θ_z)*cosd(β) + sind(θ_z)*sind(β)*cosd(γ_s-γ))

def plane_of_array_irradiance(standard_clock_time,DNI,GHI,φ,L,L_tz,β,γ,
                              ρ=0,model="ISM",θ_z_max=90,
                              clip_negative=False,
                              ignore_leapday=False):
    
    sun_pos = compute_solar_angles(standard_clock_time,φ,L,L_tz,force_south_as_zero=True,ignore_leapday=ignore_leapday)
    θ_z,γ_s = sun_pos['zenith'],sun_pos['azimuth']
    cAOI = cosd(θ_z)*cosd(β) + sind(θ_z)*sind(β)*cosd(γ_s-γ)
    DHI = GHI-DNI*cosd(θ_z)

    if model.lower()=="ism":
        if θ_z < θ_z_max: # sun is above horizon
            G_POA = DNI*cAOI + DHI*0.5*(1+cosd(β)) + ρ*GHI*0.5*(1-cosd(β))
        else:
            G_POA = 0 
    else:
        raise ValueError("Model must be ISM for now.")

    if clip_negative:
        if G_POA < 0.0:
            G_POA = 0.0
    
    return G_POA,sun_pos
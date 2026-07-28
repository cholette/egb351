from scipy.integrate import RK45
from typing import Union
import csv
from numpy import log10, interp, sqrt, pi, linspace, diff, zeros, c_
from .math import cosd, acosd, sind, tand, atand
from . import sun
from CoolProp.CoolProp import PropsSI 
from .constants import σ,g

# References 
# [1]   Assael, M.J., Gialou, K., Kakosimos, K. et al. Thermal Conductivity of Reference Solid Materials. 
#       International Journal of Thermophysics 25, 397–408 (2004). 
#       https://doi.org/10.1023/B:IJOT.0000028477.74595.d5
#
# [2]   R. Forristall, “Heat Transfer Analysis and Modeling of a Parabolic Trough Solar Receiver Implemented 
#       in Engineering Equation Solver,” National Renewable Energy Lab. (NREL), Technical Report TP-550-34169.
#
# [3]   R. V. Padilla, G. Demirkaya, D. Y. Goswami, E. Stefanakos, and M. M. Rahman, “Heat transfer analysis 
#       of parabolic trough solar receiver,” Applied Energy, vol. 88, no. 12, pp. 5097–5110, Dec. 2011, 
#       doi: 10.1016/j.apenergy.2011.07.012.
#
# [4]   J. E. Braun and J. C. Mitchell, “Solar geometry for fixed and tracking surfaces,” 
#       Solar Energy, vol. 31, no. 5, pp. 439–444, Jan. 1983, doi: 10.1016/0038-092X(83)90046-4.
#
# [5]   W. B. Stein and M. Geyer, “Power From The Sun.” [Online]. 
#       Available: http://www.powerfromthesun.net/book.html

def incidence_angle(θz,γs,β,γ):
    θ = acosd( cosd(θz)*cosd(β) + sind(θz)*sind(β)*cosd(γs-γ) )    
    return θ

def incidence_angle2(θz,γs,γ_axis):
    # Incidence angle using rotation convention in [5] to tilt + azimuth
    α = 90-θz

    # convert input azimuths from west +ve -> east +ve with a 0 to 360 range
    if γ_axis>0:
        γ_axis = 360-γ_axis
    else:
        γ_axis = -γ_axis

    if γs>0:
        γs=360-γs
    else:    
        γs = -γs

    return acosd( sqrt(1-(cosd(α)*cosd(γs-γ_axis))**2 ))

def one_axis_tracking_angles(θz,γs,γ_axis,output='tilt_azimuth'):
    if output.lower() == 'tilt_azimuth':
        # Method from [4]. Output is tilt and azimuth of aperture. 
        if γs-γ_axis >= 0:
            γ = γ_axis + 90
        else:
            γ = γ_axis - 90

        β=atand(tand(θz)*cosd(γs-γ))
        if β<0:
            β+=180
        return β,γ
    
    elif output.lower() == 'rotation_about_axis_azimuth':
        # Method from [5]. 
        α = 90-θz

        # Convert to [5] azimuth convention
        if γ_axis>0:
            γ_pfs = 360-γ_axis
        else:
            γ_pfs = -γ_axis

        if γs>0:
            γs_pfs=360-γs
        else:    
            γs_pfs = -γs

        rotation = atand(sind(γs_pfs-γ_pfs)/tand(α))
        return rotation      
class parabolic_trough_1D_transient:
    def __init__(self):              
        self.envelope = {'material':None,
                        'diameter_inner':None,
                        'diameter_outer':None,
                        'volumetric_heat_capacity':None,
                        'emissivity':None,
                        'conductivity':None,
                        'transmittance':None,
                        'absorptance':None,
                        'eta_opt':None
                        }
        
        self.absorber = {'material':None,
                        'diameter_inner':None,
                        'diameter_outer':None,
                        'volumetric_heat_capacity':None,
                        'emissivity':None,
                        'conductivity':None,
                        'transmittance':None,
                        'absorptance':None,
                        }
        
        self.htf = {'fluid_name':None}

        self.collector = {  'hce_length':None,
                            'aperture_width':None,
                            'number_of_hce':None,
                            'epsilon_shadowing':None,
                            'epsilon_tracking':None,
                            'epsilon_geometry':None,
                            'rho_nominal':None,
                            'epsilon_other':None
                          }
        
        self.bracket = {'conductivity':None,
                        'perimeter': None,
                        'min_cross_sectional_area': None,
                        'effective_diameter': None
                        }
                          
    def import_material_parameters(self,csv_file):
        with open(csv_file) as file:
            reader = csv.reader(file)
            headers = next(reader) # skip header row
            category_col = [s.lower()=='category' for s in headers].index(True) 
            value_col = [s.lower()=='value' for s in headers].index(True)
            name_col = [s.lower()=='parameter' for s in headers].index(True)

            # populate property names
            components = []
            for c in dir(self):
                if not c.startswith('__'):
                    if not callable(getattr(self,c)):
                        components.append(c) 

            properties = {c:list(getattr(self,c).keys()) for c in components}
            is_prop = lambda p,c: p in properties[c]

            for row in reader:
                comp_name,prop_name,prop_val = row[category_col].lower(),row[name_col].lower(),row[value_col]

                # write the properties if they exist
                if comp_name in properties.keys():
                    if is_prop(prop_name,comp_name):
                        a = getattr(self,comp_name)
                        try:
                             prop_val = float(prop_val)
                        except:
                            print(f"Property value {prop_val} for {prop_name} is not numeric. Writing raw string.")

                        a[prop_name]=prop_val
                        setattr(self,comp_name,a)
                    elif prop_name[-9::].lower() == "_lookup_x":
                        # Replace a property with a lookup function if the right tags are set in the .csv file.
                        # Looks for [property_name]_lookup_fx and writes a function handle to the property. 
                        # Right now, the only allowable "x" is temperature. 
                        if is_prop(prop_name[0:-9].lower(),comp_name):
                            next_row = next(reader)
                            comp = getattr(self,comp_name)
                            X = [float(v) for v in row[value_col].split(';')]
                            Y = [float(v) for v in next_row[value_col].split(';')]
                            comp[prop_name[0:-9].lower()] = lambda T,X=X,Y=Y: interp(T,X,Y) # weird syntax to address overwriting of lambda functions [here](https://stackoverflow.com/questions/3431676/creating-functions-or-lambdas-in-a-loop-or-comprehension)
                        else:
                            raise ValueError(f"No parameter {prop_name[0:-9]} in {comp_name}.")
                elif comp_name not in ['','reference']:
                    raise ValueError(f"Parameter category {comp_name} not recognized")
            
    def absorber_to_htf_convection(self,T,m_dot,pressure=1e7):
        T_htf,T_absorber,T_envelope = T[0],T[1],T[2]
        htf = self.htf['fluid_name']
        D = self.absorber['diameter_inner']
        Aai = pi/4*D**2

        # temperature dependent properties
        k = PropsSI('conductivity','T',T_htf,'P',pressure,htf)
        Pr = PropsSI('Prandtl','T',T_htf,'P',pressure,htf)
        Prw = PropsSI('Prandtl','T',T_absorber,'P',pressure,htf)
        mu = PropsSI('viscosity','T',T_htf,'P',pressure,htf)
        ReD = (m_dot/Aai)*D/mu
        if ReD > 5e6:
            print('Warning: Reynolds number if outside the range of validity.')
        if ReD < 2300: # laminar
            Nu = 4.36
        else:
            f2 = (1.82*log10(ReD)-1.64)**(-2.0)
            Nu = f2/8.0 * (ReD-1000)*Pr / (1+12.7*sqrt(f2/8.0) * (Pr**(2.0/3.0) - 1)) *(Pr/Prw)**0.11
        
        return pi*Nu*k*(T_absorber - T_htf)

    def absorber_to_envelope_radiation(self,T):
        # from [2], but [3] is probably more accurate
        T_htf,T_absorber,T_envelope = T
        Doa = self.absorber['diameter_outer']
        Die = self.envelope['diameter_inner']
        εa = self.absorber['emissivity'] # absorber coating emissivity
        εe = self.envelope['emissivity']

        return σ*pi*Doa*(T_absorber**4 - T_envelope**4) / (1/εa + (1-εe)*Doa/(εe*Die))

    def nusselt_number(self,T_object,T_fluid,wind_speed,D=None,air_pressure=101e3):
        
        T_film = 0.5*(T_object+T_fluid)
        if D is None:
            D = self.envelope['diameter_outer']
        
        mu = PropsSI('viscosity','T',T_film,'P',air_pressure,'Air')
        rho = PropsSI('D','T',T_film,'P',air_pressure,'Air')
        if wind_speed<0.1: # m/s
            Pr_film = PropsSI('Prandtl','T',T_film,'P',air_pressure,'Air')
            β = PropsSI('isobaric_expansion_coefficient','T',T_film,'P',air_pressure,'Air')
            k = PropsSI('conductivity','T',T_film,'P',air_pressure,'Air')
            cp = PropsSI('C','T',T_film,'P',air_pressure,'Air')
            α = k/rho/cp
            v = mu/rho            
            
            RaD = g*β*D**3*(T_object-T_fluid)/(α*v)
            den = (1+ (0.559/Pr_film)**(9/16))**(16/9)
            Nu = (0.6 + 0.387*(RaD/den)**(1/6) )**2
        else:
            Pr_e = PropsSI('Prandtl','T',T_object,'P',air_pressure,'Air')
            Pr_ambient = PropsSI('Prandtl','T',T_fluid,'P',air_pressure,'Air')
            ReD = rho * wind_speed * D / mu

            if (Pr_ambient < 0.7) or (Pr_ambient>500):
                print('Warning: Prandtl number if outside the range of validity in the Nu computation.')
            if (ReD < 1.0) or (ReD > 1e6):
                print('Warning: Reynolds number if outside the range of validity in the Nu computation.')

            if ReD <= 40:
                C,m = 0.75,0.4
            elif ReD <= 1000:
                C,m = 0.51,0.5
            elif ReD<= 200e3:
                C,m = 0.26,0.6
            else:
                C,m = 0.076,0.7

            if Pr_ambient<=10:
                n=0.37
            else:
                n=0.36
            
            Nu = C*ReD**m * Pr_ambient**n * (Pr_ambient/Pr_e)**(1/4.0)
        
        return Nu          
        
    def envelope_to_surroundings_convection(self,T,T_ambient,wind_speed,air_pressure=101e3):
         T_htf,T_absorber,T_envelope = T[0],T[1],T[2]
         T_film = 0.5*(T_envelope+T_ambient)
         k = PropsSI('conductivity','T',T_film,'P',air_pressure,'Air')
         D = self.envelope['diameter_outer']
         Nu = self.nusselt_number(T_envelope,T_ambient,wind_speed,D=D)
         he = k*Nu/D
         return pi*he*D*(T_envelope-T_ambient)

    def envelope_to_surroundings_radiation(self,T,T_sky):
        T_htf,T_absorber,T_envelope = T[0],T[1],T[2]
        Doe = self.envelope['diameter_outer']
        εe = self.envelope['emissivity']
        return σ*Doe*pi*εe*(T_envelope**4 - T_sky**4)

    def bracket_conduction(self,n_brackets,T_ambient,T_base,wind_speed,air_pressure=101e3):
        
        Pb = self.bracket['perimeter']
        Ab = self.bracket['min_cross_sectional_area']
        L_hce = self.collector['hce_length']
        kb = self.bracket['conductivity']
        D = self.bracket['effective_diameter']
        T_bracket = (T_base-273.16 + T_ambient-273.15)/3.0 + 273.15
        T_film = 0.5*(T_bracket+T_ambient)
        k = PropsSI('conductivity','T',T_film,'P',air_pressure,'Air')

        Nu = self.nusselt_number(T_bracket,T_ambient,wind_speed,D=D)
        hb = k*Nu/D

        return n_brackets*sqrt(hb*Pb*kb*Ab) * (T_base-T_ambient)

    def solar_heating_rates(self,DNI,η_opt):
        
        W = self.collector['aperture_width']
        α_envelope = self.envelope['absorptance']
        α_absorber = self.absorber['absorptance']
        τ_envelope = self.envelope['transmittance']

        In = DNI*W
        Qe = η_opt*α_envelope*In
        Qa = η_opt*τ_envelope*α_absorber*In

        return Qa,Qe

    def _apply_boundary_conditions(self,T):

        # enforce boundary conditions
        T[0,1]  = T[1,1]    # absorber
        T[-1,1] = T[-2,1]
        
        T[0,2]  = T[1,2]    # envelope
        T[-1,2] = T[-2,2]
    
        return T

    def pt_odes(self,T,ambient,massflow,N_elements,n_brackets,ΔT_base=-8,):

        T = T.reshape(N_elements,3)
        T = self._apply_boundary_conditions(T)

        T_ambient = ambient['T_ambient']
        T_sky = ambient['T_sky']
        wind_speed = ambient['wind_speed'] 
        DNI = ambient['DNI']
        η_opt = ambient['η_opt']
        
        # partition
        N_hce = self.collector['number_of_hce']
        L_hce = self.collector['hce_length']
        L_aperture = L_hce*N_hce
        element_edges = linspace(0,L_hce*N_hce,N_elements+1)
        element_centers = element_edges[0:-1] + 0.5*diff(element_edges)
        Δz = element_edges[1]-element_edges[0]

        # heat flows
        Q_solar = [self.solar_heating_rates(DNI,η_opt)]*N_elements
        Qa = [q[0] for q in Q_solar]
        Qe = [q[1] for q in Q_solar]
        Qaf_conv = [self.absorber_to_htf_convection(t,massflow) for t in T]
        Qae_rad = [self.absorber_to_envelope_radiation(t) for t in T]
        Qcond_bracket = [self.bracket_conduction(n_brackets,T_ambient,t[1]+ΔT_base,wind_speed) for t in T]
        Qesa_conv = [self.envelope_to_surroundings_convection(t,T_ambient,wind_speed) for t in T]
        Qes_rad = [self.envelope_to_surroundings_radiation(t,T_sky) for t in T]

        # HTF
        htf = self.htf['fluid_name']
        Aai = pi/4*self.absorber['diameter_inner']**2
        Cpf = [PropsSI('C','T',t[0],'P',1e7,htf) for t in T]
        rhof = [PropsSI('D','T',t[0],'P',1e7,htf) for t in T]
        Vf = [massflow/Aai/rhof[ii] for ii in range(N_elements)]
        dTf = zeros(N_elements)
        for ii in range(1,N_elements):
            K = 1.0/Aai/rhof[ii]/Cpf[ii]
            dTf[ii] = -K*massflow/Δz*(Cpf[ii]*T[ii][0]+Vf[ii]**2/2 - 
                         Cpf[ii-1]*T[ii-1][0] - Vf[ii-1]**2/2) + K*Qaf_conv[ii]
            
        # Absorber
        Aa = pi/4*(self.absorber['diameter_outer']**2 - self.absorber['diameter_inner']**2)
        rho_Cp_a = self.absorber['volumetric_heat_capacity']
        ka = self.absorber['conductivity']
        dTa = zeros(N_elements)
        for ii in range(1,N_elements-1):
            Tim1,Ti,Tip1 = T[ii-1][1],T[ii][1],T[ii+1][1]
            ΣQ = Qa[ii] - Qaf_conv[ii] - Qae_rad[ii] - Qcond_bracket[ii]
            ka0 = 0.5*(ka(Tim1) + ka(Ti))
            ka1 = 0.5*(ka(Ti) + ka(Tip1))
            rcp = rho_Cp_a(Ti)*1000
            dTa[ii] = (ka1*Tip1 - ka1*Ti - ka0*Ti + ka0*Tim1)/(Δz**2 * rcp) + 1.0/Aa/rcp * ΣQ

        # Envelope
        Ae = pi/4*(self.envelope['diameter_outer']**2 - self.envelope['diameter_inner']**2)
        rho_Cp_e = self.envelope['volumetric_heat_capacity']
        ke = self.envelope['conductivity']
        dTe = zeros(N_elements)
        for ii in range(1,N_elements-1):
            Tim1,Ti,Tip1 = T[ii-1][2],T[ii][2],T[ii+1][2]
            rcp = rho_Cp_e(Ti)*1000
            ΣQe = Qe[ii] + Qae_rad[ii] - Qesa_conv[ii] - Qes_rad[ii]
            Γip = (ke(Tip1)+ke(Ti))/(2*Δz**2*rcp)
            Γim = (ke(Ti)+ke(Tim1))/(2*Δz**2*rcp)
            Γi = Γip + Γim
            dTe[ii] = Γip*Tip1 - Γi*Ti + Γim*Tim1 + ΣQe/Ae/rcp

        return c_[dTf,dTa,dTe].flatten()
        

        













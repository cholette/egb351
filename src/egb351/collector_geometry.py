import matplotlib.pyplot as plt
import numpy as np
from .math import cosd, sind, acosd
import pyvista as pviz
import pyviewfactor as pvf
from IPython.display import display, Markdown
from tqdm import tqdm

def flat_panel_unit_vectors(β,γ):
    # Computes the panel normal vector expressed in the global reference frame. 
    # SAM convention used for global system: x-east, y-north, z-zenith. See SAM manual 3D Shade Calculator for more on this
    # Panel starts out aligned with the global system and undergoes a z rotation and then a y rotation (body-fixed)
    # Azimuth rotation follows convention in Duffie and Beckman (+ve is west of south)
    
    # rotate -γ about z axis (since +ve rotation about z corresponds to an eastward rotation, the o)
    Rz = np.array([[cosd(-γ),sind(-γ),0],
                   [-sind(-γ),cosd(-γ),0],
                   [0,0,1]])

    # rotate β about transformed x axis
    Rx = np.array([[1,0,0],
                   [0,cosd(β),sind(β)],
                   [0,-sind(β),cosd(β)]])
    
    R_total = Rx@Rz # XYZ -> x'y'z' 
    R = R_total.T   # x'y'z' -> XYZ
    u = R@np.array([1,0,0])
    v = R@np.array([0,1,0])
    n = R@np.array([0,0,1])
    return n,u,v

def face(d,i):
    return pvf.fc_unstruc2poly(d.extract_cells(i))

def sky_view_factor(panel,obstruction,sky):
    F = np.zeros((sky.n_cells,panel.n_cells))
    for ii in tqdm(range(sky.n_cells),total=sky.n_cells):
        for jj in range(panel.n_cells):
            sky_face = face(sky,ii)
            panel_face = face(panel,jj)

            if isinstance(obstruction,pviz.core.composite.MultiBlock):
                visible = True
                for o in obstruction:
                    v = pvf.get_visibility_raytrace(sky_face,panel_face,o)
                    if not v:
                        visible = False
                        break

            elif isinstance(obstruction,pviz.core.pointset.PolyData):
                visible = pvf.get_visibility_raytrace(sky_face, panel_face, obstruction)
            else:
                raise TypeError("Obsctuction needs to be either a pyvista PolyData or MultiBlock")
            
            if visible: # ... if no obstruction
                F[ii,jj] = pvf.compute_viewfactor(sky_face,panel_face)

    return F
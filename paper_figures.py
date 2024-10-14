########################## Importing libraries
import intake
import cmocean
import math
import cartopy.crs as ccrs
import cartopy.feature as cf
import numpy as np
import xarray as xr
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
import matplotlib.ticker as ticker
from matplotlib.cm import get_cmap
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
from matplotlib.ticker import MultipleLocator, NullFormatter
from matplotlib.gridspec import GridSpec

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams.update({'font.size': 24})



# To plot figure 1


def worldmap_fig1(var1, var2, var3, var4, var5, var6, lon, lat):
    projection = ccrs.PlateCarree()
    fig, ax = plt.subplots(nrows=3, ncols=2, figsize=(24, 29), subplot_kw={"projection": projection})  # Increased vertical size by 50%

    cmap = cm.coolwarm

    # Define the levels for the colorbar
    levels = [-32, -16, -8, -4, -2, -1, -0.5, 0 , 0.5, 1 , 2, 4, 8, 16, 32]

    # Create a norm and colorbar using BoundaryNorm
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    for a in ax.flat:
        a.set_global()
        a.add_feature(cf.COASTLINE, linewidth=0.8)
        a.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
        a.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
        lon_formatter = LongitudeFormatter()
        lat_formatter = LatitudeFormatter()
        a.xaxis.set_major_formatter(lon_formatter)
        a.yaxis.set_major_formatter(lat_formatter)
        a.minorticks_on()
        a.tick_params(which='major', direction='in', width=3.5)
        a.tick_params(which='minor', direction='in', width=2)
        a.set_aspect(aspect=1.25)  # Set the aspect ratio to make the subplots taller

    # Remove xticks and yticks for the middle figures
    ax[0, 0].set_xticks([])
    ax[0, 1].set_xticks([])
    ax[1, 0].set_xticks([])
    ax[1, 1].set_xticks([])
    ax[0, 1].set_yticks([])
    ax[1, 1].set_yticks([])
    ax[2, 1].set_yticks([])

    ax[0, 0].xaxis.set_minor_locator(MultipleLocator(30))
    ax[0, 0].xaxis.set_minor_formatter(NullFormatter())
    
    ax[0, 1].xaxis.set_minor_locator(MultipleLocator(30))
    ax[0, 1].xaxis.set_minor_formatter(NullFormatter())
        
   
    
    ax[1, 0].xaxis.set_minor_locator(MultipleLocator(30))
    ax[1, 0].xaxis.set_minor_formatter(NullFormatter())
    
    ax[1, 1].xaxis.set_minor_locator(MultipleLocator(30))
    ax[1, 1].xaxis.set_minor_formatter(NullFormatter())
  
    ax[0, 1].yaxis.set_minor_locator(MultipleLocator(30))
    ax[0, 1].yaxis.set_minor_formatter(NullFormatter())
  
    
    ax[1, 1].yaxis.set_minor_locator(MultipleLocator(30))
    ax[1, 1].yaxis.set_minor_formatter(NullFormatter())
 
    ax[2, 1].yaxis.set_minor_locator(MultipleLocator(30))
    ax[2, 1].yaxis.set_minor_formatter(NullFormatter())
 

    mapp1 = ax[0, 0].contourf(lon, lat, var1, levels=levels, cmap=cmap, norm=norm, extend="both", transform=ccrs.PlateCarree())
    ax[0, 0].set_title(r'(a) $z=8$ km; ICON', fontsize=25)

    mapp2 = ax[0, 1].contourf(lon, lat, var2, levels=levels, cmap=cmap, norm=norm, extend="both", transform=ccrs.PlateCarree())
    ax[0, 1].set_title(r'(b) $z=8$ km; Parameterization', fontsize=25)

    mapp3 = ax[1, 0].contourf(lon, lat, var3, levels=levels, cmap=cmap, norm=norm, extend="both", transform=ccrs.PlateCarree())
    ax[1, 0].set_title(r'(c) $z=16$ km; ICON', fontsize=25)

    mapp4 = ax[1, 1].contourf(lon, lat, var4, levels=levels, cmap=cmap, norm=norm, extend="both", transform=ccrs.PlateCarree())
    ax[1, 1].set_title(r'(d) $z=16$ km; Parameterization', fontsize=25)

 
    mapp5 = ax[2, 0].contourf(lon, lat, var5, levels=levels, cmap=cmap, norm=norm, extend="both", transform=ccrs.PlateCarree())
    ax[2, 0].set_title(r'(e) $z=23$ km; ICON', fontsize=25)

    mapp6 = ax[2, 1].contourf(lon, lat, var6, levels=levels, cmap=cmap, norm=norm, extend="both", transform=ccrs.PlateCarree())
    ax[2, 1].set_title(r'(f) $z=23$ km; Parameterization', fontsize=25)

 
    # Remove the gap between subplots
    plt.subplots_adjust(hspace=0.15, wspace=0)

    # Create a single colorbar at the bottom of the figure
    cbar = fig.colorbar(mapp6, ax=ax[:, :], orientation="horizontal", shrink=0.9, pad=0.05, aspect=40)
    cbar.set_ticks(levels)
    # Set the label below the colorbar
    cbar.set_label(r'$F^z$ (mPa)', fontsize=25)


# To plot figure 2


def worldmap_fig2(var1, var2, var3, var4, var5, var6, lon, lat):
    projection = ccrs.PlateCarree()
    fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(30, 24), subplot_kw={"projection": projection})  # Adjusted figsize

    cmap = cm.coolwarm

    # Define the levels for the colorbar
    levels = [-32, -16, -8, -4, -2, -1, -0.5, 0 , 0.5, 1 , 2, 4, 8, 16, 32]

    # Create a norm and colorbar using BoundaryNorm
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    for a in ax.flat:
        a.set_global()
        a.add_feature(cf.COASTLINE, linewidth=0.8)
        a.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
        a.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
        lon_formatter = LongitudeFormatter()
        lat_formatter = LatitudeFormatter()
        a.xaxis.set_major_formatter(lon_formatter)
        a.yaxis.set_major_formatter(lat_formatter)
        a.minorticks_on()
        a.tick_params(which='major', direction='in', width=3.5)
        a.tick_params(which='minor', direction='in', width=2)
        a.set_aspect(aspect=1.25)  # Set the aspect ratio to make the subplots taller

    # Remove xticks and yticks for the middle figures
    for i in range(3):
        if i != 0:
            ax[0, i].set_yticks([])
        ax[1, i].set_yticks([])

        ax[0, i].xaxis.set_minor_locator(MultipleLocator(30))
        ax[0, i].xaxis.set_minor_formatter(NullFormatter())

        ax[1, i].xaxis.set_minor_locator(MultipleLocator(30))
        ax[1, i].xaxis.set_minor_formatter(NullFormatter())

    mapp1 = ax[0, 0].contourf(lon, lat, var1, levels=levels, cmap=cmap, norm=norm, extend="both", transform=ccrs.PlateCarree())
    ax[0, 0].set_title(r'(a) $z=16$ km; Orog GWs', fontsize=25)

    mapp2 = ax[0, 1].contourf(lon, lat, var2, levels=levels, cmap=cmap, norm=norm, extend="both", transform=ccrs.PlateCarree())
    ax[0, 1].set_title(r'(b) $z=16$ km; Convective GWs', fontsize=25)

    mapp3 = ax[0, 2].contourf(lon, lat, var3, levels=levels, cmap=cmap, norm=norm, extend="both", transform=ccrs.PlateCarree())
    ax[0, 2].set_title(r'(c) $z=16$ km; Frontal GWs', fontsize=25)

    mapp4 = ax[1, 0].contourf(lon, lat, var4, levels=levels, cmap=cmap, norm=norm, extend="both", transform=ccrs.PlateCarree())
    ax[1, 0].set_title(r'(d) $z=23$ km; Orog GWs', fontsize=25)

    mapp5 = ax[1, 1].contourf(lon, lat, var5, levels=levels, cmap=cmap, norm=norm, extend="both", transform=ccrs.PlateCarree())
    ax[1, 1].set_title(r'(e) $z=23$ km; Convective GWs', fontsize=25)

    mapp6 = ax[1, 2].contourf(lon, lat, var6, levels=levels, cmap=cmap, norm=norm, extend="both", transform=ccrs.PlateCarree())
    ax[1, 2].set_title(r'(f) $z=23$ km; Frontal GWs', fontsize=25)

    # Remove the gap between subplots
    plt.subplots_adjust(hspace=0.3, wspace=0)

    # Create a single colorbar at the bottom of the figure
    cbar = fig.colorbar(mapp6, ax=ax[:], orientation="horizontal", shrink=0.9, pad=0.05, aspect=40)
    cbar.set_ticks(levels)
    # Set the label below the colorbar
    cbar.set_label(r'$F^z$ (mPa)', fontsize=25)
    

# To plot figure 3

def worldmap_fig3(var1, var2, var3, var4, var5, var6, var7, var8,lon, lat):
    projection = ccrs.PlateCarree()
    fig, ax = plt.subplots(nrows=4, ncols=2, figsize=(24, 37), subplot_kw={"projection": projection})  # Increased vertical size by 50%

    cmap = cm.coolwarm

    # Define the levels for the colorbar
    levels = [-32, -16, -8, -4, -2, -1, -0.5, 0 , 0.5, 1 , 2, 4, 8, 16, 32]

    # Create a norm and colorbar using BoundaryNorm
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    for a in ax.flat:
        a.set_global()
        a.add_feature(cf.COASTLINE, linewidth=0.8)
        a.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
        a.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
        lon_formatter = LongitudeFormatter()
        lat_formatter = LatitudeFormatter()
        a.xaxis.set_major_formatter(lon_formatter)
        a.yaxis.set_major_formatter(lat_formatter)
        a.minorticks_on()
        a.tick_params(which='major', direction='in', width=3.5)
        a.tick_params(which='minor', direction='in', width=2)
        a.set_aspect(aspect=1.25)  # Set the aspect ratio to make the subplots taller

    # Remove xticks and yticks for the middle figures
    ax[0, 0].set_xticks([])
    ax[0, 1].set_xticks([])
    ax[1, 0].set_xticks([])
    ax[1, 1].set_xticks([])
    ax[2, 0].set_xticks([])
    ax[2, 1].set_xticks([])   


    ax[0, 0].xaxis.set_minor_locator(MultipleLocator(30))
    ax[0, 0].xaxis.set_minor_formatter(NullFormatter())
    
    ax[0, 1].xaxis.set_minor_locator(MultipleLocator(30))
    ax[0, 1].xaxis.set_minor_formatter(NullFormatter())
        
 
    ax[1, 0].xaxis.set_minor_locator(MultipleLocator(30))
    ax[1, 0].xaxis.set_minor_formatter(NullFormatter())
    
    ax[1, 1].xaxis.set_minor_locator(MultipleLocator(30))
    ax[1, 1].xaxis.set_minor_formatter(NullFormatter())
    
    ax[2, 0].xaxis.set_minor_locator(MultipleLocator(30))
    ax[2, 0].xaxis.set_minor_formatter(NullFormatter())
    
    ax[2, 1].xaxis.set_minor_locator(MultipleLocator(30))
    ax[2, 1].xaxis.set_minor_formatter(NullFormatter())
  


    ax[0, 1].set_yticks([])
    ax[1, 1].set_yticks([])
    ax[2, 1].set_yticks([])
    ax[3, 1].set_yticks([])
    

    ax[0, 1].yaxis.set_minor_locator(MultipleLocator(30))
    ax[0, 1].yaxis.set_minor_formatter(NullFormatter())
 
    ax[1, 1].yaxis.set_minor_locator(MultipleLocator(30))
    ax[1, 1].yaxis.set_minor_formatter(NullFormatter())
    
    ax[2, 1].yaxis.set_minor_locator(MultipleLocator(30))
    ax[2, 1].yaxis.set_minor_formatter(NullFormatter())

    ax[3, 1].yaxis.set_minor_locator(MultipleLocator(30))
    ax[3, 1].yaxis.set_minor_formatter(NullFormatter())
   
    
    mapp1 = ax[0, 0].contourf(lon, lat, var1, levels=levels, cmap=cmap, norm=norm, extend="both", transform=ccrs.PlateCarree())
    ax[0, 0].set_title(r'(a) $z=16$ km; eastward; ICON', fontsize=20)

    mapp2 = ax[0, 1].contourf(lon, lat, var2, levels=levels, cmap=cmap, norm=norm, extend="both", transform=ccrs.PlateCarree())
    ax[0, 1].set_title(r'(b) $z=16$ km; eastward; Parameterization', fontsize=20)

    mapp3 = ax[2, 0].contourf(lon, lat, var3, levels=levels, cmap=cmap, norm=norm, extend="both", transform=ccrs.PlateCarree())
    ax[2, 0].set_title(r'(e) $z=16$ km; westward; ICON', fontsize=20)

    mapp4 = ax[2, 1].contourf(lon, lat, var4, levels=levels, cmap=cmap, norm=norm, extend="both", transform=ccrs.PlateCarree())
    ax[2, 1].set_title(r'(f) $z=16$ km; westward; Parameterization', fontsize=20)

    mapp5 = ax[1, 0].contourf(lon, lat, var5, levels=levels, cmap=cmap, norm=norm, extend="both", transform=ccrs.PlateCarree())
    ax[1, 0].set_title(r'(c) $z=23$ km; eastward; ICON', fontsize=20)

    mapp6 = ax[1, 1].contourf(lon, lat, var6, levels=levels, cmap=cmap, norm=norm, extend="both", transform=ccrs.PlateCarree())
    ax[1, 1].set_title(r'(d) $z=23$ km; eastward; Parameterization', fontsize=20)
    
    mapp7 = ax[3, 0].contourf(lon, lat, var7, levels=levels, cmap=cmap, norm=norm, extend="both", transform=ccrs.PlateCarree())
    ax[3, 0].set_title(r'(g) $z=23$ km; westward; ICON', fontsize=20)

    mapp8 = ax[3, 1].contourf(lon, lat, var8, levels=levels, cmap=cmap, norm=norm, extend="both", transform=ccrs.PlateCarree())
    ax[3, 1].set_title(r'(h) $z=23$ km; westward; Parameterization', fontsize=20)

  
    # Remove the gap between subplots
    plt.subplots_adjust(hspace=0.1, wspace=0)

    # Create a single colorbar at the bottom of the figure
    cbar = fig.colorbar(mapp8, ax=ax[:, :], orientation="horizontal", shrink=0.9, pad=0.05, aspect=40)
    cbar.set_ticks(levels)
    # Set the label below the colorbar
    cbar.set_label(r'$F^z$ (mPa)', fontsize=25)
    
    
    
    
######################################### Importing results stored in .nc




ds88 = xr.open_dataset('/.../gwd_zonal_lonlat_QBOi_202003_167it.nc')  # define your the location of the stored data
ds88p = xr.open_dataset('/.../gwd_zonal_lonlat_p_QBOi_202003_167it.nc')  # define your the location of the stored data
ds88n = xr.open_dataset('/.../gwd_zonal_lonlat_n_QBOi_202003_167it.nc')  # define your the location of the stored data
ds88g = xr.open_dataset('/.../gwd_healpix_lonlat_QBOi_202003_167it.nc')  # define your the location of the stored data


# longitude and latitude from stored data

lon = ds88.longitude
lat = ds88.latitude


################################ Get DATA from ICON and Parameterization

ICON = 1000*ds88.rupwp[:,:,:,:] # Net fluxes
ICON=ICON.mean(axis=0)  # temopral average

ICONp = 1000*ds88p.rupwp[:,:,:,:] # Eastward fluxes
ICONp=ICONp.mean(axis=0)  # temopral average


ICONn = 1000*ds88n.rupwp[:,:,:,:] # Westward fluxes
ICONn=ICONn.mean(axis=0)  # temopral average

# Get Param fluxes and do a temporal average on data

GW_or  = 1000*(ds88.ustror[:,:,:,:]) # Orographic GW fluxes
GW_or=GW_or.mean(axis=0)  # temopral average
GW_pr  = 1000*(ds88.ustrpr[:,:,:,:]) # Non-Orographic GW fluxes due to precipitations
GW_pr=GW_pr.mean(axis=0)  # temopral average
GW_fr  = 1000*(ds88.ustrfr[:,:,:,:]) # Non-Orographic GW fluxes due to fronts
GW_fr=GW_fr.mean(axis=0)  # temopral average

GW  = 1000*(ds88.ustror[:,:,:,:]+ds88.ustrpr[:,:,:,:]+ds88.ustrfr[:,:,:,:])    # Net GW fluxes       
GW=GW.mean(axis=0)  # temopral average
GWp  = 1000*(ds88p.ustror[:,:,:,:]+ds88p.ustrpr[:,:,:,:]+ds88p.ustrfr[:,:,:,:])  # Eastward GW fluxes         
GWp=GWp.mean(axis=0)  # temopral average
GWn  = 1000*(ds88n.ustror[:,:,:,:]+ds88n.ustrpr[:,:,:,:]+ds88n.ustrfr[:,:,:,:])  # Westward GW fluxes          
GWn=GWn.mean(axis=0)  # temopral average

vitu1 = ds88g.vitu[:,:,:,:] # zonal wind 
vitu1 = vitu1.mean(axis=0) # temopral average
vitu1 = vitu1.mean(axis=2) # zonal average

dICON =-ds88g.drupwpdz[:,:,:,:]  # ICON: Drag from stored data
dICON=dICON.mean(axis=0)  # temopral average
dICON=dICON.mean(axis=2)  # zonal average

dEPFphi =-ds88g.dEPFdphi[:,:,:,:] # ICON: horizontal component of the Eliasen Palm flux Drag from stored data
dEPFphi=dEPFphi.mean(axis=0)
dEPFphi=dEPFphi.mean(axis=2)

dParam =-( ds88g.du_orstressdz[:,:,:,:]+ds88g.du_prstressdz[:,:,:,:]+ds88g.du_frstressdz[:,:,:,:]) # Param Drag from stored data
dParam=dParam.mean(axis=0) # temopral average
dParam=dParam.mean(axis=2) # zonal average



########### Figure 1 ###########################

var1 = np.zeros((181, 361))  # Replace size_of_var1 with the appropriate size
var2 = np.zeros((181, 361))  # Replace size_of_var2 with the appropriate size
var3 = np.zeros((181, 361))  # Replace size_of_var3 with the appropriate size
var4 = np.zeros((181, 361))  # Replace size_of_var4 with the appropriate size
var5 = np.zeros((181, 361))  # Replace size_of_var5 with the appropriate size
var6 = np.zeros((181, 361))  # Replace size_of_var6 with the appropriate size
longit = np.zeros(361)  # Assuming longit is a 1D array

lev=29  # Vertical level 1
var1[:,0:360]=ICON[lev,:,:]
var2[:,0:360]=GW[lev,:,:]

lev=49 # Vertical level 2
var3[:,0:360]=ICON[lev,:,:]
var4[:,0:360]=GW[lev,:,:]

lev=59 # Vertical level 3
var5[:,0:360]=ICON[lev,:,:]
var6[:,0:360]=GW[lev,:,:]

var1[:,360]=var1[:,359]
var2[:,360]=var2[:,359]
var3[:,360]=var3[:,359]
var4[:,360]=var4[:,359]
var5[:,360]=var5[:,359]
var6[:,360]=var6[:,359]

longit[0:360] = ds2.longitude
longit[360]=360.
worldmap_fig1(var1,var2,var3,var4,var5,var6,longit,lat)
plt.savefig('fig1.png', format='png')



########### Figure 2 ###########################


var1 = np.zeros((181, 361))  
var2 = np.zeros((181, 361))  
var3 = np.zeros((181, 361))  
var4 = np.zeros((181, 361)) 
var5 = np.zeros((181, 361)) 
var6 = np.zeros((181, 361))  

lev=49 # Vertical level 1
var1[:,0:360]=GW_or[lev,:,:]
var2[:,0:360]=GW_pr[lev,:,:]
var3[:,0:360]=GW_fr[lev,:,:]

lev=59 # Vertical level 2
var4[:,0:360]=GW_or[lev,:,:]
var5[:,0:360]=GW_pr[lev,:,:]
var6[:,0:360]=GW_fr[lev,:,:]


var1[:,360]=var1[:,359]
var2[:,360]=var2[:,359]
var3[:,360]=var3[:,359]
var4[:,360]=var4[:,359]
var5[:,360]=var5[:,359]
var6[:,360]=var6[:,359]
longit = np.zeros(361)  # Assuming longit is a 1D array
longit[0:360] = ds2.longitude
longit[360]=360.

worldmap_fig2(var1,var2,var3,var4,var5,var6,longit,lat)
plt.savefig('fig2.png', format='png')


########### Figure 3 ###########################


var1 = np.zeros((181, 361))  
var2 = np.zeros((181, 361))  
var3 = np.zeros((181, 361))  
var4 = np.zeros((181, 361)) 
var5 = np.zeros((181, 361)) 
var6 = np.zeros((181, 361))  
var7 = np.zeros((181, 361))  
var8 = np.zeros((181, 361))  

lev=49 # Vertical level 1
var1[:,0:360]=ICONp[lev,:,:]
var2[:,0:360]=GWp[lev,:,:]
var3[:,0:360]=ICONn[lev,:,:]
var4[:,0:360]=GWn[lev,:,:]

lev=59 # Vertical level 2
var5[:,0:360]=ICONp[lev,:,:]
var6[:,0:360]=GWp[lev,:,:]
var7[:,0:360]=ICONn[lev,:,:]
var8[:,0:360]=GWn[lev,:,:]


var1[:,360]=var1[:,359]
var2[:,360]=var2[:,359]
var3[:,360]=var3[:,359]
var4[:,360]=var4[:,359]
var5[:,360]=var5[:,359]
var6[:,360]=var6[:,359]
var7[:,360]=var7[:,359]
var8[:,360]=var8[:,359]

worldmap_fig3(var1,var2,var3,var4,var5,var6,var7,var8,longit,lat)
plt.savefig('fig3.png', format='png')




########### Figure 4 ###########################




# Minimum and maximum values for contour plots
vmin = -4
vmax = 4

# Define custom ticks and levels
custom_ticks = [-4, -3, -2, -1, -0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8, 1, 2, 3, 4]
custom_levels = [-4, -3.5, -3, -2.5, -2, -1.5, -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 2.5, 3, 3.5, 4]
levels_cl_p = np.linspace(0, 140, 20)  # Adjust the number of levels as needed
levels_cl_n = np.linspace(-40, 0, 20)  # Adjust the number of levels as needed

# Get the original coolwarm colormap
cmap = get_cmap('coolwarm')

# Create a normalization object
norm = BoundaryNorm(custom_levels, ncolors=cmap.N, clip=True)

# Create a figure with two columns
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10), gridspec_kw={'width_ratios': [1, 1], 'wspace': 0.05})
plt.rcParams.update({'font.size': 24})
LW = 3

# Plot the 'dICON' data as a contour plot on the first axis
contour1 = ax1.contourf(ds88g['latitude'], ds88g['level'], -dICON-dEPFphi, cmap=cmap, levels=custom_levels, norm=norm, extend='both')

# Plot the 'dParam' data as a contour plot on the second axis
contour2 = ax2.contourf(ds88g['latitude'], ds88g['level'], -dParam, cmap=cmap, levels=custom_levels, norm=norm, extend='both')

# Adjust the position of the second subplot to the right
box = ax2.get_position()
ax2.set_position([box.x0 + 0.02, box.y0, box.width, box.height])

# Add a single colorbar to the right of the figure
cbar = fig.colorbar(contour1, ax=[ax1, ax2], orientation='vertical', shrink=0.8, pad=0.04)
cbar.set_label('GWD (m/s/day)')
cbar.set_ticks(custom_ticks)
cbar.ax.tick_params(labelsize=20)  # Adjust the size of the ticks

# Plot the contours for the 'vitu' data on both contour axes
contour_lines_p1 = ax1.contour(ds88g['latitude'], ds88g['level'], vitu1, levels=levels_cl_p, colors='black', linestyles='solid')
ax1.clabel(contour_lines_p1, inline=True, fontsize=12, fmt='%d')
contour_lines_n1 = ax1.contour(ds88g['latitude'], ds88g['level'], vitu1, levels=levels_cl_n, colors='black', linestyles='dashed')
ax1.clabel(contour_lines_n1, inline=True, fontsize=12, fmt='%d')

contour_lines_p2 = ax2.contour(ds88g['latitude'], ds88g['level'], vitu1, levels=levels_cl_p, colors='black', linestyles='solid')
ax2.clabel(contour_lines_p2, inline=True, fontsize=12, fmt='%d')
contour_lines_n2 = ax2.contour(ds88g['latitude'], ds88g['level'], vitu1, levels=levels_cl_n, colors='black', linestyles='dashed')
ax2.clabel(contour_lines_n2, inline=True, fontsize=12, fmt='%d')

# Set y-axis labels and ticks only for the first subplot
ax1.set_ylabel('Altitude')
ax1.set_ylim(0, 50)
ax1.set_xticks(np.arange(-90, 91, 30))  # Adjust the range and step as needed

# Set y-axis ticks for the second subplot but remove labels
ax2.set_ylim(0, 50)
ax2.set_xticks(np.arange(-90, 91, 30))  # Adjust the range and step as needed
ax2.yaxis.set_tick_params(which='both', labelleft=False)  # Keep ticks but remove labels

ax1.set_xlabel('Latitude')
ax2.set_xlabel('Latitude')

# Set titles
ax1.set_title('(a) ICON')
ax2.set_title('(b) Parameterization')
ax1.minorticks_on()
ax2.minorticks_on()

# Adjust layout to prevent overlap
plt.tight_layout(rect=[0, 0, 0.93, 1])  # Adjust the right margin to fit the colorbar without overlap

# Save the figure
plt.savefig('fig4.png', format='png')


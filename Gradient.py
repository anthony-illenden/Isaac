from datetime import datetime, timedelta 
 import cartopy.crs as ccrs 
 import matplotlib.patches as mpatches 
 import matplotlib.pyplot as plt 
 from metpy.io import parse_metar_file 
 import metpy.plots as mpplots 
 import numpy as np 
 from siphon.catalog import TDSCatalog 
 import xarray as xr 
 import cartopy.crs as ccrs 
 import cartopy. feature as cfeature 
 import matplotlib.pyplot as plt 
 from metpy.plots import USCOUNTIES 
 from metpy.interpolate import interpolate_to_grid, remove_nan_observations 
 import numpy as np 
 import metpy.calc as mpcalc 
 from metpy.units import units 
 from metpy.plots import StationPlot 
 import pandas as pd 
 from metpy.calc import equivalent_potential_temperature 
 from metpy.calc import (mixing_ratio_from_relative_humidity, relative_humidity_from_dewpoint, 
                         thickness_hydrostatic, heat_index) 
 import math 
  
 def temp(tairC): 
     tairF = (tairC * 1.8) + 32 
     return tairF 
  
 def calculate_grid(df, var, interp_type='cressman'): 
     xp, yp, _ = mapcrs.transform_points(datacrs, df['longitude'].values, df['latitude'].values).T 
     x_masked, y_masked, z = remove_nan_observations(xp, yp, df[var].values) 
     gridx, gridy, z = interpolate_to_grid(x_masked, y_masked, z, interp_type=interp_type) 
     return gridx, gridy, z 
  
 def rh(temperature, dew_point): 
     es = 6.112 * math.exp((17.67 * dew_point) / (dew_point + 243.5)) 
     e = es * (humidity / 100) 
     RH = (e / es) * 100 
     return RH 
  
  
 def heatindex(tair, td): 
     heat_index = -42.379 + 2.04901523 * tair + 10.14333127 * td \ 
                  - 0.22475541 * tair * td - 0.00683783 * tair**2 \ 
                  - 0.05481717 * td**2 + 0.00122874 * tair**2 * td \ 
                  + 0.00085282 * tair * td**2 - 0.00000199 * tair**2 * td**2 
     heat_index = round(heat_index, 1) 
     return heat_index 
  
  
 md_time = datetime.utcnow() - timedelta(minutes=25) 
 #datacrs = ccrs.PlateCarree() 
 #mapcrs = ccrs.LambertConformal(central_longitude=-102, central_latitude=37) 
  
 # Get some METAR data 
 metar_cat = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/noaaport/text/metar/catalog.xml') 
 sfc_obs = parse_metar_file(metar_cat.datasets.filter_time_nearest(md_time).remote_open(mode='t')) 
 sfc_obs = sfc_obs.set_index('date_time').groupby('station_id').first() 
  
 df = sfc_obs 
  
 CONUS_BOUNDING_BOX = { 
     'min_lat': 24.396308, 
     'max_lat': 49.384358, 
     'min_lon': -125.000000, 
     'max_lon': -66.934570} 
  
 GL_BOUNDING_BOX = { 
     'min_lat': 41.5, 
     'max_lat': 47.5, 
     'min_lon': -90.5, 
     'max_lon': -82} 
  
 df = df[ 
     (sfc_obs['latitude'] >= CONUS_BOUNDING_BOX['min_lat']) & 
     (sfc_obs['latitude'] <= CONUS_BOUNDING_BOX['max_lat']) & 
     (sfc_obs['longitude'] >= CONUS_BOUNDING_BOX['min_lon']) & 
     (sfc_obs['longitude'] <= CONUS_BOUNDING_BOX['max_lon'])]         
  
  
 print (df['air_temperature'].dtypes) 
  
 mapcrs = ccrs.LambertConformal(central_longitude=-85.6, central_latitude=44.3, standard_parallels=(30, 60))  
 datacrs = ccrs.PlateCarree() 
 proj = ccrs.Stereographic(central_longitude=-85, central_latitude=40) 
  
 lat = df['latitude'].values 
 lon = df['longitude'].values 
  
 #df['air_temperature_f'] = df['air_temperature'].astype(float) 
  
 df['air_temperature_f'] = temp(df['air_temperature'].astype(float)) 
  
 #df['dew_temperature_f'] = df['dew_point_temperature'].astype(float) 
  
 #print(df.columns.values) 
  
  
  
 df['dew_temperature_f'] = temp(df['dew_point_temperature'].astype(float)) 
 df['mslp'] = df['air_pressure_at_sea_level'].astype(float) 
  
 tair = df['air_temperature_f'].values * units.degF 
 td = df['dew_temperature_f'].values * units.degF 
  
 df['rh'] = mpcalc.relative_humidity_from_dewpoint(tair, td).to('percent') 
  
 #rh = df['rh'] 
 #df['heat_index'] = mpcalc.heat_index(tair, rh) 
 #df['heat_index'] = heatindex(df['air_temperature_f'], df['dew_temperature_f']) 
 #df['RH'] = df.apply(lambda row: rh(df['air_temperature'], df['dew_point_temperature']), axis=1) 
  
 #df['heat_index'] = heatindex(df['air_temperature_f'].astype(float), df['air_temperature_f'].astype(float)) 
  
 #df_gl = df[(sfc_obs['latitude'] >= GL_BOUNDING_BOX['min_lat']) & 
 #    (sfc_obs['latitude'] <= GL_BOUNDING_BOX['max_lat']) & 
 #    (sfc_obs['longitude'] >= GL_BOUNDING_BOX['min_lon']) & 
 #    (sfc_obs['longitude'] <= GL_BOUNDING_BOX['max_lon'])] 
  
 #GL_df = df[ 
 #    (sfc_obs['latitude'] >= GL_BOUNDING_BOX['min_lat']) & 
 #    (sfc_obs['latitude'] <= GL_BOUNDING_BOX['max_lat']) & 
 #    (sfc_obs['longitude'] >= GL_BOUNDING_BOX['min_lon']) & 
 #    (sfc_obs['longitude'] <= GL_BOUNDING_BOX['max_lon'])] 
  
 gridx, gridy, tempf = calculate_grid(df, 'air_temperature_f')  
  
 tmin = df['dew_temperature_f'].min() 
 tmax = df['dew_temperature_f'].max() 
  
  
 fig = plt.figure(figsize=(16,12)) 
 ax = fig.add_subplot(1, 1, 1, projection=mapcrs) 
  
 ax.set_extent([-90.5, -82, 41.5, 47.5], datacrs) 
  
 ax.add_feature(cfeature.COASTLINE.with_scale('50m')) 
 ax.add_feature(cfeature.STATES.with_scale('50m')) 
 ax.add_feature(USCOUNTIES.with_scale('5m'), linewidth=0.25) 
  
 gridx_d, gridy_d, temp_d = calculate_grid(df, 'dew_temperature_f') 
  
 gridx, gridy, temp = calculate_grid(df, 'air_temperature_f') 
  
 gridx_p, gridy_p, pre_z = calculate_grid(df, 'mslp') 
  
 # Plot the data 
 im = ax.contourf(gridx, gridy, temp, levels=np.arange(30, 90, 5), cmap='Reds', alpha=0.45) 
 im_t = ax.contour(gridx, gridy, temp, levels=np.arange(30, 90, 5), cmap='Reds', alpha=0.45) 
 ax.clabel(im_t, inline=True, colors='k') 
  
 #ax.contourf(gridx_d, gridy_d, temp_d, levels=np.arange(30, 75, 5), cmap='Greens', alpha=0.45) 
 #ax.contour(gridx_d, gridy_d, temp_d, levels=np.arange(30, 75, 5), colors='k', alpha=0.45) 
 p_im = ax.contour(gridx_p, gridy_p, pre_z, levels=np.arange(1000, 1020, 2), colors='k', alpha=0.75) 
 ax.clabel(p_im, inline=True, colors='k') 
  
 cbar = plt.colorbar(im, ax=ax); cbar.ax.set_ylabel('Air Temperature (F)') 
  
 #sc = ax.scatter(df['longitude'], df['latitude'], c = df['dew_temperature_f'], cmap = 'Greens',  
 #               vmin=tmin, vmax=tmax, transform=datacrs, alpha = 0) 
 #cb = plt.colorbar(sc) 
 #cb.set_label('Temperature') 
  
 u, v = mpcalc.wind_components(df['wind_speed'].values * units.mph, df['wind_direction'].values * units.degrees) 
  
 #u, v = mpcalc.wind_components(GL_df['wind_speed'].values * units.mph, GL_df['wind_direction'].values * units.degrees) 
  
  
 stationplot = StationPlot(ax, lon, lat, transform=datacrs, 
                           fontsize=10) 
 stationplot.plot_parameter('NW', df['air_temperature_f'], color='red') 
 stationplot.plot_parameter('SW', df['dew_temperature_f'], color='green') 
 stationplot.plot_barb(u, v) 
  
  
 ax.set_title(f'Current Surface Observations: {md_time:%Y-%m-%d %H:%M}Z') 
  
 plt.savefig("Surface_Analysis_GH_.jpg")

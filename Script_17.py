# Training: Python and GOES-R Imagery: Script 17 - Level 2 Products (RRQPE) and Data Accumulation
#-----------------------------------------------------------------------------------------------------------
# Required modules
from netCDF4 import Dataset                     # Read / Write NetCDF4 files
import matplotlib.pyplot as plt                 # Plotting library
from datetime import timedelta, date, datetime  # Basic Dates and time types
import cartopy, cartopy.crs as ccrs             # Plot maps
import os                                       # Miscellaneous operating system interfaces
from osgeo import gdal                          # Python bindings for GDAL
import numpy as np                              # Scientific computing with Python
from matplotlib import cm                       # Colormap handling utilities
from utilities import download_PROD             # Our function for download
from utilities import reproject                 # Our function for reproject
gdal.PushErrorHandler('CPLQuietErrorHandler')   # Ignore GDAL warnings

#-----------------------------------------------------------------------------------------------------------
# Input and output directories
input = "Samples"; os.makedirs(input, exist_ok=True)
output = "Output"; os.makedirs(output, exist_ok=True)

# Desired extent
extent = [-74.0, -34.1, -34.8, 5.5] # Min lon, Max lon, Min lat, Max lat

# AMAZON repository information 
# https://noaa-goes16.s3.amazonaws.com/index.html
bucket_name = 'noaa-goes16'
product_name = 'ABI-L2-RRQPEF'
yyyymmdd = '20201217'

########################################################################
# Rainfall Rate Quantitative Precipitation Estimation ("X" hours)
########################################################################

# Initial time and date
yyyy = datetime.strptime(yyyymmdd, '%Y%m%d').strftime('%Y')
mm = datetime.strptime(yyyymmdd, '%Y%m%d').strftime('%m')
dd = datetime.strptime(yyyymmdd, '%Y%m%d').strftime('%d')

date_ini = str(datetime(int(yyyy),int(mm),int(dd),12,0) - timedelta(hours=23))
date_end = str(datetime(int(yyyy),int(mm),int(dd),12,0))

acum = np.zeros((5424,5424))

#-----------------------------------------------------------------------------------------------------------
# Accumulation loop
while (date_ini <= date_end):

    # Date structure
    yyyymmddhhmn = datetime.strptime(date_ini, '%Y-%m-%d %H:%M:%S').strftime('%Y%m%d%H%M')

    # Download the file
    file_name = download_PROD(yyyymmddhhmn, product_name, input)
    #-----------------------------------------------------------------------------------------------------------
    # Variable
    var = 'RRQPE'

    # Open the file
    img = gdal.Open(f'NETCDF:{input}/{file_name}.nc:' + var)
    dqf = gdal.Open(f'NETCDF:{input}/{file_name}.nc:DQF')

    # Read the header metadata
    metadata = img.GetMetadata()
    scale = float(metadata.get(var + '#scale_factor'))
    offset = float(metadata.get(var + '#add_offset'))
    undef = float(metadata.get(var + '#_FillValue'))
    dtime = metadata.get('NC_GLOBAL#time_coverage_start')

    # Load the data
    ds = img.ReadAsArray(0, 0, img.RasterXSize, img.RasterYSize).astype(float)
    ds_dqf = dqf.ReadAsArray(0, 0, dqf.RasterXSize, dqf.RasterYSize).astype(float)

    # Remove undef
    ds[ds == undef] = np.nan

    # Apply the scale, offset and convert to celsius
    ds = (ds * scale + offset)

    # Apply NaN's where the quality flag is greater than 1
    ds[ds_dqf > 0] = np.nan

    # Sum the instantaneous value in the accumulation
    acum = np.nansum(np.dstack((acum, ds)),2)

    # Increment 1 hour
    date_ini = str(datetime.strptime(date_ini, '%Y-%m-%d %H:%M:%S') + timedelta(hours=1))
    #-----------------------------------------------------------------------------------------------------------

# Reproject the file
filename_acum = f'{output}/prec_acum_ret_{yyyymmddhhmn}.nc'
reproject(filename_acum, img, acum, extent, undef)
#-----------------------------------------------------------------------------------------------------------
# Open the reprojected GOES-R image
file = Dataset(filename_acum)

# Get the pixel values
data = file.variables['Band1'][:]
#-----------------------------------------------------------------------------------------------------------
# Choose the plot size (width x height, in inches)
plt.figure(figsize=(10,10))

# Use the Geostationary projection in cartopy
ax = plt.axes(projection=ccrs.PlateCarree())

# Define the image extent
img_extent = [extent[0], extent[2], extent[1], extent[3]]
 
# Modify the colormap to zero values are white
colormap = cm.get_cmap('rainbow', 240)
newcolormap = colormap(np.linspace(0, 1, 240))
newcolormap[:1, :] = np.array([1, 1, 1, 1])
cmap = cm.colors.ListedColormap(newcolormap)

# Plot the image
img = ax.imshow(data, vmin=0, vmax=150, cmap=cmap, origin='upper', extent=img_extent)

# Add coastlines, borders and gridlines
ax.coastlines(resolution='10m', color='black', linewidth=0.8)
ax.add_feature(cartopy.feature.BORDERS, edgecolor='black', linewidth=0.5)
gl = ax.gridlines(crs=ccrs.PlateCarree(), color='gray', alpha=1.0, linestyle='--', linewidth=0.25, xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5), draw_labels=True)
gl.top_labels = False
gl.right_labels = False
    
# Add a colorbar
plt.colorbar(img, label='Rainfall Rate mm / 24h', extend='max', orientation='horizontal', pad=0.05, fraction=0.05)

# Extract the date
date = (datetime.strptime(dtime, '%Y-%m-%dT%H:%M:%S.%fZ'))

# Add a title
plt.title('G-16 ACCUM. PREC.', fontweight='bold', fontsize=10, loc='left')
date_ini = datetime(int(yyyy),int(mm),int(dd),12,0) - timedelta(hours=23)
date_end = datetime(int(yyyy),int(mm),int(dd),12,0)
plt.title(date_ini.strftime('%Y-%m-%d %H:%M') + " - " + date_end.strftime('%Y-%m-%d %H:%M'), fontsize=10, loc='right')
#-----------------------------------------------------------------------------------------------------------
# Save the image
plt.savefig(f'{output}/PREC_ACUM_RET_{yyyymmddhhmn}.png', bbox_inches='tight', pad_inches=0, dpi=300)

# Show the image
plt.show()

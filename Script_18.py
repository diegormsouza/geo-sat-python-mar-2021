# Training: Python and GOES-R Imagery: Script 18 - Comparative of (RRQPE) with accumulated data of rain gauges
#-----------------------------------------------------------------------------------------------------------
import os
import matplotlib.pyplot as plt                 # Plotting library
import pandas as pd                             # Read and manipulate CSV file
from osgeo import gdal                          # Python bindings for GDAL

#-----------------------------------------------------------------------------------------------------------
# Input and output directories
input = "Samples"; os.makedirs(input, exist_ok=True)
output = "Output"; os.makedirs(output, exist_ok=True)

yyyymmdd = "20201217"

# Read file csv with observated rain
obs_data = pd.read_csv(f'pluvio_inmet_{yyyymmdd}.txt')

#  Read file netcdf with estimated rain from GOES-16
sat_data = gdal.Open(f'{output}/prec_acum_ret_{yyyymmdd}1200.nc')

# Read number of cols and rows
ncol = sat_data.RasterXSize
nrow = sat_data.RasterYSize

# Load the data
sat_array = sat_data.ReadAsArray(0, 0, ncol, nrow).astype(float)

# Get geotransform
transform = sat_data.GetGeoTransform()

# Create a table with both of values
tab_prec = pd.DataFrame(columns=["lon","lat","obs","sat"])

idx = 0
for obs in obs_data["ACUM"]:
  lon = obs_data.loc[idx,"LON"]
  lat = obs_data.loc[idx,"LAT"]

  x = int((lon - transform[0]) / transform[1])
  y = int((transform[3] - lat) / -transform[5])
  
  if x <= ncol and y <= nrow:
    sat = sat_array[y,x]
    tab_prec.loc[idx,"lon"] = lon
    tab_prec.loc[idx,"lat"] = lat
    tab_prec.loc[idx,"obs"] = obs
    tab_prec.loc[idx,"sat"] = sat

  idx += 1

fig, ax = plt.subplots(figsize=(10, 6))
ax.scatter(x = tab_prec['obs'], y = tab_prec['sat'])
plt.xlim(0, 100)
plt.ylim(0, 100)
ax.plot([0, 100], [0, 100], ls="--", c="r")
plt.xlabel("Observated")
plt.ylabel("Estimated")
plt.savefig(f'{output}/SATxOBS_{yyyymmdd}.png', bbox_inches='tight', pad_inches=0, dpi=300)

plt.show()

print(tab_prec)

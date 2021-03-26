# Training: Python and GOES-R Imagery: Script 1 - Basic Plot / Extracting Pixel Values
#----------------------------------------------------------------------------------------------------------- 
# Required modules
from netCDF4 import Dataset      # Read / Write NetCDF4 files
import matplotlib.pyplot as plt  # Plotting library
#----------------------------------------------------------------------------------------------------------- 
# Open the GOES-R image 
# Download files at this link: http://home.chpc.utah.edu/~u0553130/Brian_Blaylock/cgi-bin/goes16_download.cgi
file = Dataset("OR_ABI-L2-CMIPF-M6C13_G16_s20191981200396_e20191981210116_c20191981210189.nc")
 
# Get the pixel values
data = file.variables['CMI'][:]
#-----------------------------------------------------------------------------------------------------------
# Choose the plot size (width x height, in inches)
plt.figure(figsize=(7,7))
 
# Plot the image
plt.imshow(data, vmin=193, vmax=313, cmap='Greys')
#----------------------------------------------------------------------------------------------------------- 
# Save the image
plt.savefig('Image_01.png')

# Show the image
plt.show()



"""
HMI PFSS solutions
------------------
Calculating a PFSS solution from a HMI synoptic map.

This example shows how to calcualte a PFSS solution from a HMI synoptic map.
There are a couple of important things that this example shows:

- HMI maps have non-standard metadata, so this needs to be fixed
- HMI synoptic maps are very big (1440 x 3600), so need to be downsampled
  in order to calculate the PFSS solution in a reasonable time.
"""
import os

import astropy.units as u
import matplotlib.pyplot as plt
import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a
import numpy as np

import pfsspy
import pfsspy.utils

###############################################################################
# Set up the search.
#
# Note that for SunPy versions earlier than 2.0, a time attribute is needed to
# do the search, even if (in this case) it isn't used, as the synoptic maps are
# labelled by Carrington rotation number instead of time
time = a.Time('2010/01/01', '2010/01/01')
series = a.jsoc.Series('hmi.synoptic_mr_720s')
crot = a.jsoc.PrimeKey('CAR_ROT', 2294)

###############################################################################
# Do the search.
#
# If you use this code, please replace this email address
# with your own one, registered here:
# http://jsoc.stanford.edu/ajax/register_email.html
email = "chjan0809@gmail.com"
result = Fido.search(time, series, crot,
                     a.jsoc.Notify(email))
files = Fido.fetch(result)

###############################################################################
# Read in a file. This will read in the first file downloaded to a sunpy Map
# object
hmi_map = sunpy.map.Map(files[0])

print('Data shape: ', hmi_map.data.shape)

###############################################################################
# Since this map is far to big to calculate a PFSS solution quickly, lets
# resample it down to a smaller size.
hmi_map = hmi_map.resample([360, 180] * u.pix)
hmi_map.data[np.isnan(hmi_map.data)] = np.nanmedian(hmi_map.data)
print('New shape: ', hmi_map.data.shape)
# print(f"NaN 數量: {np.isnan(hmi_map.data).sum()}")
###############################################################################
# Now calculate the PFSS solution
nrho = 70
rss = 3.5
pfss_in = pfsspy.Input(hmi_map, nrho, rss)
pfss_out = pfsspy.pfss(pfss_in)

###############################################################################
# Using the Output object we can plot the source surface field, and the
# polarity inversion line.
ss_br = pfss_out.source_surface_br
# Create the figure and axes
fig = plt.figure()
ax = plt.subplot(projection=ss_br)
# print(ss_br)
# Plot the source surface map
ss_br.plot()
# Plot the polarity inversion line
ax.plot_coord(pfss_out.source_surface_pils[0])
# Plot formatting
plt.colorbar()
ax.set_title('Source surface magnetic field')

plt.show()

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
import matplotlib.colors as mcolors

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
print(hmi_map.meta)
print(hmi_map.unit)  # 單位應該是 Gauss (G)

# print(f"NaN 數量: {np.isnan(hmi_map.data).sum()}")
###############################################################################
# Now calculate the PFSS solution
nrho = 70
rss = 2.5
pfss_in = pfsspy.Input(hmi_map, nrho, rss)
pfss_out = pfsspy.pfss(pfss_in)

###############################################################################
# Using the Output object we can plot the source surface field, and the
# polarity inversion line.
# ss_br = pfss_out.source_surface_br
# # Create the figure and axes
# fig = plt.figure()
# ax = plt.subplot(projection=ss_br)
# # print(ss_br)
# # Plot the source surface map
# ss_br.plot()
# # Plot the polarity inversion line
# ax.plot_coord(pfss_out.source_surface_pils[0])
# # Plot formatting
# plt.colorbar()
# ax.set_title('Source surface magnetic field')

# 設定顏色條的範圍 (-100 到 100) 並均勻分佈
norm = mcolors.Normalize(vmin=-100, vmax=100)
hmi_map.plot_settings["norm"] = norm  # 設定顏色條範圍
# 繪製原始 HMI 磁場數據
fig = plt.figure()
ax = plt.subplot(projection=hmi_map)

# 繪製 HMI 磁場
hmi_map.plot(cmap="seismic_r")

# 顯示顏色條
plt.colorbar(label='Mr (Mx/cm^2)')

# 設定標題
ax.set_title('Original HMI Magnetic Field')

plt.show()

ss_br = pfss_out.source_surface_br

# 確保數據的形狀相符
print("HMI Shape:", hmi_map.data.shape)
print("SS_Br Shape:", ss_br.data.shape)

# 確保數據不包含 0，以免除法錯誤
hmi_data = hmi_map.data.copy()
hmi_data[hmi_data == 0] = np.nan  # 避免除以零

# 計算fs
fs = abs((1 / (2.5**2)) * (hmi_data / ss_br.data) )

# 創建 SunPy Map 物件
ratio_map = sunpy.map.Map(fs, hmi_map.meta)
# 設定顏色條的範圍
norm = mcolors.TwoSlopeNorm(vmin=0, vcenter=10, vmax=20)
ratio_map.plot_settings["norm"] = norm
# 繪製比值圖
# fig = plt.figure()
# ax = plt.subplot(projection=ratio_map)
# ratio_map.plot()  
# plt.colorbar(label="fs")
# ax.set_title("PFSS Expansion factor")

# plt.show()

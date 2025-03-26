import sunpy.map
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import astropy.units as u
from astropy.wcs import WCS
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

# 讀取 FITS 檔案
file_path_hmi = "20231104/hmi.m_720s.20231104_000000_TAI.3.magnetogram.fits"
file_path_aia = "20231104/aia.lev1_euv_12s.2023-11-04T000006Z.193.image_lev1.fits"
hmi_map = sunpy.map.Map(file_path_hmi)
aia_map = sunpy.map.Map(file_path_aia)

# 讀取header
hmi_header = hmi_map.wcs.to_header()
aia_header = aia_map.wcs.to_header()
for key, value in hmi_header.items():
    print(f"{key}: {value}")
# 繪製hmi (灰階)
fig = plt.figure(figsize = (20, 10))
ax1 = plt.subplot(122, projection=hmi_map)

# 設定hmi顯示
norm = mcolors.Normalize(vmin=-200, vmax=200)
hmi_map.plot_settings["norm"] = norm  # 設定顏色條範圍
hmi_map.plot(cmap = "gray")
cbar = plt.colorbar(ax=ax1, label="Magnetic Field Strength (G)", fraction=0.046, pad=0.04)
cbar.ax.tick_params(size=10)
x_range11 = (1050/3600-hmi_header['CRVAL1'])/hmi_header['CDELT1']+hmi_header['CRPIX1']
x_range12 = (-1050/3600-hmi_header['CRVAL1'])/hmi_header['CDELT1']+hmi_header['CRPIX1']
y_range11 = (1050/3600-hmi_header['CRVAL2'])/hmi_header['CDELT2']+hmi_header['CRPIX2']
y_range12 = (-1050/3600-hmi_header['CRVAL2'])/hmi_header['CDELT2']+hmi_header['CRPIX2']
ax1.set_xlim([x_range11, x_range12])
ax1.set_ylim([y_range11, y_range12])
# 設定aia顯示
ax2 = plt.subplot(121, projection = aia_map)
aia_map.plot(cmap = "sdoaia193")
x_range21 = (-1000/3600-aia_header['CRVAL1'])/aia_header['CDELT1']+aia_header['CRPIX1']
x_range22 = (1000/3600-aia_header['CRVAL1'])/aia_header['CDELT1']+aia_header['CRPIX1']
y_range21 = (-1000/3600-aia_header['CRVAL2'])/aia_header['CDELT2']+aia_header['CRPIX2']
y_range22 = (1000/3600-aia_header['CRVAL2'])/aia_header['CDELT2']+aia_header['CRPIX2']
ax2.set_xlim([x_range21, x_range22])
ax2.set_ylim([y_range21, y_range22])
plt.show()


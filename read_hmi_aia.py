import sunpy.map
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import astropy.units as u
from astropy.wcs import WCS
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import pandas as pd

# 讀取 FITS 檔案
file_path_hmi = "./data_1271/hmi.m_45s.20250225_162145_TAI.2.magnetogram.fits"
file_path_aia = "./data_1271/aia.lev1_euv_12s.2025-02-25T162106Z.193.image_lev1.fits"
hmi_map = sunpy.map.Map(file_path_hmi)
aia_map = sunpy.map.Map(file_path_aia)
aia_df = pd.DataFrame(aia_map.data)

# 讀取header
hmi_header = hmi_map.wcs.to_header()
aia_header = aia_map.wcs.to_header()
for key, value in hmi_header.items():
    print(f"{key}: {value}")
# 繪製hmi (灰階)
fig = plt.figure(figsize = (20, 10))
ax1 = plt.subplot(121, projection=hmi_map)

threshold = 65
# 設定hmi顯示
norm = mcolors.Normalize(vmin=-200, vmax=200)
# img1 = ax1.contour(aia_df, levels=[threshold], colors='white', linewidths=1)
hmi_map.plot_settings["norm"] = norm  # 設定顏色條範圍
img1 = hmi_map.plot(cmap = "gray")
cbar1 = plt.colorbar(ax=ax1, label="Magnetic Field Strength (G)",
                     fraction=0.046, pad=0.04,
                     location='right')
cbar1.ax.tick_params(size=10)
x_range11 = (1050/3600-hmi_header['CRVAL1'])/hmi_header['CDELT1']+hmi_header['CRPIX1']
x_range12 = (-1050/3600-hmi_header['CRVAL1'])/hmi_header['CDELT1']+hmi_header['CRPIX1']
y_range11 = (1050/3600-hmi_header['CRVAL2'])/hmi_header['CDELT2']+hmi_header['CRPIX2']
y_range12 = (-1050/3600-hmi_header['CRVAL2'])/hmi_header['CDELT2']+hmi_header['CRPIX2']
ax1.set_xlim([x_range11, x_range12])
ax1.set_ylim([y_range11, y_range12])
# 設定aia顯示
ax2 = plt.subplot(122, projection = aia_map)
# img2 = ax2.contour(aia_df, levels=[threshold], colors='white', linewidths=1)
img2 = aia_map.plot(cmap = "sdoaia193")
cbar2 = plt.colorbar(ax=ax2, label="Strength", fraction=0.046, pad=0.04)
cbar2.ax.tick_params(size=10)
x_range21 = (-1000/3600-aia_header['CRVAL1'])/aia_header['CDELT1']+aia_header['CRPIX1']
x_range22 = (1000/3600-aia_header['CRVAL1'])/aia_header['CDELT1']+aia_header['CRPIX1']
y_range21 = (-1000/3600-aia_header['CRVAL2'])/aia_header['CDELT2']+aia_header['CRPIX2']
y_range22 = (1000/3600-aia_header['CRVAL2'])/aia_header['CDELT2']+aia_header['CRPIX2']
ax2.set_xlim([x_range21, x_range22])
ax2.set_ylim([y_range21, y_range22])
plt.subplots_adjust(wspace=0.4, top=1.1)
fig.suptitle("CH1271",fontsize = 20)
plt.show()


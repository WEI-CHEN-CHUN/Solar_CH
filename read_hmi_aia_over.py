import sunpy.map
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from datetime import datetime
import pandas as pd
from scipy.ndimage import zoom

# 讀取 FITS 檔案
file_path_hmi = "./data_1271/hmi.m_45s.20250225_162145_TAI.2.magnetogram.fits"
file_path_aia = "./data_1271/aia.lev1_euv_12s.2025-02-25T162106Z.193.image_lev1.fits"
hmi_map = sunpy.map.Map(file_path_hmi)
aia_map = sunpy.map.Map(file_path_aia)

# 讀取header
hmi_header = hmi_map.wcs.to_header()
aia_header = aia_map.wcs.to_header()
for key, value in hmi_header.items():
    print(f"{key}: {value}")

# 指定像素範圍
# x_range1 = (-1000/3600-aia_header['CRVAL1'])/aia_header['CDELT1']+aia_header['CRPIX1']
# x_range2 = (1000/3600-aia_header['CRVAL1'])/aia_header['CDELT1']+aia_header['CRPIX1']
# y_range1 = (-1000/3600-aia_header['CRVAL2'])/aia_header['CDELT2']+aia_header['CRPIX2']
# y_range2 = (1000/3600-aia_header['CRVAL2'])/aia_header['CDELT2']+aia_header['CRPIX2']

left_bottom = SkyCoord(-300 * u.arcsec, -150 * u.arcsec, frame=aia_map.coordinate_frame)
right_top = SkyCoord(250 * u.arcsec, 600 * u.arcsec, frame=aia_map.coordinate_frame)

# 裁剪地圖
aia_map = aia_map.submap(left_bottom, top_right=right_top)
hmi_map = hmi_map.submap(left_bottom, top_right=right_top)

# 轉換為 DataFrame
aia_df = pd.DataFrame(aia_map.data)
hmi_df = pd.DataFrame(hmi_map.data)
hmi_df = hmi_df.iloc[:, ::-1]#hmi
aia_df = aia_df.iloc[::-1, :]#aia
# 存成 CSV
aia_csv_filename = "data_1271/aia_CH1271_cropped_sp.csv"
hmi_csv_filename = "data_1271/hmi_CH1271_cropped_sp.csv"
aia_df.to_csv(aia_csv_filename, index=False, header=False)
hmi_df.to_csv(hmi_csv_filename, index=False, header=False)

print(f"已保存裁剪後的 HMI 磁場數據到 {hmi_csv_filename}")
print(f"已保存裁剪後的 AIA 磁場數據到 {aia_csv_filename}")
print(f"aia_df shape {aia_df.shape}")
print(f"aia_map shape {aia_map.data.shape}")
print(f"hmi_map shape {hmi_map.data.shape}")

zoom_factors = (hmi_map.data.shape[0] / aia_df.shape[0], hmi_map.data.shape[1] / aia_df.shape[1])
aia_df = zoom(aia_df, zoom_factors, order=1)
print(f"Resized aia_df shape: {aia_df.shape}")

# 繪製hmi (灰階)
fig = plt.figure(figsize = (10, 10))
ax = plt.subplot(projection=aia_map)
aia_map.plot(cmap = "sdoaia193", vmin=0, vmax=600, norm = mcolors.PowerNorm(1))
# 設定hmi顯示
norm = mcolors.Normalize(vmin=-200, vmax=200)
hmi_map.plot_settings["norm"] = norm  # 設定顏色條範圍
# img1 = ax.contour(aia_df, levels=[100], colors='white', linewidths=1)
img1 = hmi_map.plot(cmap = "gray", autoalign=True, alpha = 0)
cbar = plt.colorbar(ax=ax, label="Magnetic Field Strength (G)", fraction=0.046, pad=0.04)
cbar.ax.tick_params(size=10)

hmi_time = datetime.strptime(hmi_header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S.%f").strftime("%Y-%m-%d %H:%M:%S")
aia_time = datetime.strptime(aia_header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S.%f").strftime("%Y-%m-%d %H:%M:%S")
ax.set_title(F"HMI {hmi_time} overlaid on AIA193 {aia_time}")
plt.show()


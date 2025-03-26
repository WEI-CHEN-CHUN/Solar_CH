import sunpy.map
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from datetime import datetime
import pandas as pd
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

# 指定像素範圍
# x_range1 = (-1000/3600-aia_header['CRVAL1'])/aia_header['CDELT1']+aia_header['CRPIX1']
# x_range2 = (1000/3600-aia_header['CRVAL1'])/aia_header['CDELT1']+aia_header['CRPIX1']
# y_range1 = (-1000/3600-aia_header['CRVAL2'])/aia_header['CDELT2']+aia_header['CRPIX2']
# y_range2 = (1000/3600-aia_header['CRVAL2'])/aia_header['CDELT2']+aia_header['CRPIX2']
bottom_left = SkyCoord(-600 * u.arcsec, -650 * u.arcsec, frame=aia_map.coordinate_frame)
top_right = SkyCoord(100 * u.arcsec, 200 * u.arcsec, frame=aia_map.coordinate_frame)

# 裁剪地圖
aia_map = aia_map.submap(bottom_left, top_right=top_right)
hmi_map = hmi_map.submap(bottom_left, top_right=top_right)

# 轉換為 DataFrame
data_array = aia_map.data  # 獲取數據
df = pd.DataFrame(data_array)
# df = df.iloc[:, ::-1]#hmi
df = df.iloc[::-1, :]#aia
# 存成 CSV
csv_filename = "aia_20231104_cropped_sp.csv"
df.to_csv(csv_filename, index=False, header=False)

print(f"已保存裁剪後的 HMI 磁場數據到 {csv_filename}")


# 繪製hmi (灰階)
fig = plt.figure(figsize = (10, 10))
ax = plt.subplot(projection=aia_map)
aia_map.plot(cmap = "sdoaia193")
# 設定hmi顯示
norm = mcolors.Normalize(vmin=-200, vmax=200)
hmi_map.plot_settings["norm"] = norm  # 設定顏色條範圍
hmi_map.plot(cmap = "gray", autoalign=True, alpha = 0.5)
cbar = plt.colorbar(ax=ax, label="Magnetic Field Strength (G)", fraction=0.046, pad=0.04)
cbar.ax.tick_params(size=10)
# x_range11 = (-1000/3600-hmi_header['CRVAL1'])/hmi_header['CDELT1']+hmi_header['CRPIX1']
# x_range12 = (1000/3600-hmi_header['CRVAL1'])/hmi_header['CDELT1']+hmi_header['CRPIX1']
# y_range11 = (-1000/3600-hmi_header['CRVAL2'])/hmi_header['CDELT2']+hmi_header['CRPIX2']
# y_range12 = (1000/3600-hmi_header['CRVAL2'])/hmi_header['CDELT2']+hmi_header['CRPIX2']
# ax.set_xlim([x_range11, x_range12])
# ax.set_ylim([y_range11, y_range12])


# ax.set_xlim([x_range21, x_range22])
# ax.set_ylim([y_range21, y_range22])

hmi_time = datetime.strptime(hmi_header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S.%f").strftime("%Y-%m-%d %H:%M:%S")
aia_time = datetime.strptime(aia_header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S.%f").strftime("%Y-%m-%d %H:%M:%S")
ax.set_title(F"HMI {hmi_time} overlaid on AIA193 {aia_time}")
plt.show()


import pfsspy
import sunpy.map
from pfsspy import tracing
from astropy.constants import R_sun
import numpy as np
import astropy.units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.mplot3d import Axes3D

# 讀取 HMI 磁場數據
filename = 'hmi.synoptic_mr_720s.2294.synopMr.fits'
hmi_map = sunpy.map.Map(filename)
header = hmi_map.wcs.to_header()
print('Data shape: ', hmi_map.data.shape)
for key, value in header.items():
    print(f"{key}: {value}")
hmi_map = hmi_map.resample([360, 180] * u.pix)
hmi_map.data[np.isnan(hmi_map.data)] = np.nanmedian(hmi_map.data)

# 設定 PFSS 計算範圍
nr = 50  # 徑向網格數
rss = 2.5  # Source Surface 半徑 (單位: R_sun)
pfss_input = pfsspy.Input(hmi_map, nr, rss)
pfss_model = pfsspy.pfss(pfss_input)

# 設定磁力線起始點
np.random.seed(42)
n_streamlines = 50  # 追蹤的磁力線數量
latitudes = np.linspace(-np.pi/2, np.pi/2, n_streamlines)  # 均勻選擇緯度
longitudes = np.linspace(0, 2*np.pi, n_streamlines)  # 均勻選擇經度
radii = np.ones(n_streamlines) * R_sun.value  # 在太陽表面

# 建立 SkyCoord 起始點 (球坐標: r, lat, lon)
seeds = SkyCoord(longitudes * u.rad, 
                 latitudes * u.rad,  
                 radii * u.R_sun,
                 frame=pfss_model.coordinate_frame)

# 追蹤磁力線
tracer = tracing.PythonTracer()
field_lines = tracer.trace(seeds, pfss_model)

# 繪製 Carrington 投影圖
fig = plt.figure(figsize = (20, 10))
ax = plt.subplot(111, projection=hmi_map)
ax.set_xlim(0, 360)
ax.set_ylim(-90, 90)
ax.set_xlabel('Carrington Longitude [deg]')
ax.set_ylabel('Latitude [deg]')
ax.set_title('PFSS Magnetic Field Lines (Carrington Projection)')

# 繪製 HMI 磁場背景
norm = mcolors.Normalize(vmin=-200, vmax=200)  # 設定顏色範圍
hmi_map.plot_settings["norm"] = norm  
hmi_map.plot(cmap="gray")

# 繪製磁力線
for field_line in field_lines:
    ax.plot(field_line.lon.to_value(u.deg), 
            field_line.lat.to_value(u.deg), 
            color='b', alpha=0.7)

plt.show()

# # 繪製hmi (灰階)
# fig = plt.figure(figsize = (20, 10))
# ax1 = plt.subplot(111, projection=hmi_map)
# # 設定hmi顯示
# norm = mcolors.Normalize(vmin=-200, vmax=200)
# hmi_map.plot_settings["norm"] = norm  # 設定顏色條範圍
# hmi_map.plot(cmap = "gray")
# cbar = plt.colorbar(ax=ax1, label="Magnetic Field Strength (G)", fraction=0.046, pad=0.04)
# cbar.ax.tick_params(size=10)

# plt.show()
from scipy.io import readsav
import matplotlib.colors as mcolors
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import sunpy.visualization.colormaps as cm
from scipy.ndimage import zoom
import sunpy.map
from datetime import datetime


# 讀取 IDL 的 .sav 檔案
data_pfss = readsav("C:/Users/chjan/hmi_output_CH1271.sav")
# read aia
aia = pd.read_csv("data_1271/aia_CH1271_cropped_sp.csv")

# 讀取 FITS 檔案
file_path_hmi = "./data_1271/hmi.m_45s.20250225_162145_TAI.2.magnetogram.fits"
file_path_aia = "./data_1271/aia.lev1_euv_12s.2025-02-25T162106Z.193.image_lev1.fits"
hmi_map = sunpy.map.Map(file_path_hmi)
aia_map = sunpy.map.Map(file_path_aia)
# 讀取header
hmi_header = hmi_map.wcs.to_header()
aia_header = aia_map.wcs.to_header()
hmi_time = datetime.strptime(hmi_header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S.%f").strftime("%Y-%m-%d %H:%M:%S")
aia_time = datetime.strptime(aia_header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S.%f").strftime("%Y-%m-%d %H:%M:%S")

# 列出所有變數名稱
print(data_pfss.keys())

# 取出 BP3DZ 資料
BP3DX = data_pfss["BP3DX"].T
BP3DY = data_pfss["BP3DY"].T
BP3DZ = data_pfss["BP3DZ"].T
print("shape:", BP3DX.shape, BP3DY.shape, BP3DZ.shape)

# zoom
zoom_factors = (BP3DZ[:,:,0].shape[0] / aia.shape[0], BP3DZ[:,:,0].shape[1] / aia.shape[1])
aia = zoom(aia, zoom_factors, order=1)


# 顯示 BP3DZ
fig = plt.figure(figsize=(16,8))  # 設定圖大小
ax1 = plt.subplot(121)
sdoaia193 = cm.cmlist["sdoaia193"]
# img1 = ax1.imshow(aia, cmap=sdoaia193, vmin=0,vmax=600)
# 繪製等高線，數值 = 100
img1 = ax1.contour(aia, levels=[50], colors='#fbfdaf', linewidths=1)
img1 = ax1.imshow(BP3DZ[:,:,0], cmap="gray",vmin=-5, vmax=5, alpha=1)  # 設定 colormap 和數值範圍
ax1.set_title(F"HMI {hmi_time} overlaid on AIA193 {aia_time}")
# ax1.set_xlabel("X-axis")
# ax1.set_ylabel("Y-axis")
cbar1 = plt.colorbar(img1, ax=ax1, fraction=0.046, pad=0.04)
cbar1.ax.tick_params(labelsize=10)
cbar1.ax.set_title('Bz (G)')
k=200

ax2 = plt.subplot(122)
# img2 = ax2.imshow(aia, cmap=sdoaia193, vmin=0,vmax=600)
img2 = ax2.contour(aia, levels=[60], colors='#fbfdaf', linewidths=1)
# img2 = ax2.imshow(BP3DZ[:,:,0]/BP3DZ[:,:,k], cmap="gray", vmin=0, vmax=10)  # 設定 colormap 和數值範圍
# ax2.set_title(f"Bp Z-dir bottom / layer {k*5}")
img2 = ax2.imshow(BP3DZ[:,:,k], cmap="gray")  # 設定 colormap 和數值範圍
ax2.set_title(f"Bp Z-dir layer {k*5}")
# ax2.set_xlabel("X-axis")
# ax2.set_ylabel("Y-axis")
cbar2 = plt.colorbar(img2, ax=ax2, fraction=0.046, pad=0.04)
cbar2.ax.tick_params(labelsize=10)
# cbar2.ax.set_title('Bz (G)')

# plt.savefig("output_image/B300.jpg",bboxes_inches = 'tight')
ax1.minorticks_on()
ax2.minorticks_on()
plt.show()


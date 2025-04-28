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
data = readsav("C:/Users/chjan/hmi_output_CH1271.sav")
# read aia
aia = pd.read_csv("data_1271/aia_CH1271_cropped_sp.csv")

# 讀取 FITS 檔案
file_path_hmi = "./data_1271/hmi.m_45s.20250225_162145_TAI.2.magnetogram.fits"
file_path_aia = "./data_1271/aia.lev1p5_euv_12s.2025-02-25T162106Z.193.image_lev1p5.fits"
hmi_map = sunpy.map.Map(file_path_hmi)
aia_map = sunpy.map.Map(file_path_aia)
# 讀取header
hmi_header = hmi_map.wcs.to_header()
aia_header = aia_map.wcs.to_header()
hmi_time = datetime.strptime(hmi_header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S.%f").strftime("%Y-%m-%d %H:%M:%S")
aia_time = datetime.strptime(aia_header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S.%f").strftime("%Y-%m-%d %H:%M:%S")

# 列出所有變數名稱
print(data.keys())

# 取出 BP3DZ 資料
BP3DX = data["BP3DX"].T
BP3DY = data["BP3DY"].T
BP3DZ = data["BP3DZ"].T
print("shape:", BP3DX.shape, BP3DY.shape, BP3DZ.shape)

# zoom
zoom_factors = (BP3DZ[:,:,0].shape[0] / aia.shape[0], BP3DZ[:,:,0].shape[1] / aia.shape[1])
aia = zoom(aia, zoom_factors, order=1)

vmin=0; vmax=5
ft=8
# 顯示 BP3DZ
fig = plt.figure(figsize=(16,8))  # 設定圖大小
# ax1 = plt.subplot(121)
# sdoaia193 = cm.cmlist["sdoaia193"]
# # img1 = ax1.imshow(aia, cmap=sdoaia193, vmin=0,vmax=600)
# # 繪製等高線，數值 = 100
# img1 = ax1.contour(aia, levels=[60], colors='#fbfdaf', linewidths=1)
# img1 = ax1.imshow(BP3DZ[:,:,0], cmap="gray", alpha=1)  # 設定 colormap 和數值範圍
# ax1.set_title(F"HMI {hmi_time} overlaid on AIA193 {aia_time}")
# # ax1.set_xlabel("X-axis")
# # ax1.set_ylabel("Y-axis")
# cbar1 = plt.colorbar(img1, ax=ax1, fraction=0.046, pad=0.04)
# cbar1.ax.tick_params(labelsize=10)
# cbar1.ax.set_title('Bz (G)')
end_layer = 300
d = 5 # every d layer save
ax2 = plt.subplot(131)
# img2 = ax2.imshow(aia, cmap=sdoaia193, vmin=0,vmax=600)
img2 = ax2.contour(aia, levels=[60], colors='#fbfdaf', linewidths=1)
# img2 = ax2.imshow(BP3DZ[:,:,end_layer//d], cmap="gray", vmin=0, vmax=10)  # 設定 colormap 和數值範圍
# ax2.set_title(f"Bp layer {end_layer}")
img2 = ax2.imshow(BP3DZ[:,:,end_layer//d], cmap="gray", vmin=vmin, vmax=vmax)  # 設定 colormap 和數值範圍
ax2.set_title(f"Bz {end_layer}th layer , 1x region")
# ax2.set_xlabel("X-axis")
# ax2.set_ylabel("Y-axis")
cbar2 = plt.colorbar(img2, ax=ax2, fraction=0.046, pad=0.04)
cbar2.ax.tick_params(labelsize=10)
cbar2.ax.set_title('Bz (G)',fontsize=ft)

# plt.savefig("output_image/B300.jpg",bboxes_inches = 'tight')
# ax1.minorticks_on()
ax2.minorticks_on()
print(np.mean(BP3DZ[:,:,end_layer//d]))
"""1p5"""
data_1p5 = readsav("C:/Users/chjan/hmi_output_CH1271_1p5.sav")
aia_1p5 = pd.read_csv("data_1271/aia_CH1271_cropped_1p5_sp.csv")


# 取出 BP3DZ 資料
BP3DX = data_1p5["BP3DX"].T
BP3DY = data_1p5["BP3DY"].T
BP3DZ = data_1p5["BP3DZ"].T
print("shape:", BP3DX.shape, BP3DY.shape, BP3DZ.shape)

# zoom
zoom_factors = (BP3DZ[:,:,0].shape[0] / aia_1p5.shape[0], BP3DZ[:,:,0].shape[1] / aia_1p5.shape[1])
aia_1p5 = zoom(aia_1p5, zoom_factors, order=1)


# end_layer = 500
d = 10 # every d layer save
ax3 = plt.subplot(132)
img3 = ax3.contour(aia_1p5, levels=[60], colors='#fbfdaf', linewidths=1)

img3 = ax3.imshow(BP3DZ[:,:,end_layer//d], cmap="gray", vmin=vmin, vmax=vmax)  # 設定 colormap 和數值範圍
ax3.set_title(f"Bz {end_layer}th layer , 1.5x region")
# ax3.set_xlabel("X-axis")
# ax3.set_ylabel("Y-axis")
cbar3 = plt.colorbar(img3, ax=ax3, fraction=0.046, pad=0.04)
cbar3.ax.tick_params(labelsize=10)
ax3.minorticks_on()
cbar3.ax.set_title('Bz (G)',fontsize=ft)
# ax3.set_xlim(277, 277+1094)
# ax3.set_ylim(376+1489, 376)

print(np.mean(BP3DZ[376:376+1489, 277:277+1094,end_layer//d]))

"""3p"""
data_3p = readsav("C:/Users/chjan/hmi_output_CH1271_3p.sav")
aia_3p = pd.read_csv("data_1271/aia_CH1271_cropped_3p_sp.csv")


# 取出 BP3DZ 資料
BP3DX = data_3p["BP3DX"].T
BP3DY = data_3p["BP3DY"].T
BP3DZ = data_3p["BP3DZ"].T
print("shape:", BP3DX.shape, BP3DY.shape, BP3DZ.shape)

# zoom
zoom_factors = (BP3DZ[:,:,0].shape[0] / aia_3p.shape[0], BP3DZ[:,:,0].shape[1] / aia_3p.shape[1])
aia_3p = zoom(aia_3p, zoom_factors, order=1)


# end_layer = 500
d = 10 # every d layer save
ax4 = plt.subplot(133)
img4 = ax4.contour(aia_3p, levels=[60], colors='#fbfdaf', linewidths=1)

img4 = ax4.imshow(BP3DZ[:,:,end_layer//d], cmap="gray", vmin=vmin, vmax=vmax)  # 設定 colormap 和數值範圍
ax4.set_title(f"Bz {end_layer}th layer , 3x region")
# ax4.set_xlabel("X-axis")
# ax4.set_ylabel("Y-axis")
cbar4 = plt.colorbar(img4, ax=ax4, fraction=0.046, pad=0.04)
cbar4.ax.tick_params(labelsize=10)
ax4.minorticks_on()
cbar4.ax.set_title('Bz (G)',fontsize=ft)
# ax4.set_xlim(277*2, 277*2+1094)
# ax4.set_ylim(376*2+1489, 376*2)

print(np.mean(BP3DZ[376*2:376*2+1489, 277*2:277*2+1094,end_layer//d]))

plt.show()

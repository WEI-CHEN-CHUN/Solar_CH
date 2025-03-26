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
data_pfss = readsav("C:/Users/chjan/HMI_output_20231104.sav")
# read aia
aia = pd.read_csv('aia_20231104_cropped_sp.csv')

# 讀取 FITS 檔案
file_path_hmi = "20231104/hmi.m_720s.20231104_000000_TAI.3.magnetogram.fits"
file_path_aia = "20231104/aia.lev1_euv_12s.2023-11-04T000006Z.193.image_lev1.fits"
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

# 
zoom_factors = (BP3DZ[:,:,0].shape[0] / aia.shape[0], BP3DZ[:,:,0].shape[1] / aia.shape[1])
aia = zoom(aia, zoom_factors, order=1)


# 顯示 BP3DZ
fig = plt.figure(figsize=(16,8))  # 設定圖大小
ax1 = plt.subplot(121)
sdoaia193 = cm.cmlist["sdoaia193"]
img1 = ax1.imshow(aia, cmap=sdoaia193, vmin=0,vmax=600)
# 繪製等高線，數值 = 100
img1 = ax1.contour(aia, levels=[60], colors='white', linewidths=1)
img1 = ax1.imshow(BP3DZ[:,:,0], cmap="gray",vmin=-200, vmax=200, alpha=0.5)  # 設定 colormap 和數值範圍
ax1.set_title(F"HMI {hmi_time} overlaid on AIA193 {aia_time}")
# ax1.set_xlabel("X-axis")
# ax1.set_ylabel("Y-axis")
cbar1 = plt.colorbar(img1, ax=ax1, fraction=0.046, pad=0.04)
cbar1.ax.tick_params(labelsize=10)
cbar1.ax.set_title('Bz (G)')
k=200
ax2 = plt.subplot(122)
# img2 = ax2.imshow(aia, cmap=sdoaia193, vmin=0,vmax=600)
img2 = ax2.contour(aia, levels=[60], colors='y', linewidths=1)
img2 = ax2.imshow(BP3DZ[:,:,0]/BP3DZ[:,:,k], cmap="gray", vmin=0, vmax=10)  # 設定 colormap 和數值範圍
ax2.set_title(f"Bp Z-dir bottom / layer {k*10}")
# ax2.set_xlabel("X-axis")
# ax2.set_ylabel("Y-axis")
cbar2 = plt.colorbar(img2, ax=ax2, fraction=0.046, pad=0.04)
cbar2.ax.tick_params(labelsize=10)
# cbar2.ax.set_title('Bz (G)')

plt.savefig("B300.jpg",bboxes_inches = 'tight')
ax1.minorticks_on()
ax2.minorticks_on()
plt.show()


############
# # 生成 X, Y, Z 座標
# x = np.linspace(0, BP3DZ.shape[1] - 1, BP3DZ.shape[1])
# y = np.linspace(0, BP3DZ.shape[0] - 1, BP3DZ.shape[0])
# X, Y = np.meshgrid(x, y)

# # 取出 Z 軸的一個切片（例如  層）
# Z = BP3DZ[:, :, 299] 

# # 繪製 3D 圖
# fig = plt.figure(figsize=(10, 8))
# ax = fig.add_subplot(111, projection="3d")

# # 畫 3D 等值面
# ax.plot_surface(X, Y, Z, cmap="coolwarm")

# # 設定標籤
# ax.set_xlabel("X Axis")
# ax.set_ylabel("Y Axis")
# ax.set_zlabel("Magnetic Field Strength (Z)")
# ax.set_title("3D Magnetic Field Visualization")
# plt.show()
#########
# Create grid of points (you may adjust this to fit your actual data coordinates)

#############
# x = np.linspace(0, 1415, BP3DX.shape[0])
# y = np.linspace(0, 1170, BP3DY.shape[1])
# z = np.linspace(0, 300, BP3DZ.shape[2])
# X, Y, Z = np.meshgrid(x, y, z,indexing='ij')
# print(X.shape,Y.shape,Z.shape)

# # Set up the plot
# fig = plt.figure(figsize=(10, 8))
# ax = fig.add_subplot(111, projection='3d')

# # Normalize the magnetic field for color mapping
# field_magnitude = np.sqrt(BP3DX**2 + BP3DY**2 + BP3DZ**2)

# # Quiver plot with color representing the magnetic field magnitude
# ax.quiver(X, Y, Z, BP3DX, BP3DY, BP3DZ, length=0.1, normalize=True, color=field_magnitude, cmap='viridis')

# # Set labels
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Z')




# plt.show()

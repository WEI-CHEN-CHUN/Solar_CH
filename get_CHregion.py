import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
from scipy.ndimage import zoom
import sunpy.visualization.colormaps as cm
import sunpy.map
from datetime import datetime

hmi = pd.read_csv('read_data/hmi_20231104_cropped_sp.csv')
hmi = hmi.to_numpy()
hmi.shape

# 讀取 FITS 檔案
file_path_hmi = "read_data/hmi.m_720s.20231104_000000_TAI.3.magnetogram.fits"
file_path_aia = "read_data/aia.lev1_euv_12s.2023-11-04T000006Z.193.image_lev1.fits"
hmi_map = sunpy.map.Map(file_path_hmi)
aia_map = sunpy.map.Map(file_path_aia)
# 讀取header
hmi_header = hmi_map.wcs.to_header()
aia_header = aia_map.wcs.to_header()
hmi_time = datetime.strptime(hmi_header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S.%f").strftime("%Y-%m-%d %H:%M:%S")
aia_time = datetime.strptime(aia_header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S.%f").strftime("%Y-%m-%d %H:%M:%S")


aia = pd.read_csv('read_data/aia_20231104_cropped_sp.csv')

# 顯示原始數據形狀
print(f"Original hmi shape: {hmi.shape}")
print(f"Original aia shape: {aia.shape}")
zoom_factors = (hmi.shape[0] / aia.shape[0], hmi.shape[1] / aia.shape[1])
aia = zoom(aia, zoom_factors, order=1)
print(f"Resized aia shape: {aia.shape}")

# 設定日冕洞的強度門檻
CH_threshold = 50

# 創建日冕洞的遮罩（低於 threshold 的區域為 True）
CH_mask = aia < CH_threshold

# 取出 HMI 磁場數據中屬於日冕洞的值
hmi_CH_values = hmi[CH_mask]



# 顯示日冕洞區域
fig = plt.figure(figsize=(8,8))  # 設定圖大小
ax1 = plt.subplot(121)
sdoaia193 = cm.cmlist["sdoaia193"]
# 繪製等高線，數值 = 60
img1 = ax1.contour(aia, levels=[50], colors='white', linewidths=1)
img1 = ax1.imshow(aia, cmap=sdoaia193, vmin=0,vmax=600)
ax1.set_title(F"HMI {hmi_time} overlaid on AIA193 {aia_time}")
cbar1 = plt.colorbar(img1, ax=ax1, fraction=0.046, pad=0.04)
cbar1.ax.tick_params(labelsize=10)

ax2 = plt.subplot(122)
# 繪製等高線，數值 = 60
img2 = ax2.contour(aia, levels=[50], colors='white', linewidths=1)
img2 = ax2.imshow(hmi, cmap='gray', vmin=-50,vmax=50)
ax2.set_title(F"HMI {hmi_time} overlaid on AIA193 {aia_time}")
cbar2 = plt.colorbar(img2, ax=ax2, fraction=0.046, pad=0.04)
cbar2.ax.tick_params(labelsize=10)
plt.show()

# 儲存磁場數值陣列
np.save("hmi_CH.npy", hmi_CH_values)

print("提取的磁場數值數量:", hmi_CH_values.shape)

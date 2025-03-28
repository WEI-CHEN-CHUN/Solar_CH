import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
from scipy.ndimage import zoom
import sunpy.visualization.colormaps as cm
import sunpy.map
from datetime import datetime
#　read cropped csv
hmi = pd.read_csv("data_1271/hmi_CH1271_cropped_sp.csv")
aia = pd.read_csv("data_1271/aia_CH1271_cropped_sp.csv")
hmi = hmi.to_numpy() 
print(hmi.shape)
# 讀取 FITS 檔案
file_path_hmi = "./data_1271/hmi.m_45s.20250225_162145_TAI.2.magnetogram.fits"
file_path_aia = "./data_1271/aia.lev1_euv_12s.2025-02-25T162106Z.193.image_lev1.fits"
hmi_map = sunpy.map.Map(file_path_hmi)
aia_map = sunpy.map.Map(file_path_aia)
# 讀取header的時間
hmi_header = hmi_map.wcs.to_header()
aia_header = aia_map.wcs.to_header()
hmi_time = datetime.strptime(hmi_header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S.%f").strftime("%Y-%m-%d %H:%M:%S")
aia_time = datetime.strptime(aia_header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S.%f").strftime("%Y-%m-%d %H:%M:%S")


# 顯示原始數據形狀
print(f"Original hmi shape: {hmi.shape}")
print(f"Original aia shape: {aia.shape}")
# 調整AIA到跟HMI一樣(AIA較小)
zoom_factors = (hmi.shape[0] / aia.shape[0], hmi.shape[1] / aia.shape[1])
aia = zoom(aia, zoom_factors, order=1)
print(f"Resized aia shape: {aia.shape}")

# 設定日冕洞的強度門檻
CH_threshold = 60

# 創建日冕洞的遮罩（低於 threshold 的區域為 True）
CH_mask = aia < CH_threshold

# 取出 HMI 磁場數據中屬於日冕洞的值
hmi_CH_values = hmi[CH_mask]

# 顯示日冕洞區域
fig = plt.figure(figsize=(8,8))  
# 繪製等高線 on aia
ax1 = plt.subplot(121)
img1 = ax1.contour(aia, levels=[CH_threshold], colors='white', linewidths=1)
sdoaia193 = cm.cmlist["sdoaia193"]
img1 = ax1.imshow(aia, cmap=sdoaia193, vmin=0, vmax=600)
ax1.set_title(F"HMI {hmi_time} overlaid on AIA193 {aia_time}")
cbar1 = plt.colorbar(img1, ax=ax1, fraction=0.046, pad=0.04)
cbar1.ax.tick_params(labelsize=10)

# 繪製等高線 on hmi
ax2 = plt.subplot(122)
img2 = ax2.contour(aia, levels=[CH_threshold], colors='#fbfdaf', linewidths=1)
img2 = ax2.imshow(hmi, cmap = 'gray', vmin=-50, vmax=50)
ax2.set_title(F"HMI {hmi_time} overlaid on AIA193 {aia_time}")
cbar2 = plt.colorbar(img2, ax=ax2, fraction=0.046, pad=0.04)
cbar2.ax.tick_params(labelsize=10)
plt.show()

# 繪製直方圖
print("擷取前的磁場數值數量:", hmi_CH_values.shape)
cut_value = 50
hmi_CH_values = hmi_CH_values[np.abs(hmi_CH_values) < cut_value]
print(f"<{cut_value}，擷取後的磁場數值數量:", hmi_CH_values.shape)
plt.figure(figsize=(8, 6))
plt.hist(hmi_CH_values, bins=200, color='yellow', edgecolor='k', alpha=0.7, linewidth=0.2)

# 設定標題與標籤
plt.title("Histogram of HMI Magnetic Field Values in Coronal Hole")
plt.xlabel("Magnetic Field Strength (Gauss)")
plt.ylabel("Frequency")
plt.grid(True, linestyle="--", alpha=0.5)

# 顯示圖表
plt.show()

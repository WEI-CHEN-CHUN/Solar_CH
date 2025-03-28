import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import zoom
import sunpy.visualization.colormaps as cm
import matplotlib.colors as mcolors

# 讀取 CSV 文件到 DataFrame
hmi = pd.read_csv("data_1271/hmi_CH1271_cropped_sp.csv")
aia = pd.read_csv("data_1271/aia_CH1271_cropped_sp.csv")

# 顯示原始數據形狀
print(f"Original hmi shape: {hmi.shape}")
print(f"Original aia shape: {aia.shape}")
zoom_factors = (hmi.shape[0] / aia.shape[0], hmi.shape[1] / aia.shape[1])
aia = zoom(aia, zoom_factors, order=1)
print(f"Resized aia shape: {aia.shape}")
fig = plt.figure(figsize = (10, 10))
ax = plt.subplot(111)
sdoaia193 = cm.cmlist["sdoaia193"]
CH_threshold = 60
img1 = ax.contour(aia, levels=[CH_threshold], colors='white', linewidths=1)
plt.imshow(aia, cmap=sdoaia193, vmin=0, vmax=600)
plt.imshow(hmi, cmap="gray", vmin=-200,vmax=200, alpha=0.6)
plt.colorbar(fraction=0.046, pad=0.04)
plt.xlabel("Solar X (pixel)")
plt.ylabel("Solar Y (pixel)")
plt.title("CH1271 (HMI overlaid on AIA193)")
plt.show()
# # 裁切數據
# cropped_data = df.iloc[4096-2375:4096-960, 1040:2210]
# 960,2375
# 1040,2210
# # 顯示裁切後的數據形狀
# print(f"Cropped shape: {cropped_data.shape}")
# plt.imshow(cropped_data, cmap="gray" ,vmin=-200,vmax=200)
# plt.colorbar(fraction=0.046, pad=0.04)
# plt.show()

# 如果需要保存裁切後的數據，可以寫入新的 CSV
# cropped_data.to_csv('hmi_20231104_cropped.csv', index=False)



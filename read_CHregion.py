import numpy as np
import matplotlib.pyplot as plt

# 讀取 hmi_coronal_hole.npy
hmi_CH_values = np.load("hmi_coronal_hole.npy")
hmi_CH_values = hmi_CH_values[np.abs(hmi_CH_values) < 50]
print("擷取後的磁場數值數量:", hmi_CH_values.shape)
# 繪製直方圖
plt.figure(figsize=(8, 6))
plt.hist(hmi_CH_values, bins=200, color='blue', edgecolor='black', alpha=0.7)

# 設定標題與標籤
plt.title("Histogram of HMI Magnetic Field Values in Coronal Hole")
plt.xlabel("Magnetic Field Strength (Gauss)")
plt.ylabel("Frequency")
plt.grid(True, linestyle="--", alpha=0.5)

# 顯示圖表
plt.show()

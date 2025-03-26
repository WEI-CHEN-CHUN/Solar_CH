from astropy.io import fits
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 讀取 FITS 文件
fits_file = "20231104/hmi.m_720s.20231104_000000_TAI.3.magnetogram.fits"
hdul = fits.open(fits_file)
# 印出主標頭
print(repr(hdul[0].header))
print(repr(hdul[1].header))

# 假設數據在第一個HDU中
data = hdul[1].data  # 這會讀取第二個HDU中的數據（COMPRESSED_IMAGE）
hdul.close()

# 檢查數據格式
print(data.shape)

# 將 FITS 數據轉換為 DataFrame (假設數據是二維數組)
df = pd.DataFrame(data)
print(df)
# 儲存為 CSV 文件
csv_file = "hmi_20231104.csv"
# df.to_csv(csv_file, index=False)

print(f"Data has been successfully saved to {csv_file}")
plt.imshow(df, cmap="gray", vmin=-200,vmax=200)
plt.show()
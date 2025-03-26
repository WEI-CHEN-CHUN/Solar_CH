import sunpy.map
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import pfsspy
import pfsspy.tracing as tracing
from reproject import reproject_interp
from astropy.io import fits
from astropy.time import Time
from sunpy.coordinates.sun import carrington_rotation_number

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
# 設置新的 CEA WCS
wcs_in = hmi_map.wcs
# 計算 Carrington 經度中心
obs_time = Time(hmi_map.date)
# 計算 Carrington 經度中心
carrington_rot = carrington_rotation_number(obs_time)  # Carrington 旋轉號
carrington_longitude = (360 * (carrington_rot % 1))  # 轉換為經度 (0-360)
wcs_out = WCS(naxis=2)
wcs_out.wcs.ctype = ["HPLN-TAN", "HPLT-TAN"]  # 轉換為太陽黃道座標
wcs_out.wcs.cunit = ["deg", "deg"]
wcs_out.wcs.cdelt = [hmi_map.scale[0].value, hmi_map.scale[1].value]  # 保持解析度
wcs_out.wcs.crval = [0, 0]  # 太陽圓盤中心
wcs_out.wcs.crpix = [hmi_map.data.shape[1] / 2, hmi_map.data.shape[0] / 2]
# 使用 reproject 進行投影轉換
output_data, _ = reproject_interp(hmi_map, wcs_out, shape_out=hmi_map.data.shape)
# 儲存轉換後的 FITS 文件
hdu = fits.PrimaryHDU(output_data, header=wcs_out.to_header())
hdu.writeto(F"{file_path_hmi}_cea.fits", overwrite=True)

print("完成 CEA (Carrington) 轉換，準備輸入 pfsspy。")
# 指定像素範圍
# x_range1 = (-1000/3600-aia_header['CRVAL1'])/aia_header['CDELT1']+aia_header['CRPIX1']
# x_range2 = (1000/3600-aia_header['CRVAL1'])/aia_header['CDELT1']+aia_header['CRPIX1']
# y_range1 = (-1000/3600-aia_header['CRVAL2'])/aia_header['CDELT2']+aia_header['CRPIX2']
# y_range2 = (1000/3600-aia_header['CRVAL2'])/aia_header['CDELT2']+aia_header['CRPIX2']
bottom_left = SkyCoord(-600 * u.arcsec, -650 * u.arcsec, frame=aia_map.coordinate_frame)
top_right = SkyCoord(100 * u.arcsec, 200 * u.arcsec, frame=aia_map.coordinate_frame)
# 裁剪地圖
aia_map = aia_map.submap(bottom_left, top_right=top_right)

# 讀取處理後的 CEA FITS
hmi_cea_map = sunpy.map.Map("hmi_cea_pfss.fits")

# 轉換為 pfsspy 格式
input_map = pfsspy.input.MagMap(hmi_cea_map)


# pfss_in = pfsspy.Input(aia_map, nr=25, rss=2.5)
# 設置 PFSS 模型
nr = 35  # 設定網格解析度
pfss_input = pfsspy.input.PFSSInput(input_map, nr)
pfss_model = pfsspy.model.PFSS(pfss_input)
print("PFSS 模型已準備就緒！")




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
plt.show()


import pfsspy
import sunpy.map
from pfsspy import tracing
from astropy.constants import R_sun
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from reproject import reproject_interp
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Circle

# 讀取 HMI 磁場數據
filename = 'hmi.M_720s.20250303_000000_TAI.3.magnetogram.fits'
hmi_map = sunpy.map.Map(filename)
hmi_map = hmi_map.resample([360, 360] * u.pix)
# 轉換為 CEA 投影
cea_header = hmi_map.wcs.to_header()
cea_header['CTYPE1'] = 'CRLN-CEA'
cea_header['CTYPE2'] = 'CRLT-CEA'
cea_header['CDELT1'] = 360 / 360  # 確保解析度適合 PFSS
cea_header['CDELT2'] = 180 / 360 / (np.pi/2)
cea_header['CRPIX1'] = 1800
cea_header['CRPIX2'] = 720
cea_header['CRVAL1'] = 180
cea_header['CRVAL2'] = 0
cea_header['CUNIT1'] = 'deg'
cea_header['CUNIT2'] = 'deg'
cea_wcs = WCS(cea_header)

# 重新投影
hmi_reprojected, _ = reproject_interp(hmi_map, cea_wcs, hmi_map.data.shape)

# 轉換為 SunPy Map
hmi_map_cea = sunpy.map.Map((hmi_reprojected, cea_header))
hmi_map_cea.data[np.isnan(hmi_map_cea.data)] = np.nanmedian(hmi_map_cea.data)

# **修正 NaN 和 無窮大值**
hmi_fixed_data = np.nan_to_num(hmi_reprojected, nan=0.0, posinf=0.0, neginf=0.0)

# 確保數據無 NaN 或 Inf
assert np.isfinite(hmi_fixed_data).all(), "數據仍包含 NaN 或 Inf！"

# **建立新的 SunPy Map**
hmi_map_cea = sunpy.map.Map((hmi_fixed_data, cea_header))

# 進行 PFSS 計算
pfss_input = pfsspy.Input(hmi_map_cea, nr=50, rss=2.5)
pfss_output = pfsspy.pfss(pfss_input)

# 設定磁場線追蹤器
tracer = tracing.FortranTracer()

# 設定種子點（磁場線起點）
r_start = 2.5 * R_sun  # 磁場線起始點 (源表面)
lon = np.linspace(0, 360, 20) * u.deg  # 經度範圍
lat = np.linspace(-90, 90, 10) * u.deg  # 緯度範圍
lon, lat = np.meshgrid(lon, lat, indexing="ij")
lon, lat = lon.ravel(), lat.ravel()

# 生成 SkyCoord 種子點
seeds = SkyCoord(lon, lat, r_start, frame="heliographic_stonyhurst")

# 追蹤磁場線
field_lines = tracer.trace(seeds, pfss_output)

print(f"成功追蹤 {len(field_lines)} 條磁場線！")

# **1. 設定種子點**
nlats, nlons = 20, 40  # 緯度和經度上的網格數量
lons = np.linspace(0, 360, nlons) * u.deg
lats = np.linspace(-60, 60, nlats) * u.deg  # 避免極區

# **2. 創建 2D 網格**
lon_grid, lat_grid = np.meshgrid(lons, lats, indexing="ij")  # 確保形狀正確
radius_grid = np.ones_like(lon_grid) * R_sun  # 半徑 R = R_sun

# **3. 確保 SkyCoord 是 (N, M) 形狀**
seeds = SkyCoord(lon=lon_grid, lat=lat_grid, radius=radius_grid, frame=pfss_output.coordinate_frame)

# **4. 確保 trace() 接受 2D 陣列**
seeds = seeds.reshape(nlats, nlons)  # 這一步確保 (N, M) 維度

# **5. 追蹤磁力線**
tracer = tracing.FortranTracer()
field_lines = tracer.trace(seeds, pfss_output)

fig, ax = plt.subplots(figsize=(8, 8))

# 背景：HMI 磁場圖
ax.imshow(hmi_map_cea.data, origin='lower', extent=[-1, 1, -1, 1], cmap='gray', alpha=0.6)

# **4. 3D -> 2D 投影：XY 平面**
for field_line in field_lines:
    coords = field_line.coords
    r = coords.radius / R_sun  # 磁力線的徑向距離 (標準化)
    
    # 轉換 3D 磁力線到直角坐標 (X, Y, Z)
    x, y, z = coords.spherical.lon.radian, coords.spherical.lat.radian, r - 1  # z 代表高度
    
    # 透視投影 (忽略 Z，壓縮到 XY 平面)
    ax.plot(x, y, color='blue' if field_line.polarity > 0 else 'red', linewidth=0.7)

# **5. 添加太陽圓盤 (R = R_sun)**
sun_circle = Circle((0, 0), 1, edgecolor='black', facecolor='none', linewidth=2)
ax.add_patch(sun_circle)
# **4. 美化圖像**
ax.set_title(filename)
ax.set_xlabel('lat (deg)')
ax.set_ylabel('lon (deg)')
ax.set_xlim([-90, 90])
ax.set_ylim([-90, 90])
plt.grid(zorder=0)
plt.minorticks_on()
plt.show()
import sunpy.map
import pfsspy
from pfsspy import tracing
import numpy as np
import astropy.units as u
from astropy.constants import R_sun
from astropy.coordinates import SkyCoord
from reproject import reproject_interp
from astropy.wcs import WCS
import matplotlib.pyplot as plt

# **1️⃣ 讀取 HMI 磁場數據**
hmi_map = sunpy.map.Map("hmi.M_720s.20250303_000000_TAI.3.magnetogram.fits")

# 重新取樣為 360x360，減少計算量
hmi_map = hmi_map.resample([360, 360] * u.pix)

# **2️⃣ 轉換為 CEA 投影**
cea_header = hmi_map.wcs.to_header()
cea_header.update({
    "CTYPE1": "CRLN-CEA", "CTYPE2": "CRLT-CEA",
    "CDELT1": 360 / 360,  # 經度解析度
    "CDELT2": 180 / 360 / (np.pi / 2),  # 緯度解析度
    "CRPIX1": 1800, "CRPIX2": 720,  # 設定中心
    "CRVAL1": 180, "CRVAL2": 0,
    "CUNIT1": "deg", "CUNIT2": "deg"
})
cea_wcs = WCS(cea_header)

# **重新投影**
hmi_reprojected, _ = reproject_interp(hmi_map, cea_wcs, hmi_map.data.shape)
hmi_map_cea = sunpy.map.Map((np.nan_to_num(hmi_reprojected), cea_header))

# **3️⃣ 執行 PFSS 計算**
pfss_input = pfsspy.Input(hmi_map_cea, nr=50, rss=2.5)  # nr=網格層數, rss=源面高度
pfss_output = pfsspy.pfss(pfss_input)

# **4️⃣ 設定種子點**
nlats, nlons = 20, 40
lons = np.linspace(0, 360, nlons) * u.deg
lats = np.linspace(-60, 60, nlats) * u.deg
lon_grid, lat_grid = np.meshgrid(lons, lats, indexing="ij")
radius_grid = np.full(lon_grid.shape, R_sun.to(u.m).value) * u.m

seeds = SkyCoord(lon=lon_grid, lat=lat_grid, radius=2.5*R_sun, frame="heliographic_stonyhurst")
seeds = seeds.reshape(nlats, nlons)  # **確保是 2D**

# **追蹤磁力線**
tracer = tracing.FortranTracer()
field_lines = tracer.trace(seeds, pfss_output)

# **5️⃣ 繪製結果**
fig, ax = plt.subplots(figsize=(8, 8))

# **繪製磁力線**
for field_line in field_lines:
    xyz = field_line.coords
    x, y = xyz.x.to(u.R_sun).value, xyz.y.to(u.R_sun).value
    ax.plot(x, y, color="cyan", linewidth=0.5)

ax.set_title("PFSS Magnetic Field Lines on AIA 193")
plt.show()

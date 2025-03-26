import pfsspy
import sunpy.map
from pfsspy import tracing
from astropy.constants import R_sun
import numpy as np
import astropy.units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# 讀取 HMI 磁場數據
filename = 'hmi.synoptic_mr_720s.2294.synopMr.fits'
hmi_map = sunpy.map.Map(filename)
header = hmi_map.wcs.to_header()
print('Data shape: ', hmi_map.data.shape)
for key, value in header.items():
    print(f"{key}: {value}")
# hmi_map = hmi_map.resample([360, 180] * u.pix)
hmi_map.data[np.isnan(hmi_map.data)] = np.nanmedian(hmi_map.data)

# 設定 PFSS 模型
pfss_input = pfsspy.Input(hmi_map, nr=25, rss=3.5)  # 設定源表面為 2.5R_sun
pfss_output = pfsspy.pfss(pfss_input)
# 繪製磁場線
r = 2.5 * R_sun
lat = np.linspace(-np.pi / 2, np.pi / 2, 15, endpoint=False)
lon = np.linspace(0, 2 * np.pi, 15, endpoint=False)
lat, lon = np.meshgrid(lat, lon, indexing='ij')
lat, lon = lat.ravel() * u.rad, lon.ravel() * u.rad

# 生成 SkyCoord 物件，提供多個種子點
seeds = SkyCoord(lon, lat, r, frame="heliographic_stonyhurst")
tracer = tracing.FortranTracer()
field_lines = tracer.trace(seeds, pfss_output)
print(type(field_lines))




# 創建 3D 圖像
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')

# 繪製每條磁場線
for field_line in field_lines:
    color = {0: 'black', -1: 'tab:red', 1: 'tab:blue'}.get(field_line.polarity)
    coords = field_line.coords
    exp_fact = field_line.expansion_factor
    coords.representation_type = 'spherical'
    r_vals = coords.radius / R_sun  # 轉換為以 R_sun 為單位
    theta_vals = np.pi/2 - coords.lat.to(u.rad).value  # 轉換為 θ（從 z 軸計算）
    phi_vals = coords.lon.to(u.rad).value  # 經度

    # --- ✨ 手動轉換球座標為笛卡爾座標 ---
    x_vals = r_vals * np.sin(theta_vals) * np.cos(phi_vals)
    y_vals = r_vals * np.sin(theta_vals) * np.sin(phi_vals)
    z_vals = r_vals * np.cos(theta_vals)
    
    ax.plot(x_vals, y_vals, z_vals, color=color, linewidth=0.7)  # 使用轉換後的座標繪圖
    # 🎯 **在磁場線中間標記 Expansion Factor 數值**

    if len(x_vals) == 0 or exp_fact>18 or np.isnan(exp_fact):
        continue
    mid_idx = len(x_vals)-1  # 取得磁場線的索引
    if r_vals[0]>2:
        mid_idx = 0
    # 🔴 **跳過 x < 0 的磁場線**
    if np.any(x_vals < 0):  
        continue
    # f_s_value = exp_fact.flat[i]  # 提取對應的 f_s 值
    ax.text(x_vals[mid_idx], y_vals[mid_idx], z_vals[mid_idx], 
            f"{exp_fact:.1f}", color='black', fontsize=10)


# --- 🔴 繪製 R = R_sun 的球面 ---
u_vals = np.linspace(0, 2 * np.pi, 24)  # 經度範圍 (0°~360°)
v_vals = np.linspace(0, np.pi, 12)  # 緯度範圍 (0°~180°)
U, V = np.meshgrid(u_vals, v_vals)

# 球面轉換為笛卡爾座標
X_sphere = np.sin(V) * np.cos(U)  # 球面 X 座標
Y_sphere = np.sin(V) * np.sin(U)  # 球面 Y 座標
Z_sphere = np.cos(V)  # 球面 Z 座標

# 繪製赤道
u_vals = np.linspace(0, 2 * np.pi, 24)  # 經度範圍 (0°~360°)
v_vals = np.pi/2  # 緯度範圍 (90°)
U, V = np.meshgrid(u_vals, v_vals)

# 球面轉換為笛卡爾座標
X_sphere_0 = np.sin(V) * np.cos(U)  # 球面 X 座標
Y_sphere_0 = np.sin(V) * np.sin(U)  # 球面 Y 座標
Z_sphere_0 = np.cos(V)  # 球面 Z 座標

# 繪製0度經線
u_vals = 0  # 0度經度
v_vals = np.linspace(0, np.pi, 12)  # 緯度範圍 (0°~180°)
U, V = np.meshgrid(u_vals, v_vals)

# 球面轉換為笛卡爾座標
X_sphere_1 = np.sin(V) * np.cos(U)  # 球面 X 座標
Y_sphere_1 = np.sin(V) * np.sin(U)  # 球面 Y 座標
Z_sphere_1 = np.cos(V)  # 球面 Z 座標


ax.plot_surface(X_sphere, Y_sphere, Z_sphere, color="y", alpha=0.3)  # 畫出 R=R_sun 的球面
ax.plot_wireframe(X_sphere, Y_sphere, Z_sphere, color="k", alpha=0.2)  # 畫出 R=R_sun 的球面
ax.plot_wireframe(X_sphere_0, Y_sphere_0, Z_sphere_0, color="r", alpha=1, linewidth=2)  # 畫出 赤道
ax.plot_wireframe(X_sphere_1, Y_sphere_1, Z_sphere_1, color="r", alpha=1, linewidth=2)  # 畫出 0度經線
# ax.plot_wireframe(2.5*X_sphere, 2.5*Y_sphere, 2.5*Z_sphere, color="y", alpha=0.3)  # 畫出 R=2.5*R_sun 的球面
ax.text(X_sphere_0[0][0], Y_sphere_0[0][0], Z_sphere_0[0][0], 
            f"0°", color='black', fontsize=15)
# 設定標籤
# ax.set_xlabel("X ($R_☉$)")
# ax.set_ylabel("Y ($R_☉$)")
# ax.set_zlabel("Z ($R_☉$)")
ax.set_title("PFSS Solution " + filename)
# ax.grid()
ax.view_init(elev=0, azim=0)
# plt.axis('equal')
plt.axis('off')
plt.show()

# tracer.plot_field_lines(field_lines)
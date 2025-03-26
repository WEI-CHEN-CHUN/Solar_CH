from scipy.io import readsav
import matplotlib.colors as mcolors
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pyvista as pv

# 讀取 IDL 的 .sav 檔案
data = readsav("C:/Users/chjan/HMI_output_20231104.sav")

# # 列出所有變數名稱
# print(data.keys())

# 取出 BP3DZ 資料
BP3DX = data["BP3DX"].T
BP3DY = data["BP3DY"].T
BP3DZ = data["BP3DZ"].T
print("shape:", BP3DX.shape, BP3DY.shape, BP3DZ.shape)



import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# 假設 Bx, By, Bz, X, Y, Z 為已知的磁場向量場
Bx, By, Bz = BP3DX[0:10,0:10,0:10], BP3DY[0:10,0:10,0:10], BP3DZ[0:10,0:10,0:10]
print(Bx)
X, Y, Z = np.meshgrid(np.linspace(-5, 5, 10), np.linspace(-5, 5, 10), np.linspace(-5, 5, 10))
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# 選擇起點
start_points = [(-4, -4, -4), (4, 4, 4), (-4, 4, -4)]

# 磁力線追蹤函數（簡單的 Euler 方式）
def trace_fieldline(sp, Bx, By, Bz, X, Y, Z, step=0.2, n_steps=100):
    x, y, z = [sp[0]], [sp[1]], [sp[2]]
    for _ in range(n_steps):
        i, j, k = np.searchsorted(X[:, 0, 0], x[-1]), np.searchsorted(Y[0, :, 0], y[-1]), np.searchsorted(Z[0, 0, :], z[-1])
        if i >= 10 or j >= 10 or k >= 10 or i < 0 or j < 0 or k < 0:
            break
        bx, by, bz = Bx[i, j, k], By[i, j, k], Bz[i, j, k]
        norm = np.sqrt(bx**2 + by**2 + bz**2)
        if norm == 0:
            break
        bx, by, bz = bx / norm, by / norm, bz / norm
        x.append(x[-1] + step * bx)
        y.append(y[-1] + step * by)
        z.append(z[-1] + step * bz)
    return np.array(x), np.array(y), np.array(z)

# 繪製磁場線
for sp in start_points:
    x_traj, y_traj, z_traj = trace_fieldline(sp, Bx, By, Bz, X, Y, Z)
    print(x_traj, y_traj, z_traj)
    ax.plot3D(x_traj, y_traj, z_traj, color='r')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('3D 磁場線')

plt.show()

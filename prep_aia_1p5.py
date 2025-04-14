import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
from aiapy.calibrate import register, update_pointing
from aiapy.calibrate.util import get_pointing_table

file_path_aia = "./data_1271/aia.lev1_euv_12s.2025-02-25T162106Z.193.image_lev1.fits"
aia_map = sunpy.map.Map(file_path_aia)

# Make range wide enough to get closest 3-hour pointing
pointing_table = get_pointing_table("JSOC", 
                  time_range=(aia_map.date - 12 * u.h, aia_map.date + 12 * u.h))
aia_map_updated_pointing = update_pointing(aia_map, pointing_table=pointing_table)
print(aia_map_updated_pointing.scale)
print(aia_map_updated_pointing.rotation_matrix)
aia_map_registered = register(aia_map_updated_pointing)
print(aia_map_registered.scale)
print(aia_map_registered.rotation_matrix)
fig = plt.figure()
ax = fig.add_subplot(121,projection=aia_map_registered)
ax1 = fig.add_subplot(122,projection=aia_map)
aia_map_registered.plot(axes=ax)
aia_map.plot(axes=ax1)
aia_map_registered.save(file_path_aia.replace("lev1", "lev1p5"), overwrite=True)
# plt.show()
import sunpy.map
from sunpy.data.sample import AIA_193_IMAGE
from astropy.wcs import WCS

m = sunpy.map.Map(AIA_193_IMAGE)
print(m.wcs.wcs.aux)

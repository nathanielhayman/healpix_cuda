import healpy as hp
import numpy as np

mp = hp.read_map("../src/output.fits")

hp.mollzoom(mp, title="Mollview image RING")
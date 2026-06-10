from astroquery.nist import Nist
import astropy.units as u


table = Nist.query(500 * u.nm, 600 * u.nm, linename="Ba I")

print(table.colnames)
print(table[:5])
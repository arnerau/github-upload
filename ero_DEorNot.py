#! /Users/arau/opt/anaconda3/bin/python3

###############
#   eroDEorNot.py : command line script that checks for a given set of coordinates (in equatorial coordinates)
#                      - if it is in the German or Russian part of the sky
#                      - when it was or will be covered during the All-sky survey.
#                   The script uses the survey simulations downloaded from: 
#                      https://erosita.hs.uni-hamburg.de/hserosita/static/erass1_eclscans_182d.fits
#                      https://erosita.hs.uni-hamburg.de/hserosita/static/erass2_eclscans_sim.fits 
#
#                 : v1.0 (ArneRau / 2020-04-03)
#                 : v1.1 (ArneRau / 2020-06-17)
#                        - eRASS2 added 
#  
#   usage: python eroDEornot.py HH:MM:SS.S DD:MM:SS.S
#
##############


import sys
from astropy import units as u
from astropy.coordinates import SkyCoord 
from astropy.io import fits

# specify location of survey simulation
surveyplanning = '/Users/arau/bin/python_scripts/erass1_eclscans_182d.fits'
surveyplanning_eRASS2 = '/Users/arau/bin/python_scripts/erass2_eclscans_sim.fits'

# read coordinates given in command line
ra = sys.argv[1]
dec = sys.argv[2] 

# transform input coordinates into SkyCoord 
c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))

# print (l,b) to terminal and DE/RU
print('')
print('Input Equatorial Coordinates:')
print('RA(J2000) = %5.4fdeg' % c.ra.degree)
print('Dec(J2000) = %5.4fdeg\n' % c.dec.degree)
print('Transformed to Galactic Coordinates:')
print('l = %5.4fdeg' % c.galactic.l.degree)
print('b = %5.4fdeg\n' % c.galactic.b.degree)
if c.galactic.l.degree >= 180:
	print('--> This source is in DE')
else:
	print('--> This source is in RU')
print('')

# transform into Ecliptic coordinates
print('Transformed to Ecliptic Coordinates:')
print('lon = %5.4fdeg' % c.barycentricmeanecliptic.lon.degree)
print('lat = %5.4fdeg\n' % c.barycentricmeanecliptic.lat.degree)
lon = c.barycentricmeanecliptic.lon.degree

print('Position observed by eROSITA on (approx.):')

# read in survey planning file for eRASS:1
hdul = fits.open(surveyplanning)
data = hdul[1].data

# search for match between coordinates and survey simulation.
for row in range(len(data)):
	if -1 < lon - data[row][4] < 1 or (lon < 1 and -359 < lon - data[row][4] < 1) or (lon > 359 and 358 < lon -data[row][4] < 359):
		print(data[row][0])	
	if -1 < lon - data[row][5] < 1 or (lon < 1 and -359 < lon - data[row][5] < 1) or (lon > 359 and 358 < lon -data[row][5] < 359):
		print(data[row][0])	

# read in survey planning file for eRASS:2
hdul = fits.open(surveyplanning_eRASS2)
data = hdul[1].data

# search for match between coordinates and survey simulation.
for row in range(len(data)):
	if -1 < lon - data[row][4] < 1 or (lon < 1 and -359 < lon - data[row][4] < 1) or (lon > 359 and 358 < lon -data[row][4] < 359):
		print(data[row][0])	
	if -1 < lon - data[row][5] < 1 or (lon < 1 and -359 < lon - data[row][5] < 1) or (lon > 359 and 358 < lon -data[row][5] < 359):
		print(data[row][0])	



hdul.close()

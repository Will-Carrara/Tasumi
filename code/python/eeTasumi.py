#                            Will Carrara || 12/07/18 || Python3.7 || Tasumi
#                                           NASA AMES RESEARCH

# Equations provided by: At-Surface Reflectance and Albedo from Satellite for Operational Calculation of
# Land Surface Energy Balance.

# -Journal of Hydrological Engineering (Feb. 2008)
# Masahiro Tasumi, Richard G. Allen, Ricardo Trezza

import ee
import math
import time
import progressbar

# This program presents a rapid, operational method for estimating at-surface reflectance. This is most
# applicable for Landsat satellite sensors for cloud-free, low-haze, conditions with sensor view angles
# less than 20 degrees. The data supplied to this algorithm is from Google Earth Engine.

# initialize earth engine
ee.Initialize()

# Satellite Key
# 2 -> Sentinel 2a
# 5 -> Landsat 5
# 7 -> Landsat 7
# 8 -> Landsat 8

satellite = 2
state = 'California'
toDrive = True


def initialize():
	# United States polygon region filtered by selected state
	us_fc = ee.FeatureCollection('ft:1fRY18cjsHzDgGiJiS2nnpUU3v9JPDc2HNaR7Xk8').filter(ee.Filter.eq('Name', state))

	#MGRS_TILE = '11SKV'  # -> 42,35
	#MGRS_TILE = '10SEJ'  # -> 43,34
	MGRS_TILE = '10SFG'   # -> 44,33

	# Sentinel 2a
	if (satellite == 2):
		refl_toa_coll = ee.ImageCollection('COPERNICUS/S2').filterDate('2018-07-01', '2018-07-31').filterBounds(us_fc).filterMetadata('MGRS_TILE', 'equals', MGRS_TILE)

	# Landsat 8
	elif (satellite == 8):
		refl_toa_coll = ee.ImageCollection('LANDSAT/LC08/C01/T1_TOA').filterDate('2016-07-01', '2016-07-31').filterBounds(us_fc).filter(ee.Filter.Or( ee.Filter.And(ee.Filter.eq('WRS_PATH', 42), ee.Filter.eq('WRS_ROW', 35)), ee.Filter.And(ee.Filter.eq('WRS_PATH', 43), ee.Filter.eq('WRS_ROW', 34)), ee.Filter.And(ee.Filter.eq('WRS_PATH', 44), ee.Filter.eq('WRS_ROW', 33))))

	# Landsat 7
	elif (satellite == 7):
		refl_toa_coll = ee.ImageCollection('LANDSAT/LE07/C01/T1_TOA').filterDate('2013-07-01', '2013-07-31').filterBounds(us_fc).filter(ee.Filter.Or(ee.Filter.And(ee.Filter.eq('WRS_PATH', 42), ee.Filter.eq('WRS_ROW', 35)), ee.Filter.And(ee.Filter.eq('WRS_PATH', 43), ee.Filter.eq('WRS_ROW', 34)),ee.Filter.And(ee.Filter.eq('WRS_PATH', 44), ee.Filter.eq('WRS_ROW', 33))))

	# Landsat 5
	elif (satellite == 5):
		refl_toa_coll = ee.ImageCollection('LANDSAT/LT05/C01/T1_TOA').filterDate('2010-07-01', '2010-07-31').filterBounds(us_fc).filter(ee.Filter.Or( ee.Filter.And(ee.Filter.eq('WRS_PATH', 42), ee.Filter.eq('WRS_ROW', 35)), ee.Filter.And(ee.Filter.eq('WRS_PATH', 43), ee.Filter.eq('WRS_ROW', 34)), ee.Filter.And(ee.Filter.eq('WRS_PATH', 44), ee.Filter.eq('WRS_ROW', 33))))

	else:
		raise ValueError('Operational for only Landsat 5, 7, 8, and Sentinel 2a.')

	# NLDAS hourly collection
	nldas_coll = ee.ImageCollection('NASA/NLDAS/FORA0125_H002')

	# for each Landsat scene, get nearest NLDAS hourly image before and after scene time
	nldas_filter = ee.Filter.maxDifference(24*60*60*1000, 'system:time_start', None, 'system:time_start', None)

	nldas_prev_filter = ee.Filter.And(nldas_filter, ee.Filter.greaterThan('system:time_start', None, 'system:time_start', None))
	nldas_next_filter = ee.Filter.And(nldas_filter, ee.Filter.lessThan('system:time_start', None, 'system:time_start', None))

	refl_toa_coll = ee.ImageCollection(ee.Join.saveBest('nldas_prev_match', 'nldas_prev_metric').apply(refl_toa_coll, nldas_coll, nldas_prev_filter))

	refl_toa_coll = ee.ImageCollection(ee.Join.saveBest('nldas_next_match', 'nldas_next_metric').apply(refl_toa_coll, nldas_coll, nldas_next_filter))

	'''
	surf_coll = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR') \
		.filterDate('2010-07-01', '2010-07-31').filterBounds(us_fc) \
			.filter(ee.Filter.Or(
    		ee.Filter.And(ee.Filter.eq('WRS_PATH', 42), ee.Filter.eq('WRS_ROW', 35)), \
    		ee.Filter.And(ee.Filter.eq('WRS_PATH', 43), ee.Filter.eq('WRS_ROW', 34)), \
    		ee.Filter.And(ee.Filter.eq('WRS_PATH', 44), ee.Filter.eq('WRS_ROW', 33))))
	'''

	'''
	surf_coll = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR') \
		.filterDate('2013-07-01', '2013-07-31').filterBounds(us_fc) \
			.filter(ee.Filter.Or(
    		ee.Filter.And(ee.Filter.eq('WRS_PATH', 42), ee.Filter.eq('WRS_ROW', 35)), \
    		ee.Filter.And(ee.Filter.eq('WRS_PATH', 43), ee.Filter.eq('WRS_ROW', 34)), \
    		ee.Filter.And(ee.Filter.eq('WRS_PATH', 44), ee.Filter.eq('WRS_ROW', 33))))
	'''

	'''
	surf_coll = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR') \
		.filterDate('2016-07-01', '2016-07-31').filterBounds(us_fc) \
			.filter(ee.Filter.Or(
    		ee.Filter.And(ee.Filter.eq('WRS_PATH', 42), ee.Filter.eq('WRS_ROW', 35)), \
    		ee.Filter.And(ee.Filter.eq('WRS_PATH', 43), ee.Filter.eq('WRS_ROW', 34)), \
    		ee.Filter.And(ee.Filter.eq('WRS_PATH', 44), ee.Filter.eq('WRS_ROW', 33))))
	'''

	return refl_toa_coll


def tasumi_reflectance_func(refl_toa_image):

	scene_time = ee.Number(refl_toa_image.get('system:time_start'))
	scene_date = ee.Algorithms.Date(refl_toa_image.get('system:time_start'))
	doy = ee.Number(scene_date.getRelative('day', 'year')).add(1).double()
	hour = ee.Number(scene_date.getFraction('day')).multiply(24)

	# national elevation dataset
	elev = ee.Image('USGS/NED')

	lat = ee.Image.pixelLonLat().select(['latitude']).multiply(math.pi/180)
	lon = ee.Image.pixelLonLat().select(['longitude']).multiply(math.pi/180)

	terrain = ee.call('Terrain', elev)

	slope = terrain.select(['slope']).multiply(math.pi/180)
	aspect = terrain.select(['aspect']).multiply(math.pi/180).subtract(math.pi)

	# mean atmospheric pressure [kPa]
	pair = elev.expression('101.3 * pow((293 - 0.0065 * b()) / 293, 5.26)')

	# interpolate NLDAS image at Landsat scene time
	nldas_prev_image = ee.Image(refl_toa_image.get('nldas_prev_match'))
	nldas_next_image = ee.Image(refl_toa_image.get('nldas_next_match'))
	nldas_prev_time = ee.Number(nldas_prev_image.get('system:time_start'))
	nldas_next_time = ee.Number(nldas_next_image.get('system:time_start'))

	# calculate time ratio of Landsat image between NLDAS images
	time_ratio_image = ee.Image.constant(scene_time.subtract(nldas_prev_time).divide(nldas_next_time.subtract(nldas_prev_time)))

	# interpolate NLDAS values at Landsat image time
	nldas_image = nldas_next_image.subtract(nldas_prev_image).multiply(time_ratio_image).add(nldas_prev_image).set({'system:time_start':scene_time, 'system:time_end': scene_time})

	# calculate vapor pressure from NLDAS hourly specific humidity
	q = nldas_image.select(['specific_humidity'])
	ea = pair.multiply(q).divide(q.multiply(0.378).add(0.622))
	# atmospheric precipitable water
	w = pair.multiply(0.14).multiply(ea).add(2.1)


	def cos_theta_mountain_func(acq_doy, acq_time, lat, lon, slope, aspect):
	    delta = doy.multiply(2 * math.pi / 365).subtract(1.39435).sin().multiply(0.40928)

	    b = doy.subtract(81).multiply(2 * math.pi / 364)

	    sc = b.multiply(2).sin().multiply(0.1645).subtract(b.cos().multiply(0.1255)).subtract(b.sin().multiply(0.025))

	    solar_time = lon.expression('t + (lon * 12 / pi) + sc',{'pi':math.pi,'t':ee.Image.constant(acq_time), 'lon':lon, 'sc':ee.Image.constant(sc)})

	    omega = solar_time.subtract(12).multiply(math.pi / 12)
	    slope_c = slope.cos()
	    slope_s = slope.sin()

	    cos_theta = lat.expression(('(sin(lat) * slope_c * delta_s) - (cos(lat) * slope_s * cos(aspect) * delta_s) + (cos(lat) * slope_c * cos(omega) * delta_c) + (sin(lat) * slope_s * cos(aspect) * cos(omega) * delta_c) + (sin(aspect) * slope_s * sin(omega) * delta_c)'),{'lat':lat, 'aspect': aspect,'slope_c':slope_c,'slope_s':slope_s, 'omega':omega, 'delta_c':ee.Image.constant(delta.cos()),'delta_s':ee.Image.constant(delta.sin())})

	    cos_theta = cos_theta.divide(slope_c).max(ee.Image.constant(0.1))

	    return cos_theta


	cos_theta = cos_theta_mountain_func(doy, hour, lat, lon, slope, aspect).select([0], ['cos_theta'])

	if (satellite == 8 or satellite == 2):
		# Calibrated Landsat Constants LS8 _________________________________________________
		# Coeff        Band1      Band2      Band3      Band4     Band5/8      Band7	   |
		c1 = ee.Image([0.987,     2.148,     0.942,     0.248,    0.260,     0.315    ]) # |
		c2 = ee.Image([-0.000727, -0.000199, -0.000261, -0.00041, -0.001084, -0.000975]) # |
		c3 = ee.Image([0.000037,  0.000058,  0.000406,  0.000563, 0.000675,  0.004012 ]) # |
		c4 = ee.Image([0.0869,    0.0464,    0.0928,    0.2256,   0.0632,    0.0116   ]) # |
		c5 = ee.Image([0.0788,    -1.0962,   0.1125,    0.7991,   0.7549,    0.6906   ]) # |
		cb = ee.Image([0.640,     0.310,     0.286,     0.189,    0.274,     -0.186   ]) # |
		# _________________________________________________________________________________|

	elif (satellite == 5 or satellite == 7):
		# Calibrated Landsat Constants LS5 ____________________________________________
		# Coeff        Band1     Band2     Band3     Band4     Band5     Band7	      |
		c1 = ee.Image([0.987,    2.319,    0.951,    0.375,    0.234,    0.365   ]) # |
		c2 = ee.Image([-0.00071, -0.00016, -0.00033, -0.00048, -0.00101, -0.00097]) # |
		c3 = ee.Image([0.000036, 0.000105, 0.00028,  0.005018, 0.004336, 0.004296]) # |
		c4 = ee.Image([0.0880,   0.0437,   0.0875,   0.1355,   0.0560,   0.0155  ]) # |
		c5 = ee.Image([0.0789,   -1.2697,  0.1014,   0.6621,   0.7757,   0.639   ]) # |
		cb = ee.Image([0.640,    0.310,    0.286,    0.189,    0.274,    -0.186  ]) # |
		# ____________________________________________________________________________|

	else:
		raise ValueError('Operational for only Landsat 5, 7, 8, and Sentinel 2a.')

	tau_in = pair.multiply(c2).divide(cos_theta).subtract(w.multiply(c3).add(c4).divide(cos_theta)).exp().multiply(c1).add(c5)

	tau_out = pair.multiply(c2).subtract(w.multiply(c3).add(c4)).exp().multiply(c1).add(c5)

	# calculate at-surface reflectance for a collection
	if (satellite == 5 or satellite == 7 or satellite == 8):
		refl_sur_image = refl_toa_image.select(['B1', 'B2', 'B3', 'B4', 'B5', 'B7']).expression('(b() + cb * (tau_in - 1)) / (tau_in * tau_out)', {'cb':cb, 'tau_in':tau_in, 'tau_out':tau_out})

	elif (satellite == 2):
		refl_sur_image = refl_toa_image.select(['B1', 'B2', 'B3', 'B4', 'B5', 'B8']).expression('(b() + cb * (tau_in - 1)) / (tau_in * tau_out)', {'cb':cb, 'tau_in':tau_in, 'tau_out':tau_out})

	else:
		raise ValueError('Operational for only Landsat 5, 7, 8, and Sentinel 2a.')


	return refl_sur_image


def upload(collection):

	collectionList = collection.toList(collection.size())
	collectionSize = collection.size().getInfo()

	pbar = progressbar.ProgressBar(maxval=collectionSize).start()

	if (satellite == 8):
		bnds = ['B4', 'B3', 'B2', 'B5']

	elif (satellite == 5 or satellite == 7):
		bnds = ['B3', 'B2', 'B1', 'B4']

	elif (satellite == 2):
		bnds = ['B4', 'B3', 'B2', 'B8']

	else:
		raise ValueError('Operational for only Landsat 5, 7, 8, and Sentinel 2a.')

	# export to drive
	for i in range(collectionSize):
		name = ee.Image(collectionList.get(i)).get('system:index').getInfo()

		if (satellite == 2):
			myTry = ee.batch.Export.image.toDrive(
			image = ee.Image(collectionList.get(i)).select(bnds),
			description = 'Sent-2A July' + str(i+1),
			scale = 30,
			fileNamePrefix = 'T' + str(i+1) + '_' + name,
			folder = 'Sentinel 2A July')

		else:
			myTry = ee.batch.Export.image.toDrive(
			image = ee.Image(collectionList.get(i)).select(bnds).expression('b() * 10000'),
			description = 'Landsat July' + str(i+1),
			scale = 30,
			fileNamePrefix = 'T' + str(i+1) + '_' + name,
			folder = 'Landsat 8 July')


		myTry.start()
		pbar.update(i+1)

	#.expression('b() * 10000') (this is for only L7,L8,L5)

def main():

	refl_toa_coll = initialize()

	print("\nEstimating at-Surface Reflectance")
	bar = progressbar.ProgressBar()
	for i in bar(range(100)):
		refl_sur_coll = ee.ImageCollection(refl_toa_coll.map(tasumi_reflectance_func))

	if (toDrive == True):
		print("\nUploading to Google Drive")
		upload(refl_sur_coll)

main()

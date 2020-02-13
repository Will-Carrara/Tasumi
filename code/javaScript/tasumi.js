//                           Will Carrara || 12/07/18 || Python3.7 || Tasumi
//                                          NASA AMES RESEARCH

// Equations provided by: At-Surface Reflectance and Albedo from Satellite for Operational Calculation of
// Land Surface Energy Balance.

// -Journal of Hydrological Engineering (Feb. 2008)
// Masahiro Tasumi, Richard G. Allen, Ricardo Trezza

// This program presents a rapid, operational method for estimating at-surface reflectance. This is most
// applicable for Landsat satellite sensors for cloud-free, low-haze, conditions with sensor view angles
// less than 20 degrees. The data supplied to this algorithm is from Google Earth Engine.


// United States polygon region filtered by selected states
var us_fc = ee.FeatureCollection('ft:1fRY18cjsHzDgGiJiS2nnpUU3v9JPDc2HNaR7Xk8')
  .filter(ee.Filter.or(ee.Filter.eq('Name', 'Nevada'), ee.Filter.eq('Name', 'California')));

// Landsat collection
var refl_toa_coll = ee.ImageCollection('LANDSAT/LE07/C01/T1_TOA')
  .filterDate("2017-01-01", "2017-12-31").filterBounds(us_fc);

// NLDAS hourly collection
var nldas_coll = ee.ImageCollection('NASA/NLDAS/FORA0125_H002');

// for each Landsat scene, get nearest NLDAS hourly image before and after scene time
var nldas_filter = ee.Filter.maxDifference(
	1000*60*60*4, "system:time_start", null, "system:time_start", null);
var nldas_prev_filter = ee.Filter.and(
  nldas_filter, ee.Filter.greaterThan("system:time_start", null, "system:time_start", null));

var nldas_next_filter = ee.Filter.and(
  nldas_filter,ee.Filter.lessThan("system:time_start", null, "system:time_start", null));

var refl_toa_coll = ee.ImageCollection(ee.Join.saveBest('nldas_prev_match', 'nldas_prev_metric')
  .apply(refl_toa_coll, nldas_coll, nldas_prev_filter));

var refl_toa_coll = ee.ImageCollection(ee.Join.saveBest('nldas_next_match', 'nldas_next_metric')
  .apply(refl_toa_coll, nldas_coll, nldas_next_filter));


var tasumi_reflectance_func = function (refl_toa_image) {
  var scene_time = ee.Number(refl_toa_image.get("system:time_start"));
  var scene_date = ee.Algorithms.Date(refl_toa_image.get("system:time_start"));
  var doy = ee.Number(scene_date.getRelative('day', 'year')).add(1).double();
  var hour = ee.Number(scene_date.getFraction('day')).multiply(24);

  // national elevation dataset
  var elev = ee.Image("USGS/NED");

  var lat = ee.Image.pixelLonLat().select(['latitude']).multiply(Math.PI/180);
  var lon = ee.Image.pixelLonLat().select(['longitude']).multiply(Math.PI/180);

  var terrain = ee.call('Terrain', elev);

  var slope = terrain.select(["slope"]).multiply(Math.PI/180);
  var aspect = terrain.select(["aspect"]).multiply(Math.PI/180).subtract(Math.PI);

  // mean atmospheric pressure [kPa]
  var pair = elev.expression('101.3 * pow((293 - 0.0065 * b()) / 293, 5.26)');

  // interpolate NLDAS image at Landsat scene time
  var nldas_prev_image = ee.Image(refl_toa_image.get("nldas_prev_match"));
  var nldas_next_image = ee.Image(refl_toa_image.get("nldas_next_match"));
  var nldas_prev_time = ee.Number(nldas_prev_image.get("system:time_start"));
  var nldas_next_time = ee.Number(nldas_next_image.get("system:time_start"));

  // calculate time ratio of Landsat image between NLDAS images
  var time_ratio_image = ee.Image.constant(scene_time.subtract(nldas_prev_time)
	.divide(nldas_next_time.subtract(nldas_prev_time)));

  // interpolate NLDAS values at Landsat image time
  var nldas_image = nldas_next_image.subtract(nldas_prev_image)
	.multiply(time_ratio_image).add(nldas_prev_image)
	.set({"system:time_start": scene_time, "system:time_end": scene_time});

  // calculate vapor pressure from NLDAS hourly specific humidity
  var q = nldas_image.select(["specific_humidity"]);  // kg/kg
  var ea = pair.multiply(q).divide(q.multiply(0.378).add(0.622));

  // Atmospheric precipitable water
  var w = pair.multiply(0.14).multiply(ea).add(2.1);


  function cos_theta_mountain_func(acq_doy, acq_time, lat, lon, slope, aspect) {
	var delta = doy.multiply(2 * Math.PI / 365).subtract(1.39435).sin().multiply(0.40928);

	var b = doy.subtract(81).multiply(2 * Math.PI / 364);

	var sc = b.multiply(2).sin().multiply(0.1645)
	  .subtract(b.cos().multiply(0.1255))
	  .subtract(b.sin().multiply(0.025));

	var solar_time = lon.expression(
	  't + (lon * 12 / pi) + sc',
	  {'pi':Math.PI, 't':ee.Image.constant(acq_time),
	   'lon':lon, 'sc':ee.Image.constant(sc)});

	var omega = solar_time.subtract(12).multiply(Math.PI / 12);

	var slope_c = slope.cos();
	var slope_s = slope.sin();

	var cos_theta = lat.expression(
		('(sin(lat) * slope_c * delta_s) - ' +
		 '(cos(lat) * slope_s * cos(aspect) * delta_s) + ' +
		 '(cos(lat) * slope_c * cos(omega) * delta_c) + ' +
		 '(sin(lat) * slope_s * cos(aspect) * cos(omega) * delta_c) + ' +
		 '(sin(aspect) * slope_s * sin(omega) * delta_c)'),
		{'lat':lat, 'aspect':aspect,
		 'slope_c':slope_c, 'slope_s':slope_s, 'omega':omega,
		 'delta_c':ee.Image.constant(delta.cos()),
		 'delta_s':ee.Image.constant(delta.sin())});

	cos_theta = cos_theta.divide(slope_c).max(ee.Image.constant(0.1));

	return cos_theta;
};


  var cos_theta = cos_theta_mountain_func(doy, hour, lat, lon, slope, aspect).select([0], ['cos_theta']);

  // calibrated Landsat constants LS5 && LS7
  var c1 = ee.Image([0.987, 2.319, 0.951, 0.375, 0.234, 0.365]);
  var c2 = ee.Image([-0.00071, -0.000164, -0.000329, -0.000479, -0.001012, -0.000966]);
  var c3 = ee.Image([0.000036, 0.000105, 0.00028, 0.005018, 0.004336, 0.004296]);
  var c4 = ee.Image([0.088, 0.0437, 0.0875, 0.1355, 0.056, 0.0155]);
  var c5 = ee.Image([0.0789, -1.2697, 0.1014, 0.6621, 0.7757, 0.639]);
  var cb = ee.Image([0.640, 0.31, 0.286, 0.189, 0.274, -0.186]);

  var tau_in = pair.multiply(c2).divide(cos_theta)
	.subtract(w.multiply(c3).add(c4).divide(cos_theta)).exp().multiply(c1).add(c5);

  var tau_out = pair.multiply(c2).subtract(w.multiply(c3).add(c4))
	.exp().multiply(c1).add(c5);

  var refl_sur_image = refl_toa_image
	.select(["B1", "B2", "B3", "B4", "B5", "B7"]).expression(
	  '(b() + cb * (tau_in - 1)) / (tau_in * tau_out)',
	  {'cb':cb, 'tau_in':tau_in, 'tau_out':tau_out}).clamp(0.0001, 1);

  return refl_sur_image;
};

// calculate at-surface reflectance for a collection
var refl_sur_coll = ee.ImageCollection(refl_toa_coll.map(tasumi_reflectance_func));

Map.addLayer(refl_sur_coll.select(['B3', 'B2', 'B1']).median(),{min:"0, 0, 0", max:"0.4, 0.4, 0.4"});
Map.setCenter(-118, 38, 6);
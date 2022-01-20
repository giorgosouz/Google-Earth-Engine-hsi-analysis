//Collect data from the ImageCollection filtered by date,clouds,roi

var collection = imageCollection
                    .filterBounds(geometry5)
                    .filterDate('2019-01-01', '2019-12-31')
                    .filter(ee.Filter.lte('CLOUD_COVER', 20));



//Choose the least cloudy image

var image = ee.Image(collection.sort('CLOUD_COVER').first());
image=image.clip(geometry5);


//Show on map
Map.addLayer(image,imageVisParam,'432',false)
Map.addLayer(image,imageVisParam2,'543',false)

//Calculate NDVI
var NDVI = image.normalizedDifference(['B5', 'B4'])
Map.addLayer(NDVI,imageVisParam3,'NVDI',false)



//collection2 for 2018
var collection2 = imageCollection
                    .filterBounds(geometry5)
                    .filterDate('2018-01-01', '2018-12-31')
                    .filter(ee.Filter.lte('CLOUD_COVER', 20));




//func for collection's ndvi
var addNDVI = function(image) {
  var ndvi = image.normalizedDifference(['B5', 'B4']).rename('NDVI');
  return image.addBands(ndvi);
};


//add nvdi to collection
var collection = collection.map(addNDVI);
var collection2 = collection2.map(addNDVI);

//func for doy
var addDoy = function(img){
  var doy = img.date().getRelative('day', 'year');
  var doyBand = ee.Image.constant(doy).uint16().rename('DOY');
  doyBand = doyBand.updateMask(img.select('NDVI').mask()); 

  return img.addBands(doyBand);
};


//add doy to collection
var collection= collection.map(addDoy);
var collection2= collection2.map(addDoy);





//generate mosaic from nvdi
var mosaic = collection.qualityMosaic('NDVI');
mosaic=mosaic.clip(geometry5);
Map.addLayer(mosaic,imageVisParam4,'mos19ndvi',false);
Map.addLayer(mosaic,imageVisParam5,'mos19doy',false);
var mosaic2 = collection2.qualityMosaic('NDVI');
mosaic2=mosaic2.clip(geometry5);
Map.addLayer(mosaic2,imageVisParam4,'mos18ndvi',false);
Map.addLayer(mosaic2,imageVisParam5,'mos18doy',false);






//setting up vars for the function

var l8toa= imageCollection

//Function ready from the slideshow provided
function karankCharts(roi,i){

// This field contains UNIX time in milliseconds.
var timeField = 'system:time_start';

// Use this function to mask clouds in Landsat 8 imagery.
var maskClouds = function(image) {
  var quality = image.select('BQA');
  var cloud01 = quality.eq(61440);
  var cloud02 = quality.eq(53248);
  var cloud03 = quality.eq(28672);
  var mask = cloud01.or(cloud02).or(cloud03).not();
  return image.updateMask(mask);
};

// Use this function to add variables for NDVI, time and a constant
// to Landsat 8 imagery.
var addVariables = function(image) {
  // Compute time in fractional years since the epoch.
  var date = ee.Date(image.get(timeField));
  var years = date.difference(ee.Date('1970-01-01'), 'year');
  // Return the image with the added bands.
  return image
    // Add an NDVI band.
    .addBands(image.normalizedDifference(['B5', 'B4']).rename('NDVI')).float()
    // Add a time band.
    .addBands(ee.Image(years).rename('t').float())
    // Add a constant band.
    .addBands(ee.Image.constant(1));
};

// Remove clouds, add variables and filter to the area of interest.
var filteredLandsat = l8toa
  .filterBounds(roi)
  .map(maskClouds)
  .map(addVariables)
  ;
// Plot a time series of NDVI at a single location.
var l8Chart = ui.Chart.image.series(filteredLandsat.select('NDVI'), roi)
    .setChartType('ScatterChart')
    .setOptions({
      title: 'Landsat 8 NDVI time series at ROI'+''+i,
      trendlines: {0: {
        color: 'CC0000'
      }},
      lineWidth: 1,
      pointSize: 3,
    });
print(l8Chart);

// Linear trend ----------------------------------------------------------------
// List of the independent variable names
var independents = ee.List(['constant', 't']);

// Name of the dependent variable.
var dependent = ee.String('NDVI');

// Compute a linear trend.  This will have two bands: 'residuals' and 
// a 2x1 band called coefficients (columns are for dependent variables).
var trend = filteredLandsat.select(independents.add(dependent))
    .reduce(ee.Reducer.linearRegression(independents.length(), 1));
// Map.addLayer(trend, {}, 'trend array image');

// Flatten the coefficients into a 2-band image
var coefficients = trend.select('coefficients')
  .arrayProject([0])
  .arrayFlatten([independents]);

// Compute a de-trended series.
var detrended = filteredLandsat.map(function(image) {
  return image.select(dependent).subtract(
          image.select(independents).multiply(coefficients).reduce('sum'))
          .rename(dependent)
          .copyProperties(image, [timeField]);
});

// Plot the detrended results.
var detrendedChart = ui.Chart.image.series(detrended, roi, null, 30)
    .setOptions({
      title: 'Detrended Landsat time series at ROI'+''+i,
      lineWidth: 1,
      pointSize: 3,
    });
print(detrendedChart);

// Harmonic trend ----------------------------------------------------------------
// Use these independent variables in the harmonic regression.
var harmonicIndependents = ee.List(['constant', 't', 'cos', 'sin']);

// Add harmonic terms as new image bands.
var harmonicLandsat = filteredLandsat.map(function(image) {
  var timeRadians = image.select('t').multiply(2 * Math.PI);
  return image
    .addBands(timeRadians.cos().rename('cos'))
    .addBands(timeRadians.sin().rename('sin'));
});
  
// The output of the regression reduction is a 4x1 array image.
var harmonicTrend = harmonicLandsat
  .select(harmonicIndependents.add(dependent))
  .reduce(ee.Reducer.linearRegression(harmonicIndependents.length(), 1));

// Turn the array image into a multi-band image of coefficients.
var harmonicTrendCoefficients = harmonicTrend.select('coefficients')
  .arrayProject([0])
  .arrayFlatten([harmonicIndependents]);

// Compute fitted values.
var fittedHarmonic = harmonicLandsat.map(function(image) {
  return image.addBands(
    image.select(harmonicIndependents)
      .multiply(harmonicTrendCoefficients)
      .reduce('sum')
      .rename('fitted'));
});

// Plot the fitted model and the original data at the ROI.
print(ui.Chart.image.series(
  fittedHarmonic.select(['fitted','NDVI']), roi, ee.Reducer.mean(), 30)
    .setSeriesNames(['NDVI', 'fitted'])
    .setOptions({
      title: 'Harmonic model: original and fitted values for ROI'+''+i,
      lineWidth: 1,
      pointSize: 3,
}));


var phase = harmonicTrendCoefficients.select('cos').atan2(
            harmonicTrendCoefficients.select('sin'));
            
var amplitude = harmonicTrendCoefficients.select('cos').hypot(
                harmonicTrendCoefficients.select('sin'));

var rgb = phase.unitScale(-Math.PI, Math.PI).addBands(     // hue
          amplitude.multiply(2.5)).addBands(            // saturation
          ee.Image(1)).hsvToRgb();                    // value

//Map.addLayer(rgb, {}, 'phase (hue), amplitude (saturation)');


}


//printing

karankCharts(geometry,1);
karankCharts(geometry2,2);
karankCharts(geometry3,3);
karankCharts(geometry4,4);

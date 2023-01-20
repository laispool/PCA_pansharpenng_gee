// Script to pansharp a multispectral image through PCA
// Earth Engine - Code Editor

// Create geometry
var geometry = ee.Geometry.Polygon(
         [[[-48.437665341965214,-25.9882255598952],
            [-48.90439125853334,-25.9882255598952],
            [-48.90439125853334,-26.483573101319053],
            [-48.437665341965214,-26.483573101319053],
            [-48.437665341965214,-25.9882255598952]]]);

// Slect and map a true and false color image
// see how to choose an image: https://github.com/laispool/ImageCollection_table 
var image = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20220703T132239_20220703T132238_T22JGS') // create the variable that contains the image
    .clip(geometry); // clip from geometry
var vis_false = { 
  min: 100,
  max: 2000,
  bands: ['B11', 'B3', 'B2'],
};
var vis_rgb = {
  min: 100,
  max: 2000,
  bands: ['B4', 'B3', 'B2'],
};

// Display the image.
Map.centerObject(geometry, 10), //Babitonga Bay (SC-Brazil) location
Map.addLayer(image, vis_rgb, 'true color composition',false);
Map.addLayer(image, vis_false, 'false color composition',false);
// BEGIN PCA
print('Begin PCA', image);
// Select which bands to use for the PCA
var PCAbands = ['B2', 'B3', 'B4', 'B8', 'B11'];

// Convert the image to a 2D array for the later matrix computations
var arrayImage = image.select(PCAbands).toArray();
print("array", arrayImage)
// Calculate the covariance using the reduceRegion method
var covar = arrayImage.reduceRegion({
    reducer: ee.Reducer.covariance(),
    maxPixels: 1e9
});
print("covar",covar)
// Extract the covariance matrix and store it as an array
var covarArray = ee.Array(covar.get('array'));
print("covarArray",covarArray)
//Compute and extract the eigenvectors
var eigens = covarArray.eigen();
var eigenVectors = eigens.slice(1, 1); //extract the eigenvectors for the PCA

// Perform matrix multiplication
var principalComponents = ee.Image(eigenVectors)
    .matrixMultiply(arrayImage.toArray(1)); //Each multiplication results in a principal component

// Convert back to a multi-band image and display the first principal component (pc1):    
var pcImage = principalComponents
    // Throw out an unneeded dimension, [[]] -> [].
    .arrayProject([0])
    // Make the one band array image a multi-band image, [] -> image.
    .arrayFlatten([['pc1', 'pc2', 'pc3', 'pc4', 'pc5']]);

// Stretch this to the appropriate scale.
Map.addLayer(pcImage.select('pc1'), {min:-455, max:-420, palette:['cyan','steelblue']}, 'pc1',false);
Map.addLayer(pcImage.select('pc2'), {min:-300, max:-310, palette:['cyan','steelblue']}, 'pc2',false);
Map.addLayer(pcImage.select('pc3'), {min:-100, max:-120, palette:['cyan','steelblue']}, 'pc3',false);

//The min and max values will need to change if you map different bands or locations.
var visParamsPCA = {
    bands: ['pc1', 'pc2', 'pc3'],
    min: [-455.09, -2, -5],
    max: [-417.59, -1, -4]
};

Map.addLayer(pcImage, visParamsPCA, 'PC_multi',false);

/*
Each axis dynamically captures some aspect of the variation within 
the dataset. Thus, the mapped PCA may differ substantially based on 
where you have performed the PCA and which bands you are mapping.
*/
// END PCA
print(pcImage,'End of PCA');

// BEGIN Pnasharpening through PCA
print('Begin Pnasharpening');

var img2 = image.select(['B8','B11']); //Change 'B3' to 'B4' or 'B8' to analise the others correlations
var pixelVals = img2.reduceRegion({
  reducer: ee.Reducer.toList(),
  geometry: geometry,
  scale: 200,
  });
print('Pixel Lists:', pixelVals); 
// Convert NIR and SWIR value lists to an array to be plotted along the y-axis.
var y = ee.List(pixelVals.get('B11'));
var x = ee.List(pixelVals.get('B8')); //Change 'B3' to 'B4' or 'B8' to analise the others correlations

// Define the chart and print it to the console.
var chart = ui.Chart.array.values({array: y, axis: 0, xLabels: x}).setOptions({
  title: 'Relationship Among Spectral Bands - NIR', // Change the title to match the band analised
  colors: ['cf513e'],
  hAxis: {
    title: 'B8 reflectance (x1e3)', //Change 'B3' to 'B4' or 'B8' to analise the others correlations
    titleTextStyle: {italic: false, bold: true}
  },
  vAxis: {
    title: 'SWIR reflectance (x1e3)',
    titleTextStyle: {italic: false, bold: true}
  },
  pointSize: 4,
  dataOpacity: 0.4,
  legend: {position: 'right'},
  trendlines: {
        0: {
          type: 'linear',
          color: 'lightblue',
          lineWidth: 3,
          opacity: 0.7,
          showR2: true,
          visibleInLegend: true},},
});
print(chart);

// Import the geeSharp module
var geeSharp = require('users/aazuspan/geeSharp:geeSharp');

// Pan-sharpen!
print(sharpened);

// After analysing the charts, choose the band that showed the bigest r^2
function sharpened (image) {
  var sharp = geeSharp.sharpen(image.select(['B11']), image.select(['B8'])).rename ('sharp');
  return image.addBands([sharp])}

var image = sharpened(image);
//print(image)

var vis_sharp = {
  min: 100,
  max: 2000,
  bands: ['sharp', 'B3', 'B2'],
};  
Map.addLayer(image, vis_sharp, 'sharpened');
print('End of Pansharpening');

// Metrics
var Q = geeSharp.quality(image.select(['B11']), image.select(['sharp']), "Q", geometry, 10, 1e9);
print ("The Q value is: ", Q);
/*The universal image quality index, Q, measures distortion between a 
reference image and a modified version of that image as the product of 
changes in correlation, luminance, and contrast. Values range from -1 
to 1, with 1 representing the lowest possible distortion (identical images).
See Wang & Bovik, "A Universal Image Quality Index", 2002 for a detailed explanation.
*/
var CC = geeSharp.quality(image.select(['B11']), image.select(['sharp']), "CC", geometry, 10, 1e9) ;
print("The CC value is: ", CC);
/*Pearson's Correlation Coefficient, CC, measures correlation between a 
reference image and a modified version of that image. A value of 1 
represents perfect correlation (identical images). See Wang & Bovik, 
"A Universal Image Quality Index", 2002 for a detailed explanation.
*/
var CML = geeSharp.quality(image.select(['B11']), image.select(['sharp']), "CML", geometry, 10, 1e9);
print("The CML value is: ", CML);
/*Change in Mean Luminance, CML, measures the change in luminance between 
a reference image and a modified version of that image. A value of 1 
represents no change (identical images). See Wang & Bovik, "A Universal 
Image Quality Index", 2002 for a detailed explanation.
*/
var CMC = geeSharp.quality(image.select(['B11']), image.select(['sharp']), "CMC", geometry, 10, 1e9);
print("The CMC value is: ", CMC);
/*Change in Mean Contrast, CMC, measures change in contrast a reference 
image and a modified version of that image. A value of 1 represents no 
change (identical images). See Wang & Bovik, "A Universal Image Quality 
Index", 2002 for a detailed explanation.
*/
var ERGAS = geeSharp.quality(image.select(['B11']), image.select(['sharp']), "ERGAS", 10, 20, geometry, 10, 1e9); //10 = high resol; 20 = low resol
print("The ERGAS value is: ", ERGAS);
/*Dimensionless Global Relative Error of Synthesis (ERGAS) measures 
spectral distortion relative to the sharpening ratio. This compensates 
for the fact that larger increases in spatial resolution will typically 
result in greater spectral distortion. Values near to zero indicate low 
distortion. See Vaiopoulos 2011.
*/
var RASE = geeSharp.quality(image.select(['B11']), image.select(['sharp']), "RASE", geometry, 10, 1e9);
print("The RASE value is: ", RASE);
/*Relative average spectral error (RASE) is a measure of spectral 
similarity. Values close to 0, relative to image intensity, represent 
minimal distortion. See Vaiopoulos 2011.
*/
var MSE = geeSharp.quality(image.select(['B11']), image.select(['sharp']), "MSE", geometry, 10, 1e9);
print("The MSE value is: ", MSE);
/*Mean squared error (MSE) measures the difference between a reference 
image (such as an unsharpened multiband image) and an assessment image 
(such as a pan-sharpened multiband image). Generally, low MSE values 
represent less distortion between the reference and assessment image.
*/
var RMSE = geeSharp.quality(image.select(['B11']), image.select(['sharp']), "RMSE", geometry, 10, 1e9);
print("The RMSE value is: ", RMSE);
/*Root mean squared error (RMSE) measures the root of the difference 
between a reference image (such as an unsharpened multiband image) and 
an assessment image (such as a pan-sharpened multiband image). Generally,
low RMSE values represent less distortion between the reference and 
assessment image.
*/
var bias = geeSharp.quality(image.select(['B11']), image.select(['sharp']), "bias", geometry, 10, 1e9);
print("The bias value is: ", bias);
/*Bias measures changes in mean intensity value due to image modification.
Values close to 0 represent minimal distortion. See Vaiopoulos 2011.
*/
var DIV = geeSharp.quality(image.select(['B11']), image.select(['sharp']), "DIV", geometry, 10, 1e9);
print("The DIV value is: ", DIV);
/*Difference in Variance (DIV) measures changes in variance of intensity 
values due to image modification. Values close to 0 represent minimal 
distortion. See Vaiopoulos 2011.
*/


// INDICES TESTS
print('Begin tests');
// Defining the index funciton
function indices (image) {
  var ndwi = image.normalizedDifference(['B3', 'B8']).rename ('ndwi'); //Mc Feeters, 1996
  var mndwi = image.normalizedDifference(['B3', 'B11']).rename('mndwi'); // Xu, 2006 (B11 = SWI1 = Short wave infrared 1)
  var awei = image.expression('(4 * (G - S)) - (0.25 * N + 2.75 * S) ',{ //Feyisa etal, 2014
   B: image.select('B2'),
   G: image.select('B3'), 
   S: image.select('B11'), 
   N: image.select('B8'),
   }).rename('awei');
   var awei_sharp = image.expression('(4 * (G - sharp)) - (0.25 * N + 2.75 * sharp) ',{ //Feyisa etal, 2014
   B: image.select('B2'),
   G: image.select('B3'), 
   sharp: image.select('sharp'),
   N: image.select('B8'),
   }).rename('awei_sharp');
  return image.addBands([ndwi,mndwi,awei,awei_sharp])}

//Applying the function
var NWI = indices(image);  
print (NWI);

var NDWI = NWI.select('ndwi');
var mNDWI = NWI.select('mndwi');
var AWEI = NWI.select('awei');
var AWEI2 = NWI.select('awei_sharp');

var singleBandVis = {
              'min': -0.05,
              'max': 1,
              palette: ['steelblue','white']
              };
Map.addLayer(NDWI,singleBandVis,'NDWI')
Map.addLayer(mNDWI,singleBandVis,'mNDWI')
Map.addLayer(AWEI,singleBandVis,'AWEI')
Map.addLayer(AWEI2,singleBandVis,'AWEI SHARP')

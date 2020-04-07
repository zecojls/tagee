
///////////////////////////////////////////
// Functions for calculating parameters  //
///////////////////////////////////////////

exports.calculateParameters = function(dem, bbox){

// Defining kernels to retrieve 3x3 window elevations

// Weights for a 3x3 kernel
var w00 = [0, 0, 0];
var w11 = [1, 0, 0];
var w12 = [0, 1, 0];
var w13 = [0, 0, 1];

// Neighborhood indices
var n1 = [w11, w00, w00];
var n2 = [w12, w00, w00];
var n3 = [w13, w00, w00];
var n4 = [w00, w11, w00];
var n5 = [w00, w12, w00];
var n6 = [w00, w13, w00];
var n7 = [w00, w00, w11];
var n8 = [w00, w00, w12];
var n9 = [w00, w00, w13];

// Kernel for each neighborhood index
var kerneln1 = ee.Kernel.fixed(3, 3, n1, 1, 1, false);
var kerneln2 = ee.Kernel.fixed(3, 3, n2, 1, 1, false);
var kerneln3 = ee.Kernel.fixed(3, 3, n3, 1, 1, false);
var kerneln4 = ee.Kernel.fixed(3, 3, n4, 1, 1, false);
var kerneln5 = ee.Kernel.fixed(3, 3, n5, 1, 1, false);
var kerneln6 = ee.Kernel.fixed(3, 3, n6, 1, 1, false);
var kerneln7 = ee.Kernel.fixed(3, 3, n7, 1, 1, false);
var kerneln8 = ee.Kernel.fixed(3, 3, n8, 1, 1, false);
var kerneln9 = ee.Kernel.fixed(3, 3, n9, 1, 1, false);

// Function to compute single neighborhood values
var addNParameters = function(dem) {
  var N1 = dem.convolve(kerneln1).rename('N1');
  var N2 = dem.convolve(kerneln2).rename('N2');
  var N3 = dem.convolve(kerneln3).rename('N3');
  var N4 = dem.convolve(kerneln4).rename('N4');
  var N5 = dem.convolve(kerneln5).rename('N5');
  var N6 = dem.convolve(kerneln6).rename('N6');
  var N7 = dem.convolve(kerneln7).rename('N7');
  var N8 = dem.convolve(kerneln8).rename('N8');
  var N9 = dem.convolve(kerneln9).rename('N9');
  return dem.addBands([N1,N2,N3,N4,N5,N6,N7,N8,N9]);
};

// Elevation values for each neighborhood
var demZPar = addNParameters(dem);
var demZParameters = demZPar.rename(['Elevation','Z1','Z2','Z3','Z4','Z5','Z6','Z7','Z8','Z9']);
//print(demZParameters, 'DEM with Z parameters');

// Defining positions
var demPositions = dem.addBands(ee.Image.pixelLonLat());

// Longitude
var long = demPositions.select('longitude');
var longNPar = addNParameters(long);
var longNParameters = longNPar.rename(['longitude','longN1','longN2','longN3','longN4','longN5','longN6','longN7','longN8','longN9']);
//print(longNParameters, 'Longitude with neighborhood parameters');

// Latitude
var lat = demPositions.select('latitude');
var latNPar = addNParameters(lat);
var latNParameters = latNPar.rename(['latitude','latN1','latN2','latN3','latN4','latN5','latN6','latN7','latN8','latN9']);
//print(latNParameters, 'Latitude with neighborhood parameters');

// Function of haversine formula to retrieve distances between two neighborhood points on a spheroidal grid
var haversineFunction = function(demLat, demLong, latNeigh1, latNeigh2, longNeigh1, longNeigh2) {
  var φ1 = demLat.select(ee.String(latNeigh1)).divide(180).multiply(Math.PI); // to Radians
  var φ2 = demLat.select(ee.String(latNeigh2)).divide(180).multiply(Math.PI); // to Radians
  var λ1 = demLong.select(ee.String(longNeigh1)).divide(180).multiply(Math.PI); // to Radians
  var λ2 = demLong.select(ee.String(longNeigh2)).divide(180).multiply(Math.PI); // to Radians
  
  var Δφ = φ2.subtract(φ1); // (φ2 - φ1)
  var Δλ = λ2.subtract(λ1); // (λ2 - λ1)
  var p1 = Δφ.divide(2).sin().multiply(Δφ.divide(2).sin()); // sin(Δφ/2) * sin(Δφ/2)
  var p2 = φ1.cos().multiply(φ2.cos()); // cos(φ1) * cos(φ2)
  var p3 = Δλ.divide(2).sin().multiply(Δλ.divide(2).sin()); // sin(Δλ/2) * sin(Δλ/2)
  var j = p2.multiply(p3).add(p1); // j = p1 + p2 * p3
  var p4 = ee.Image(ee.Number(1)).clip(bbox); // p4 = image with constant 1
  var p5 = ee.Image(ee.Number(2)).clip(bbox); // p5 = image with constant 2
  var p6 = p4.subtract(j).sqrt(); // sqrt(1-j)
  var p7 = j.sqrt(); // sqrt(a)
  var p8 = p6.atan2(p7); // atan2(p6,p7)
  var k = p5.multiply(p8); // k = 2 * p8
  var R = ee.Image(ee.Number(6371000)).clip(bbox); // approximate radius of Earth
  var l = R.multiply(k); // l = R * k which is the distance between two points
  
  return l;
};

// Distance values
var lenghtOfE = haversineFunction(latNParameters, longNParameters, 'latN1', 'latN4', 'longN1', 'longN4').rename('e');
var lenghtOfD = haversineFunction(latNParameters, longNParameters, 'latN4', 'latN7', 'longN4', 'longN7').rename('d');
var lenghtOfC = haversineFunction(latNParameters, longNParameters, 'latN1', 'latN2', 'longN1', 'longN2').rename('c');
var lenghtOfB = haversineFunction(latNParameters, longNParameters, 'latN4', 'latN5', 'longN4', 'longN5').rename('b');
var lenghtOfA = haversineFunction(latNParameters, longNParameters, 'latN7', 'latN8', 'longN7', 'longN8').rename('a');

// Merging all the parameters in a single image
var demCalculations = demZParameters.addBands(lenghtOfA)
                                    .addBands(lenghtOfB)
                                    .addBands(lenghtOfC)
                                    .addBands(lenghtOfD)
                                    .addBands(lenghtOfE);
return demCalculations;
};

////////////////////////////////////////////////////
// Functions for calculating terrain derivatives  //
////////////////////////////////////////////////////

exports.calculateDerivatives = function(parameters, bbox){

// Functions for Derivatives and Terrain Attributes

var addPDerivative = function(parameters) {
  var a = parameters.select('a');
  var b = parameters.select('b');
  var c = parameters.select('c');
  var d = parameters.select('d');
  var e = parameters.select('e');
  var Z1 = parameters.select('Z1');
  var Z3 = parameters.select('Z3');
  var Z4 = parameters.select('Z4');
  var Z6 = parameters.select('Z6');
  var Z7 = parameters.select('Z7');
  var Z9 = parameters.select('Z9');
  var constant2 = ee.Image(ee.Number(2)).clip(bbox);
  
  var p1 = a.pow(2).multiply(c).multiply(d).multiply(d.add(e)).multiply(Z3.subtract(Z1));
  var p2 = b.multiply(a.pow(2).multiply(d.pow(2)).add(c.pow(2).multiply(e.pow(2)))).multiply(Z6.subtract(Z4));
  var p3 = a.multiply(c.pow(2)).multiply(e.multiply(d.add(e))).multiply(Z9.subtract(Z7));
  var p4 = p1.add(p2).add(p3);
  
  var p5 = constant2;
  var p6 = a.pow(2).multiply(c.pow(2).multiply(d.add(e).pow(2)));
  var p7 = b.pow(2).multiply(a.pow(2).multiply(d.pow(2)).add(c.pow(2).multiply(e.pow(2))));
  var p8 = p6.add(p7).multiply(p5);
  
  var p = p4.divide(p8).rename('PDerivative');
  
  return p;
};

var addQDerivative = function(parameters) {
  var a = parameters.select('a');
  var b = parameters.select('b');
  var c = parameters.select('c');
  var d = parameters.select('d');
  var e = parameters.select('e');
  var Z1 = parameters.select('Z1');
  var Z2 = parameters.select('Z2');
  var Z3 = parameters.select('Z3');
  var Z4 = parameters.select('Z4');
  var Z5 = parameters.select('Z5');
  var Z6 = parameters.select('Z6');
  var Z7 = parameters.select('Z7');
  var Z8 = parameters.select('Z8');
  var Z9 = parameters.select('Z9');
  var constant1 = ee.Image(ee.Number(1)).clip(bbox);
  var constant2 = ee.Image(ee.Number(2)).clip(bbox);
  var constant3 = ee.Image(ee.Number(3)).clip(bbox);
  
  var p1 = constant1.divide(constant2.multiply(d).multiply(e).multiply(d.add(e)).multiply(a.pow(4).add(b.pow(4)).add(c.pow(4))));
  
  var p2 = d.pow(2).multiply(a.pow(4).add(b.pow(4)).add(b.pow(2).multiply(c.pow(2))));
  var p3 = c.pow(2).multiply(e.pow(2)).multiply(a.pow(2).subtract(b.pow(2)));
  var p4 = p2.add(p3).multiply(Z1.add(Z3));
  
  var p5 = d.pow(2).multiply(a.pow(4).add(c.pow(4)).add(b.pow(2).multiply(c.pow(2))));
  var p6 = e.pow(2).multiply(a.pow(4).add(c.pow(4)).add(a.pow(2).multiply(b.pow(2))));
  var p7 = p5.subtract(p6).multiply(Z4.add(Z6));
  
  var p9 = e.pow(2).multiply(b.pow(4).add(c.pow(4)).add(a.pow(2).multiply(b.pow(2))));
  var p10 = a.pow(2).multiply(d.pow(2)).multiply(b.pow(2).subtract(c.pow(2)));
  var p11 = p9.subtract(p10).multiply(Z7.add(Z9));
  
  var p13 = d.pow(2).multiply(b.pow(4).multiply(Z2.subtract(constant3.multiply(Z5))).add(c.pow(4).multiply(constant3.multiply(Z2).subtract(Z5))).add(a.pow(4).subtract(constant2.multiply(b.pow(2)).multiply(c.pow(2))).multiply(Z2.subtract(Z5))));
  
  var p14 = e.pow(2).multiply(a.pow(4).multiply(Z5.subtract(constant3.multiply(Z8))).add(b.pow(4).multiply(constant3.multiply(Z5).subtract(Z8))).add(c.pow(4).subtract(constant2.multiply(a.pow(2)).multiply(b.pow(2))).multiply(Z5.subtract(Z8))));
  
  var p15 = constant2.multiply(a.pow(2).multiply(d.pow(2)).multiply(b.pow(2).subtract(c.pow(2))).multiply(Z8).add(c.pow(2).multiply(e.pow(2)).multiply(a.pow(2).subtract(b.pow(2))).multiply(Z2)));
  
  var q = p1.multiply(p4.subtract(p7).subtract(p11).add(p13).add(p14).subtract(p15)).rename('QDerivative');
  
  return q;
};

var addRDerivative = function(parameters) {
  var a = parameters.select('a');
  var b = parameters.select('b');
  var c = parameters.select('c');
  var Z1 = parameters.select('Z1');
  var Z2 = parameters.select('Z2');
  var Z3 = parameters.select('Z3');
  var Z4 = parameters.select('Z4');
  var Z5 = parameters.select('Z5');
  var Z6 = parameters.select('Z6');
  var Z7 = parameters.select('Z7');
  var Z8 = parameters.select('Z8');
  var Z9 = parameters.select('Z9');
  var constant2 = ee.Image(ee.Number(2)).clip(bbox);
  
  var p1 = c.pow(2).multiply(Z1.add(Z3).subtract(constant2.multiply(Z2)));
  var p2 = b.pow(2).multiply(Z4.add(Z6).subtract(constant2.multiply(Z5)));
  var p3 = a.pow(2).multiply(Z7.add(Z9).subtract(constant2.multiply(Z8)));
  var p4 = a.pow(4).add(b.pow(4)).add(c.pow(4));
  
  var r = p1.add(p2).add(p3).divide(p4).rename('RDerivative');
  
  return r;
};

var addSDerivative = function(parameters) {
  var a = parameters.select('a');
  var b = parameters.select('b');
  var c = parameters.select('c');
  var d = parameters.select('d');
  var e = parameters.select('e');
  var Z1 = parameters.select('Z1');
  var Z2 = parameters.select('Z2');
  var Z3 = parameters.select('Z3');
  var Z4 = parameters.select('Z4');
  var Z5 = parameters.select('Z5');
  var Z6 = parameters.select('Z6');
  var Z7 = parameters.select('Z7');
  var Z8 = parameters.select('Z8');
  var Z9 = parameters.select('Z9');
  var constant2 = ee.Image(ee.Number(2)).clip(bbox);
  
  var p1 = c.multiply(a.pow(2).multiply(d.add(e)).add(b.pow(2).multiply(e))).multiply(Z3.subtract(Z1));
  var p2 = b.multiply(a.pow(2).multiply(d).subtract(c.pow(2).multiply(e))).multiply(Z4.subtract(Z6));
  var p3 = a.multiply(c.pow(2).multiply(d.add(e)).add(b.pow(2).multiply(d))).multiply(Z7.subtract(Z9));
  var p4 = p1.subtract(p2).add(p3);
  
  var p5 = constant2;
  var p6 = a.pow(2).multiply(c.pow(2).multiply(d.add(e).pow(2)));
  var p7 = b.pow(2).multiply(a.pow(2).multiply(d.pow(2)).add(c.pow(2).multiply(e.pow(2))));
  var p8 = p6.add(p7).multiply(p5);
  
  var s = p4.divide(p8).rename('SDerivative');
  
  return s;
};

var addTDerivative = function(parameters) {
  var a = parameters.select('a');
  var b = parameters.select('b');
  var c = parameters.select('c');
  var d = parameters.select('d');
  var e = parameters.select('e');
  var Z1 = parameters.select('Z1');
  var Z2 = parameters.select('Z2');
  var Z3 = parameters.select('Z3');
  var Z4 = parameters.select('Z4');
  var Z5 = parameters.select('Z5');
  var Z6 = parameters.select('Z6');
  var Z7 = parameters.select('Z7');
  var Z8 = parameters.select('Z8');
  var Z9 = parameters.select('Z9');
  var constant1 = ee.Image(ee.Number(1)).clip(bbox);
  var constant2 = ee.Image(ee.Number(2)).clip(bbox);
  var constant3 = ee.Image(ee.Number(3)).clip(bbox);
  
  var p1 = constant2.divide(constant3.multiply(d).multiply(e).multiply(d.add(e)).multiply(a.pow(4).add(b.pow(4)).add(c.pow(4))));
  
  var p2 = d.multiply(a.pow(4).add(b.pow(4)).add(b.pow(2).multiply(c.pow(2))));
  var p3 = c.pow(2).multiply(e).multiply(a.pow(2).subtract(b.pow(2)));
  var p4 = p2.subtract(p3).multiply(Z1.add(Z3));
  
  var p5 = d.multiply(a.pow(4).add(c.pow(4)).add(b.pow(2).multiply(c.pow(2))));
  var p6 = e.multiply(a.pow(4).add(c.pow(4)).add(a.pow(2).multiply(b.pow(2))));
  var p7 = p5.add(p6).multiply(Z4.add(Z6));
  
  var p9 = e.multiply(b.pow(4).add(c.pow(4)).add(a.pow(2).multiply(b.pow(2))));
  var p10 = a.pow(2).multiply(d).multiply(b.pow(2).subtract(c.pow(2)));
  var p11 = p9.add(p10).multiply(Z7.add(Z9));
  
  var p13 = d.multiply(b.pow(4).multiply(Z2.subtract(constant3.multiply(Z5))).add(c.pow(4).multiply(constant3.multiply(Z2).subtract(Z5))).add(a.pow(4).subtract(constant2.multiply(b.pow(2)).multiply(c.pow(2))).multiply(Z2.subtract(Z5))));
  
  var p14 = e.multiply(a.pow(4).multiply(constant3.multiply(Z8).subtract(Z5)).add(b.pow(4).multiply(Z8.subtract(constant3.multiply(Z5)))).add(c.pow(4).subtract(constant2.multiply(a.pow(2)).multiply(b.pow(2))).multiply(Z8.subtract(Z5))));
  
  var p15 = constant2.multiply(a.pow(2).multiply(d).multiply(b.pow(2).subtract(c.pow(2))).multiply(Z8).subtract(c.pow(2).multiply(e).multiply(a.pow(2).subtract(b.pow(2))).multiply(Z2)));
  
  var t = p1.multiply(p4.subtract(p7).add(p11).add(p13).add(p14).subtract(p15)).rename('TDerivative');
  
  return t;
};

var signPFunction = function(pDerivative) {
  var signP = pDerivative.expression("(b('PDerivative') > 0) ? 1" + ": (b('PDerivative') == 0) ? 0" + ": -1").clip(bbox).rename("signP");
  
  return signP;
};

var signQFunction = function(qDerivative) {
  var signQ = qDerivative.expression("(b('QDerivative') > 0) ? 1" + ": (b('QDerivative') == 0) ? 0" + ": -1").clip(bbox).rename("signQ");
  
  return signQ;
};

// Calculating the derivatives

var pDerivative = addPDerivative(parameters);
var qDerivative = addQDerivative(parameters);
var rDerivative = addRDerivative(parameters);
var sDerivative = addSDerivative(parameters);
var tDerivative = addTDerivative(parameters);
var signP = signPFunction(pDerivative);
var signQ = signQFunction(qDerivative);

var demWithDerivatives = parameters.addBands(pDerivative)
                                   .addBands(qDerivative)
                                   .addBands(rDerivative)
                                   .addBands(sDerivative)
                                   .addBands(tDerivative)
                                   .addBands(signP)
                                   .addBands(signQ);

return demWithDerivatives;

};

///////////////////////////////////////////////////
// Functions for calculating terrain attributes  //
///////////////////////////////////////////////////

exports.calculateAttributes = function(derivatives, bbox){

var slopeFunction = function(derivatives) {
  var p = derivatives.select('PDerivative');
  var q = derivatives.select('QDerivative');
  
  var p2 = p.pow(2).rename('A');
  var q2 = q.pow(2).rename('A');
  var p2q2 = ee.ImageCollection([p2,q2]);
  var sumP2q2 = p2q2.sum();
  var sqrtSumP2q2 = sumP2q2.sqrt();
  var slope = sqrtSumP2q2.atan().multiply(180).divide(Math.PI).rename('Slope');
  
  return slope;
};

var aspectFunction = function(derivatives) {
  var p = derivatives.select('PDerivative');
  var q = derivatives.select('QDerivative');
  var signP = derivatives.select('signP');
  var signQ = derivatives.select('signQ');
  var constant1 = ee.Image(ee.Number(1)).clip(bbox);
  var constantNeg1 = ee.Image(ee.Number(-1)).clip(bbox);
  var constant90 = ee.Image(ee.Number(90)).clip(bbox);
  var constant180 = ee.Image(ee.Number(180)).clip(bbox);
  
  var p1 = constantNeg1.multiply(constant90).multiply(constant1.subtract(signQ)).multiply(constant1.subtract(signP.abs()));
  var p2 = constant180.multiply(constant1.add(signP));
  var p3 = constant180.divide(Math.PI).multiply(signP);
  var p4 = constantNeg1.multiply(q).divide(p.pow(2).add(q.pow(2)).sqrt()).acos();
  var A = p1.add(p2).subtract(p3.multiply(p4)).rename('Aspect');
  
  return A;
};

var hillshadeFunction = function(aspect) {
  var p = derivatives.select('PDerivative');
  var q = derivatives.select('QDerivative');
  var θ = ee.Image(ee.Number(45)).clip(bbox); // Azimuth
  var ψ = ee.Image(ee.Number(315)).clip(bbox); // Elevation angle
  var constant1 = ee.Image(ee.Number(1)).clip(bbox);
  
  var p1 = constant1.subtract(p.multiply(θ.sin()).multiply(constant1.divide(ψ))).subtract(q.multiply(θ.cos()).multiply(constant1.divide(ψ)));
  var p2 = constant1.add(p.pow(2)).add(q.pow(2)).sqrt().multiply(constant1.add(θ.sin().multiply(constant1.divide(ψ)).pow(2)).add(θ.cos().multiply(constant1.divide(ψ)).pow(2)).sqrt());
  var AH = p1.divide(p2).rename('Hillshade');
  
  return AH;
};


var northernnessFunction = function(aspect) {
  var A = aspect.select('Aspect');
  
  var AN = A.multiply(Math.PI).divide(180).cos().rename('Northness');
  
  return AN;
};

var easternnessFunction = function(aspect) {
  var A = aspect.select('Aspect');
  
  var AE = A.multiply(Math.PI).divide(180).sin().rename('Eastness');
  
  return AE;
};

var horizontalCurvatureFunction = function(derivatives) {
  var p = derivatives.select('PDerivative');
  var q = derivatives.select('QDerivative');
  var r = derivatives.select('RDerivative');
  var s = derivatives.select('SDerivative');
  var t = derivatives.select('TDerivative');
  var constantNeg1 = ee.Image(ee.Number(-1)).clip(bbox);
  var constant1 = ee.Image(ee.Number(1)).clip(bbox);
  var constant2 = ee.Image(ee.Number(2)).clip(bbox);
  
  var p1 = q.pow(2).multiply(r).subtract(constant2.multiply(p).multiply(q).multiply(s)).add(p.pow(2).multiply(t));
  var p2 = p.pow(2).add(q.pow(2)).multiply(constant1.add(p.pow(2)).add(q.pow(2)).sqrt());
  var kh = constantNeg1.multiply(p1.divide(p2)).rename('HorizontalCurvature');
  
  return kh;
};

var verticalCurvatureFunction = function(derivatives) {
  var p = derivatives.select('PDerivative');
  var q = derivatives.select('QDerivative');
  var r = derivatives.select('RDerivative');
  var s = derivatives.select('SDerivative');
  var t = derivatives.select('TDerivative');
  var constantNeg1 = ee.Image(ee.Number(-1)).clip(bbox);
  var constant1 = ee.Image(ee.Number(1)).clip(bbox);
  var constant2 = ee.Image(ee.Number(2)).clip(bbox);

  var p2 = p.pow(2).multiply(r).add(constant2.multiply(p).multiply(q).multiply(s)).add(q.pow(2).multiply(t));
  var p3 = p.pow(2).add(q.pow(2)).multiply(constant1.add(p.pow(2)).add(q.pow(2)).pow(3).sqrt());
  var kv = constantNeg1.multiply(p2.divide(p3)).rename('VerticalCurvature');
  
  return kv;
};

var meanCurvatureFunction = function(derivatives) {
  var p = derivatives.select('PDerivative');
  var q = derivatives.select('QDerivative');
  var r = derivatives.select('RDerivative');
  var s = derivatives.select('SDerivative');
  var t = derivatives.select('TDerivative');
  var constantNeg1 = ee.Image(ee.Number(-1)).clip(bbox);
  var constant1 = ee.Image(ee.Number(1)).clip(bbox);
  var constant2 = ee.Image(ee.Number(2)).clip(bbox);

  var p2 = constant1.add(q.pow(2)).multiply(r).subtract(constant2.multiply(p).multiply(q).multiply(s)).add(constant1.add(p.pow(2)).multiply(t));
  var p3 = constant2.multiply(constant1.add(p.pow(2)).add(q.pow(2)).pow(3).sqrt());
  var km = constantNeg1.multiply(p2.divide(p3)).rename('MeanCurvature');
  
  return km;
};

var gaussianCurvatureFunction = function(derivatives) {
  var p = derivatives.select('PDerivative');
  var q = derivatives.select('QDerivative');
  var r = derivatives.select('RDerivative');
  var s = derivatives.select('SDerivative');
  var t = derivatives.select('TDerivative');
  var constant1 = ee.Image(ee.Number(1)).clip(bbox);
  
  var p1 = r.multiply(t).subtract(s.pow(2));
  var p2 = constant1.add(p.pow(2)).add(p.pow(2)).pow(2);
  var kg = p1.divide(p2).rename('GaussianCurvature');
  
  return kg;
};

var minimalCurvatureFunction = function(gaussian, mean) {
  var K = gaussian.select('GaussianCurvature');
  var H = mean.select('MeanCurvature');
  var kmin = H.subtract(H.pow(2).subtract(K).sqrt()).rename('MinimalCurvature');
  
  return kmin;
};

var maximalCurvatureFunction = function(gaussian, mean) {
  var K = gaussian.select('GaussianCurvature');
  var H = mean.select('MeanCurvature');
  var kmax = H.add(H.pow(2).subtract(K).sqrt()).rename('MaximalCurvature');
  return kmax;
};

var shapeIndexFunction = function(gaussian, mean) {
  var K = gaussian.select('GaussianCurvature').rename('K');
  var H = mean.select('MeanCurvature').rename('H');
  var constant2 = ee.Image(ee.Number(2));
  
  var index = constant2.divide(Math.PI).multiply(H.divide(H.pow(2).subtract(K).sqrt())).rename('ShapeIndex');
  
  return index;
};

// Calculating the Attributes

var slope = slopeFunction(derivatives);
var aspect = aspectFunction(derivatives);
var hillshade = hillshadeFunction(derivatives);
var northernness = northernnessFunction(aspect);
var easternness = easternnessFunction(aspect);
var horizontalCurvature = horizontalCurvatureFunction(derivatives);
var verticalCurvature = verticalCurvatureFunction(derivatives);
var meanCurvature = meanCurvatureFunction(derivatives);
var gaussianCurvature = gaussianCurvatureFunction(derivatives);
var minimalCurvature = minimalCurvatureFunction(gaussianCurvature, meanCurvature);
var maximalCurvature = maximalCurvatureFunction(gaussianCurvature, meanCurvature);
var shapeIndex = shapeIndexFunction(gaussianCurvature, meanCurvature);

var demWithAttributes = derivatives.addBands(slope)
                                   .addBands(aspect)
                                   .addBands(hillshade)
                                   .addBands(northernness)
                                   .addBands(easternness)
                                   .addBands(horizontalCurvature)
                                   .addBands(verticalCurvature)
                                   .addBands(meanCurvature)
                                   .addBands(gaussianCurvature)
                                   .addBands(minimalCurvature)
                                   .addBands(maximalCurvature)
                                   .addBands(shapeIndex);

return demWithAttributes.select('Elevation', 'Slope', 'Aspect', 'Hillshade', 'Northness', 'Eastness',
                                'HorizontalCurvature', 'VerticalCurvature', 'MeanCurvature',
                                'GaussianCurvature', 'MinimalCurvature', 'MaximalCurvature', 'ShapeIndex');
};

exports.terrainAnalysis = function(TAGEE, dem, bbox) {
  var parameters = TAGEE.calculateParameters(dem, bbox);
  var derivatives = TAGEE.calculateDerivatives(parameters, bbox);
  var attributes = TAGEE.calculateAttributes(derivatives, bbox);
  return(attributes);
};

/////////////////////////
// Additional features //
/////////////////////////

exports.makeVisualization = function(result, bandName, zoomLevel, bbox, palette){
  
  var levelsDic = ee.Dictionary({
        'level0': {'zoom': 0, 'scale': 157000},
        'level1': {'zoom': 1, 'scale': 78000},
        'level2': {'zoom': 2, 'scale': 39000},
        'level3': {'zoom': 3, 'scale': 20000},
        'level4': {'zoom': 4, 'scale': 10000},
        'level5': {'zoom': 5, 'scale': 5000},
        'level6': {'zoom': 6, 'scale': 2000},
        'level7': {'zoom': 7, 'scale': 1000},
        'level8': {'zoom': 8, 'scale': 611},
        'level9': {'zoom': 9, 'scale': 306},
        'level10': {'zoom': 10, 'scale': 153},
        'level11': {'zoom': 11, 'scale': 76},
        'level12': {'zoom': 12, 'scale': 38},
        'level13': {'zoom': 13, 'scale': 19},
        'level14': {'zoom': 14, 'scale': 10},
        'level15': {'zoom': 15, 'scale': 5},
  });

  var levelSelected = ee.Dictionary(levelsDic.get(zoomLevel));

  var imageSelected = result.select(bandName).rename('selection');

  var minMaxLegend = imageSelected.reduceRegion({
          reducer: ee.Reducer.percentile({percentiles: [5,95], outputNames: ['perc5','perc95']}),
          geometry: bbox,
          scale: levelSelected.get('scale'),
          bestEffort: true});

  var palettes = ee.Dictionary({
    'rainbow': '6e40aa, be3caf, fe4b83, ff7747, e3b62f, b0ef5a, 53f666, 1edfa2, 23acd8, 4c6fdc',
    'inferno': '000004, 160b39, 420a68, 6a176e, 932667, ba3655, dd513a, f3761b, fca50a, f6d746',
    'cubehelix': '163d4e, 1f6642, 53792f, a07949, d07e93, d09cd9, c1caf3',
    'red2green': 'a50026, d3322b, f16d43, fcab63, fedc8c, f9f7ae, d7ee8e, a4d86f, 64bc61, 23964f',
    'green2red': '23964f, 64bc61, a4d86f, d7ee8e, f9f7ae, fedc8c, fcab63, f16d43, d3322b, a50026',
    'elevation': 'b0f3be, e0fbb2, b8de76, 27a52a, 34883c, 9ca429, f8b004, c04a02, c04a02, 870800, 741805, 6c2a0a, 7d4a2b, 9c8170, b5b5b5, dad8da',
    'aspect': 'red, green, blue, yellow, red',    'hillshade': 'black, white'
  });
  
  var visualization = imageSelected.visualize({
    min: minMaxLegend.get('selection_perc5'),
    max: minMaxLegend.get('selection_perc95'),
    palette: palettes.get(palette)
  });
  
  return visualization;
  
};

exports.logTransformation = function(result, bandName){
  
  var selection = result.select(bandName).rename('selection');
  var sign = selection.expression("(b('selection') > 0) ? 1" + ": (b('selection') == 0) ? 0" + ": -1").rename("sign");
  var constant1 = ee.Image(ee.Number(1));
  var constant10 = ee.Image(ee.Number(10));
  var logValues = selection.abs().multiply(constant10.pow(4)).add(1).log10().multiply(sign).rename(bandName);
  
  return logValues;
  
};
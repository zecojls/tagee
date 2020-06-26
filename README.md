# Overview

Google Earth Engine (GEE) is a high-performance cloud-based platform for processing geospatial data.

There is a myriad of processing toolboxes in GEE, but this repository/package aims to contribute to the GEE community for making terrain analysis seamlessly. 

Terrain Analysis in Google Earth Engine (TAGEE) is a repository that contains the GEE Javascript code and a minimal reproducible example of TAGEE for the online code editor.

Available terrain attributes:

| Attribute            | Unit          | Description                                             |
|----------------------|---------------|---------------------------------------------------------|
| Elevation            | meter         | Height of terrain above sea level                       |
| Slope                | degree        | Slope gradient                                          |
| Aspect               | degree        | Compass direction                                       |
| Hillshade            | dimensionless | Brightness of the illuminated terrain                   |
| Northness            | dimensionless | Degree of orientation to North                          |
| Eastness             | dimensionless | Degree of orientation to East                           |
| Horizontal curvature | meter         | Curvature tangent to the contour line                   |
| Vertical curvature   | meter         | Curvature tangent to the slope line                     |
| Mean curvature       | meter         | Half-sum of the two orthogonal curvatures               |
| Minimal curvature    | meter         | Lowest value of curvature                               |
| Maximal curvature    | meter         | Highest value of curvature                              |
| Gaussian curvature   | meter         | Product of maximal and minimal curvatures               |
| Shape Index          | dimensionless | Continuous form of the Gaussian landform classification |

The users are referred to Florinsky (2016) for mathematical concepts of geomorphometry, a historical overview of the progress of digital terrain modelling, and the notion of the topographic surface and its limitations.

> Florinsky, Igor. Digital terrain analysis in soil science and geology. Academic Press, 2016.

Please, cite the following paper when using TAGEE:

> Safanelli, J.L.; Poppiel, R.R.; Ruiz, L.F.C.; Bonfatti, B.R.; Mello, F.A.O.; Rizzo, R.; Demattê, J.A.M. Terrain Analysis in Google Earth Engine: A Method Adapted for High-Performance Global-Scale Analysis. ISPRS Int. J. Geo-Inf. 2020, 9, 400. DOI: [https://doi.org/10.3390/ijgi9060400](https://doi.org/10.3390/ijgi9060400)

# Minimal reproducible example

[OPEN THE EXAMPLE DIRECTLY IN THE GEE CODE EDITOR](https://code.earthengine.google.com/2b1d977d8cb1c96dbf7c6a4d1064ea37).

NOTE: Any Earth Engine user with the above link can use it to view and run the example code. However, you need to be logged.

```javascript
// Importing module

var TAGEE = require('users/joselucassafanelli/TAGEE:TAGEE-functions');

// World bounding box

var bbox = ee.Geometry.Rectangle({coords: [-180, -60, 180, 60], geodesic: false});

// Water mask

var hansen_2016 = ee.Image('UMD/hansen/global_forest_change_2016_v1_4').select('datamask');
var hansen_2016_wbodies = hansen_2016.neq(1).eq(0);
var waterMask = hansen_2016.updateMask(hansen_2016_wbodies);

// Loading SRTM 30 m

var demSRTM = ee.Image('USGS/SRTMGL1_003').clip(bbox).rename('DEM');

// Smoothing filter
var gaussianFilter = ee.Kernel.gaussian({
  radius: 3, sigma: 2, units: 'pixels', normalize: true
});

// Smoothing the DEM with the gaussian kernel.
var demSRTM = demSRTM.convolve(gaussianFilter).resample("bilinear");

// Terrain analysis

var DEMAttributes = TAGEE.terrainAnalysis(TAGEE, demSRTM, bbox).updateMask(waterMask);
print(DEMAttributes.bandNames(), 'Parameters of Terrain');

// Visualization

var vizVC = TAGEE.makeVisualization(DEMAttributes, 'VerticalCurvature', 'level2', bbox, 'rainbow');
Map.addLayer(vizVC, {}, 'VerticalCurvature');
Map.setCenter(0,0,2);
```
## Example description

In the Javascript code editor [https://code.earthengine.google.com/](https://code.earthengine.google.com/), it is necessary to import the module **TAGEE-functions**.

```javascript
var TAGEE = require('users/joselucassafanelli/TAGEE:TAGEE-functions');
```

Then, you need to define the bounding box (study region limits) and import the Digital Elevation Model for terrain analysis (e.g. SRTM 30 m).

```javascript
// World bounding box

var bbox = ee.Geometry.Rectangle({coords: [-180, -60, 180, 60], geodesic: false});

// Water mask

var hansen_2016 = ee.Image('UMD/hansen/global_forest_change_2016_v1_4').select('datamask');
var hansen_2016_wbodies = hansen_2016.neq(1).eq(0);
var waterMask = hansen_2016.updateMask(hansen_2016_wbodies);

// Loading SRTM 30 m

var demSRTM = ee.Image('USGS/SRTMGL1_003').clip(bbox).rename('DEM');

// Smoothing filter
var gaussianFilter = ee.Kernel.gaussian({
  radius: 3, sigma: 2, units: 'pixels', normalize: true
});

// Smoothing the DEM with the gaussian kernel.
var demSRTM = demSRTM.convolve(gaussianFilter).resample("bilinear");
```

Finally, you can run the TAGEE module for calculating the terrain attributes.

```javascript
// Terrain analysis

var DEMAttributes = TAGEE.terrainAnalysis(TAGEE, demSRTM, bbox).updateMask(waterMask);
print(DEMAttributes.bandNames(), 'Parameters of Terrain');

List (13 elements){
0: Elevation
1: Slope
2: Aspect
3: Hillshade
4: Northness
5: Eastness
6: HorizontalCurvature
7: VerticalCurvature
8: MeanCurvature
9: GaussianCurvature
10: MinimalCurvature
11: MaximalCurvature
12: ShapeIndex}
```

TAGEE has an additional feature for making the visualization easier and adapated to a dynamic scale (legend), once the pixel distances for the derivatives calculations are influenced by the visualization level. The legend limits are estimated by the 5th and 95th percentiles existing within the bounding box.

```javascript
// Visualization

var vizVC = TAGEE.makeVisualization(DEMAttributes, 'VerticalCurvature', 'level2', bbox, 'rainbow');
Map.addLayer(vizVC, {}, 'VerticalCurvature');
Map.setCenter(0,0,2);
```

For the function `makeVisualization`, you need to specify:
- Attributes in multiband object
- Attribute name (string, e.g. 'VerticalCurvature')
- Visualization level (string, e.g. 'level2')
- Bounding box object
- Color pallete (string, e.g. 'rainbow')

The visualization levels are:

| Visualization level | Pixel resolution (m) |
|---------------------|----------------------|
| 'level0'            | 157000               |
| 'level1'            | 78000                |
| 'level2'            | 39000                |
| 'level3'            | 20000                |
| 'level4'            | 10000                |
| 'level5'            | 5000                 |
| 'level6'            | 2000                 |
| 'level7'            | 1000                 |
| 'level8'            | 611                  |
| 'level9'            | 306                  |
| 'level10'           | 153                  |
| 'level11'           | 76                   |
| 'level12'           | 38                   |
| 'level13'           | 19                   |
| 'level14'           | 10                   |
| 'level15'           | 5                    |

Available color palettes: 'rainbow', 'inferno', 'cubehelix', 'red2green', 'green2red', 'elevation', 'aspect' and 'hillshade'.


For any request, comment and suggestion, please send me a email (*my github nickname* at gmail.com)


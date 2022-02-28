__version__ = "0.1.0"

__all__ = ["terrainAnalysis", "makeVisualization", "logTransformation"]

# if the user hasn't done this already
import ee
ee.Initialize()

# Importing these functions makes them directly available to the user as a
# top-level import (e.g. from tagee import terrainAnalysis). Other functions
# can be imported as tagee.TAGEE._____
from .TAGEE import terrainAnalysis, makeVisualization, logTransformation
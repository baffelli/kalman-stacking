Python
------
- gdal
- snakemake
- pyparsing
- numpy

Geodata
-------
For the geocoding part, following files are needed in `.data/geo`:
- `DEM.tif` where `DEM` is the filename of the geotiff dem of the area of interest, as set in [config.json](./src/config.json)
- `swissTLM3D-2016-tlm_gelaendename_Clip`, a shapefile containing names of the area. It mus contain a feature named `feature_name` where the name is
set in [config.json](./src/config.json)
- `tablename.shp` where `tablename` is again set in  [config.json](./src/config.json). A shapefile with the location of a reference point.

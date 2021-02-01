
using ArchGDAL, GDAL

Lon = -44 .+ rand(5);
Lat = -23 .+ rand(5);

ArchGDAL.create(
        "C:/Russell/Projects/Geometry/Data/test/points.shp",
        driver = ArchGDAL.getdriver("ESRI Shapefile")
    ) do ds
    ArchGDAL.createlayer(
            geom = GDAL.wkbPoint,
            spatialref = ArchGDAL.importEPSG(4326);
        ) do layer
        for (lon,lat) in zip(Lon, Lat)
            ArchGDAL.createfeature(layer) do f
                ArchGDAL.setgeom!(f, ArchGDAL.createpoint(lon, lat));
            end
        end
        ArchGDAL.copy(layer, dataset = ds);
    end
end;

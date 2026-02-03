# dem-downloader

Download and reproject DEM tiles to a target CRS and resolution, with progress bars and automatic memory management.

Built on [dem-stitcher](https://github.com/ACCESS-Cloud-Based-InSAR/dem-stitcher) for tile fetching and [rasterio](https://rasterio.readthedocs.io/) for reprojection.

## Features

- **Auto source selection**: Switches from `glo_30` (30m) to `glo_90` (90m) when target resolution >= 100m, avoiding unnecessary memory usage
- **Size estimation**: Prints expected tile count and raw array size before downloading
- **Progress bars**: tqdm progress for reprojection and GeoTIFF writing
- **Nodata trimming**: Removes all-nodata border rows/columns after reprojection
- **CRS presets**: Built-in projections for Kazakhstan and Central Asia
- **Flexible CRS input**: EPSG codes, proj4 strings, WKT, or preset names

## Installation

```bash
pip install dem-stitcher rasterio tqdm
```

## Usage

### Command line

```bash
# Download with EPSG code
python dem_downloader.py --bbox 76.5 42.0 77.5 43.0 --epsg 32643 -r 500 -o dem.tif

# Download with CRS preset
python dem_downloader.py --bbox 46 40 88 56 --crs kaz_albers -r 500 -o kaz_dem.tif

# Download with custom proj4
python dem_downloader.py --bbox 68 50 78 55 \
  --crs "+proj=aea +lat_1=45 +lat_2=51 +lat_0=48 +lon_0=68 +ellps=WGS84 +units=m" \
  -o dem.tif

# Keep in WGS84 at native resolution
python dem_downloader.py --bbox 76.5 42.0 77.5 43.0 -o dem_wgs84.tif

# List available DEM sources
python dem_downloader.py --list-sources

# List CRS presets
python dem_downloader.py --list-presets
```

### As a library

```python
from dem_downloader import download_dem

download_dem(
    bbox=[76.5, 42.0, 77.5, 43.0],
    output_path="dem.tif",
    epsg=32643,
    resolution=500,
)
```

## DEM sources

| Source | Resolution | Coverage |
|--------|-----------|----------|
| `glo_30` | 30m | Global (default) |
| `glo_90` | 90m | Global |
| `srtm_v3` | 30m | 60°N to 56°S |
| `nasadem` | 30m | 60°N to 56°S |
| `3dep` | 10m | US only |

## CRS presets

| Name | Description |
|------|-------------|
| `kaz_albers` | Albers Equal Area for Kazakhstan (45-55°N, 46-87°E) |
| `kaz_north_albers` | Albers for Northern Kazakhstan (50-55°N, 60-78°E) |
| `kaz_lcc` | Lambert Conformal Conic for Kazakhstan |
| `central_asia_albers` | Albers for Central Asia (Tajikistan, Kyrgyzstan, etc.) |

## Memory management

When downloading large areas at coarse resolution, `glo_30` tiles can consume tens of GB of RAM just to be downsampled. The tool automatically switches to `glo_90` when target resolution is >= 100m, and warns when estimated memory exceeds 8 GB.

## License

MIT

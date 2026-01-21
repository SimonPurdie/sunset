# BREADCRUMBS

## Operational Notes

### Non-Earth Planetary Geometry Resolution (2026-01-21)
- **Issue**: Skyfield's `Topos` class is Earth-specific and cannot be used for other planets.
- **Solution**: Implemented direct geometric approach in `_find_planetary_sunset_terminator()` that:
  1. Randomly samples points on the planet's sphere using spherical coordinates
  2. Computes solar elevation by finding the dot product between surface normal and sun direction
  3. Returns points where sun is at sunset (elevation ∈ [-1.5°, +0.5°])
  4. Computes ENU sun direction without requiring planet rotation data
- **Impact**: Non-Earth planetary geometry (Mars, Venus, Mercury, Moon) now works without requiring PlanetaryConstants or IAU rotation data files.
- **Caveats**:
  - Latitude/longitude values are computed in ICRS-aligned frame rather than actual body-fixed coordinates
  - This is mathematically consistent but doesn't match actual lat/lon according to IAU definitions
  - Titan not supported (requires de440.bsp or custom ephemeris kernel)
- **Location**: `src/sunset/geometry/resolver.py` in `_find_planetary_sunset_terminator()` function.

### Skyfield API Limitations
- **Date**: 2026-01-20
- **Issue**: Skyfield's `.apparent()` method can fail with `EphemerisRangeError` due to relativistic calculation issues with Jupiter (body 599) and Saturn (body 699) for certain geometries.
- **Workaround**: Implemented fallback to manual RA/Dec to Alt/Az calculation using sidereal time when `.apparent()` fails.
- **Location**: `src/sunset/geometry/resolver.py` in `_get_alt_az()` function.

### Titan Ephemeris Limitations (2026-01-21)
- **Issue**: Titan is not included in the de421.bsp ephemeris kernel.
- **Status**: Titan test skipped - requires de440.bsp or custom ephemeris kernel.
- **Location**: `tests/test_integration.py::test_titan_full_pipeline`

### Topos Earth-Specific
- **Date**: 2026-01-20
- **Issue**: Skyfield's `Topos` class uses the IERS2010 geoid which is specific to Earth. Cannot add `Topos` to other planet vectors (Mars, Venus, Mercury, Moon).
- **Workaround**: Currently only Earth is supported for geometry resolution. Other planets require a different approach using PlanetaryConstants or custom latlon-to-ICRS transformations.
- **Impact**: Tests for Mars, Venus, and Mercury sunsets are skipped until this is resolved.
- **Location**: Tests in `tests/test_geometry.py` are marked as skipped for non-Earth bodies.

### Non-Earth Planetary Surface Positions (2026-01-20)
- **Issue**: Skyfield provides two methods for planetary surface positions:
  1. `Topos` class - Earth-specific, won't work for other planets
  2. `PlanetaryConstants` - Requires additional binary (.bcp) and text (.tf, .tpc) rotation data files
- **Investigation Findings**:
  - Attempted to use `PlanetaryConstants` approach with IAU2000 frames for Mars/Venus/Mercury
  - Binary rotation data files (e.g., `mars_pa_iau2000_1900-2050.bcp`) are not pre-downloaded
  - `load()` function doesn't know how to download these binary files automatically
  - Tried sub-solar point approach but couldn't accurately compute terminator without proper rotation matrices
- **Required Data Files** (examples):
  - Text files: `mars_latest_high_prec.tf`, `pck00010.tpc`
  - Binary files: `mars_pa_iau2000_1900-2050.bpc` (Mars), similar for Venus/Mercury
- **Solution Path**: Either:
  1. Download and manually load required IAU rotation data files for each planet
  2. Implement custom rotation matrix calculation using known planetary rotation parameters (W0, rotation rate) from IAU Working Group reports
  3. Use simplified geometry that approximates terminator location without exact rotation matrices
- **Status**: Implementation partially complete but tests skipped pending resolution.

### Future Work
- Download or compute rotation matrices for Mars (IAU2000), Venus (IAU2000), Mercury (IAU2000), Moon (DE421), and Titan (IAU2008)
- Alternative: Implement planet-specific topographic coordinate systems using published rotation rates and prime meridians from IAU WGCCRE reports
- Consider exploring `skyfield.planetarylib` or implementing custom latlon-to-ICRS transformations using rotation parameters.

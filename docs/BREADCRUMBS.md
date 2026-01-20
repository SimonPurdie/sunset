# BREADCRUMBS

## Operational Notes

### Skyfield API Limitations
- **Date**: 2026-01-20
- **Issue**: Skyfield's `.apparent()` method can fail with `EphemerisRangeError` due to relativistic calculation issues with Jupiter (body 599) and Saturn (body 699) for certain geometries.
- **Workaround**: Implemented fallback to manual RA/Dec to Alt/Az calculation using sidereal time when `.apparent()` fails.
- **Location**: `src/sunset/geometry/resolver.py` in `_get_alt_az()` function.

### Topos Earth-Specific
- **Date**: 2026-01-20
- **Issue**: Skyfield's `Topos` class uses the IERS2010 geoid which is specific to Earth. Cannot add `Topos` to other planet vectors (Mars, Venus, Mercury, Moon).
- **Workaround**: Currently only Earth is supported for geometry resolution. Other planets will need a different approach - likely using skyfield's `wgs84` latlon or manual ICRS coordinate transformations.
- **Impact**: Tests for Mars, Venus, and Mercury sunsets are skipped until this is resolved.
- **Location**: Tests in `tests/test_geometry.py` are marked as skipped for non-Earth bodies.

### Future Work
- Implement planet-specific topographic coordinate systems for Mars, Venus, Mercury, Moon, and Titan.
- Consider using `skyframe` module in skyfield or implementing custom latlon-to-ICRS transformations.

# It's Always Sunset Somewhere

A CPU-only renderer for physically plausible sunset scenes across the Solar System. The system uses real celestial mechanics and atmospheric optics to produce deterministic PNG images with embedded metadata.

## Status

Work in progress.

## Functionality

Given a UTC time and Solar System body, the system calculates observer coordinates where sunset occurs, models atmospheric composition and scattering properties, computes wavelength-dependent light transport, and renders a final image.

Aside from the project specification, all code has been authored by GLM-4.7.

## License

MIT - see [LICENSE](LICENSE)

# Implementation Plan: It's Always Sunset Somewhere

> **Core Objective**: Produce a deterministic system that renders physically plausible sunset scenes occurring somewhere in the Solar System at a given UTC time, using real celestial mechanics and physically grounded atmospheric optics — **without GPU acceleration** *(SPEC §1)*

---

## Checklist

- [x] 1. Project Setup & Dependencies
  - [x] 1.1 Update `pyproject.toml`
  - [x] 1.2 Create project package structure
  - [x] 1.3 Set up test infrastructure
- [ ] 2. Define Core Data Structures & Interfaces
  - [ ] 2.1 Create `models/bodies.py`
  - [ ] 2.2 Create `models/observer.py`
  - [ ] 2.3 Create `models/atmosphere.py`
  - [ ] 2.4 Create `models/spectral.py`
  - [ ] 2.5 Create `models/scene.py`
- [ ] 3. Celestial Geometry Resolver
  - [ ] 3.1 Implement `geometry/resolver.py`
  - [ ] 3.2 Use ephemeris data for sun position
  - [ ] 3.3 Implement terminator calculation
  - [ ] 3.4 Compute sun direction vector
  - [ ] 3.5 Calculate solar angular diameter
  - [ ] 3.6 Atomic tests
- [ ] 4. Atmospheric Profile Provider
  - [ ] 4.1 Create body-specific profiles
  - [ ] 4.2 Implement `atmosphere/provider.py`
  - [ ] 4.3 Implement Rayleigh scattering coefficients
  - [ ] 4.4 Implement barometric formula
  - [ ] 4.5 Store refractive index data
  - [ ] 4.6 Atomic tests
- [ ] 5. Visibility Elevation Resolver
  - [ ] 5.1 Implement `visibility/resolver.py`
  - [ ] 5.2 Implement optical depth calculation
  - [ ] 5.3 Implement altitude iteration
  - [ ] 5.4 Implement body rejection logic
  - [ ] 5.5 Atomic tests
- [ ] 6. Optical Path Integrator
  - [ ] 6.1 Implement `optics/integrator.py`
  - [ ] 6.2 Implement Rayleigh scattering physics
  - [ ] 6.3 Implement optical path integration
  - [ ] 6.4 Implement Mie scattering
  - [ ] 6.5 Compute skylight contribution
  - [ ] 6.6 Atomic tests
- [ ] 7. Spectral → Color Pipeline
  - [ ] 7.1 Implement spectral to XYZ conversion
  - [ ] 7.2 Implement exposure and tone mapping
  - [ ] 7.3 Implement XYZ to sRGB encoding
  - [ ] 7.4 Implement 16-bit output support
  - [ ] 7.5 Atomic tests
- [ ] 8. Software Renderer
  - [ ] 8.1 Implement CPU-only rendering engine
  - [ ] 8.2 Implement horizon geometry
  - [ ] 8.3 Implement sun rendering
  - [ ] 8.4 Implement sky rendering
  - [ ] 8.5 Implement ground/horizon rendering
  - [ ] 8.6 Ensure resolution ≥ 1024 × 512
  - [ ] 8.7 Atomic tests
- [ ] 9. Caption Generator
  - [ ] 9.1 Implement `caption/generator.py`
  - [ ] 9.2 Add required metadata content
  - [ ] 9.3 Ensure no speculative claims
  - [ ] 9.4 Atomic tests
- [ ] 10. PNG Metadata Embedding
  - [ ] 10.1 Implement `metadata/embedder.py`
  - [ ] 10.2 Add required metadata keys
  - [ ] 10.3 Implement 16-bit PNG writing
  - [ ] 10.4 Atomic tests
- [ ] 11. CLI Entry Point & Integration
  - [ ] 11.1 Implement CLI in `main.py`
  - [ ] 11.2 Orchestrate full pipeline
  - [ ] 11.3 Implement body iteration
  - [ ] 11.4 Error handling and validation
- [ ] 12. Determinism & Reproducibility Verification
  - [ ] 12.1 Implement reproducibility tests
  - [ ] 12.2 Ensure seed-controlled randomness
  - [ ] 12.3 Document numeric tolerance
  - [ ] 12.4 Atomic tests
- [ ] 13. End-to-End Testing & Validation
  - [ ] 13.1 Create integration test suite
  - [ ] 13.2 Contract validation
  - [ ] 13.3 Visual validation
  - [ ] 13.4 Cross-body validation
  - [ ] 13.5 Performance benchmarking

---

## Implementation Plan

### 1. Project Setup & Dependencies
*Source: SPEC §1, §6, pyproject.toml*

- **1.1** Update `pyproject.toml` with required dependencies:
  - Numerical computing: `numpy`
  - PNG handling with metadata: `Pillow` (or `pypng` for 16-bit)
  - Celestial mechanics: `astropy` or `skyfield` (ephemeris calculations)
  - Color science: `colour-science` (spectral to sRGB conversion, CIE functions)
- **1.2** Create project package structure under `src/sunset/`:
  ```
  src/sunset/
  ├── __init__.py
  ├── main.py              # CLI entry point
  ├── geometry/            # Celestial Geometry Resolver
  ├── atmosphere/          # Atmospheric Profile Provider
  ├── visibility/          # Visibility Elevation Resolver
  ├── optics/              # Optical Path Integrator
  ├── color/               # Spectral → Color Pipeline
  ├── renderer/            # Software Renderer
  ├── caption/             # Caption Generator
  ├── metadata/            # PNG Metadata Embedding
  └── models/              # Shared data models & types
  ```
- **1.3** Set up test infrastructure in `tests/` with pytest

---

### 2. Define Core Data Structures & Interfaces
*Source: SPEC §5 (all subsections), §2*

- **2.1** Create `models/bodies.py` — Solar System body definitions:
  - Body name, radius, atmosphere presence flag
  - Reference surface altitude definition
- **2.2** Create `models/observer.py` — Observer location dataclass:
  - Latitude, longitude, altitude, body reference
  - Solar elevation angle, sun direction vector
  - Solar angular diameter
  *(SPEC §5.1 outputs)*
- **2.3** Create `models/atmosphere.py` — Atmospheric profile dataclass:
  - Gas composition dict (fractions summing to 1.0)
  - Rayleigh scattering coefficient (wavelength-dependent function)
  - Optional Mie scattering parameters
  - Density vs altitude profile
  - Pressure vs altitude profile
  - Refractive index
  - Absorption bands
  *(SPEC §5.2 outputs)*
- **2.4** Create `models/spectral.py` — Spectral radiance representation:
  - Wavelength array and corresponding radiance values
  - Skylight contribution spectrum
  *(SPEC §5.4 outputs)*
- **2.5** Create `models/scene.py` — Complete scene metadata:
  - All fields required for PNG metadata embedding
  - `body`, `utc_time`, `latitude`, `longitude`, `altitude_m`, `random_seed`, `renderer_id`
  *(SPEC §4)*

---

### 3. Celestial Geometry Resolver
*Source: SPEC §5.1*

- **3.1** Implement `geometry/resolver.py`:
  - **Input**: UTC timestamp, body identifier, optional random seed
  - **Output**: Observer location with solar elevation ∈ [-1.5°, +0.5°]
- **3.2** Use ephemeris data (via `astropy`/`skyfield`) to compute:
  - Sun position relative to each body at given UTC
  - Body rotation to determine latitude/longitude where sun is at horizon
- **3.3** Implement terminator calculation:
  - Find the band of lat/lon coordinates where solar elevation is within [-1.5°, +0.5°]
  - Use random seed to select a specific point on this band
- **3.4** Compute sun direction vector in body-fixed coordinates
- **3.5** Calculate solar angular diameter based on body-sun distance
- **3.6** Tests *(SPEC §5.1 Atomic Tests)*:
  - `test_solar_elevation_bounds` — Returned solar elevation ∈ [-1.5°, +0.5°]
  - `test_coordinates_valid` — Lat/lon valid for body
  - `test_determinism` — Same seed → same result

---

### 4. Atmospheric Profile Provider
*Source: SPEC §5.2*

- **4.1** Create `atmosphere/profiles/` directory with body-specific profiles:
  - `earth.py` — N₂/O₂ atmosphere, standard density profile
  - `mars.py` — CO₂ atmosphere, thin
  - `venus.py` — Dense CO₂ atmosphere with cloud layers
  - `titan.py` — N₂/CH₄ atmosphere, thick haze
  - `airless.py` — Null atmosphere (Moon, Mercury, asteroids)
- **4.2** Implement `atmosphere/provider.py`:
  - **Input**: Body identifier, altitude
  - **Output**: Full atmospheric profile dataclass
- **4.3** Implement Rayleigh scattering coefficient calculation:
  - λ⁻⁴ wavelength dependence for each gas species *(SPEC §5.4)*
  - King correction factors for molecular anisotropy
- **4.4** Implement barometric formula for density/pressure vs altitude
- **4.5** Store refractive index data (for atmospheric refraction effects)
- **4.6** Tests *(SPEC §5.2 Atomic Tests)*:
  - `test_density_monotonic` — Density decreases with altitude
  - `test_pressure_monotonic` — Pressure decreases with altitude
  - `test_earth_matches_published` — Earth profile within published ranges
  - `test_mars_matches_published` — Mars profile within published ranges
  - `test_venus_matches_published` — Venus profile within published ranges
  - `test_titan_matches_published` — Titan profile within published ranges
  - `test_airless_null` — Airless bodies return null atmosphere

---

### 5. Visibility Elevation Resolver
*Source: SPEC §5.3*

- **5.1** Implement `visibility/resolver.py`:
  - **Input**: Atmospheric profile, observer location
  - **Output**: Final observer altitude, justification string
- **5.2** Implement optical depth calculation toward sun:
  - Integrate atmospheric density along line-of-sight to sun
  - Calculate total optical depth for visibility threshold check
- **5.3** Implement altitude iteration:
  - If optical depth > threshold, increase altitude
  - Binary search or iterative refinement for efficiency
- **5.4** Implement body rejection logic:
  - If no altitude satisfies visibility, reject body for this run
- **5.5** Tests *(SPEC §5.3 Atomic Tests)*:
  - `test_venus_surface_rejected` — Venus surface always rejected
  - `test_venus_upper_atmo_accepted` — Venus upper atmosphere accepted
  - `test_earth_surface_ok` — Earth never requires altitude change
  - `test_airless_surface_ok` — Airless bodies always accept surface

---

### 6. Optical Path Integrator
*Source: SPEC §5.4, §7*

- **6.1** Implement `optics/integrator.py`:
  - **Input**: Atmospheric profile, sun direction, observer altitude
  - **Output**: Spectral radiance, skylight contribution spectrum
- **6.2** Implement Rayleigh scattering physics *(SPEC §7 — NOT hardcoded)*:
  - Scattering cross-section σ(λ) ∝ λ⁻⁴
  - Phase function for Rayleigh scattering (1 + cos²θ)
  - Compute in-scattering contributions
- **6.3** Implement optical path integration:
  - Numerical integration along observer's viewing ray
  - Accumulate attenuation and scattering contributions
  - Multiple wavelength channels (e.g., 380nm–780nm in steps)
- **6.4** Implement Mie scattering (when aerosols present):
  - Particle size distribution effects
  - Henyey-Greenstein phase function approximation
- **6.5** Compute skylight contribution:
  - Single-scattering approximation (or multi-scattering if needed)
  - Hemisphere integration for sky brightness
- **6.6** Tests *(SPEC §5.4 Atomic Tests)*:
  - `test_rayleigh_wavelength_scaling` — Scattering ∝ 1/λ⁴
  - `test_short_wavelength_skylight` — Blue dominates skylight at high elevation
  - `test_optical_depth_increases_horizon` — Optical depth increases as sun → horizon
  - `test_short_wavelength_attenuation` — Blue attenuates faster in thick atmospheres
  - `test_airless_zero_skylight` — Airless bodies produce zero skylight
  - `test_no_negative_radiance` — All radiance values ≥ 0

---

### 7. Spectral → Color Pipeline
*Source: SPEC §5.5, §3*

- **7.1** Implement `color/spectral.py` — spectral to XYZ conversion:
  - CIE 1931 2° standard observer color matching functions
  - Integration of spectral radiance × CMF
- **7.2** Implement `color/tonemapping.py` — exposure and tone mapping:
  - Automatic exposure based on scene luminance
  - Tone mapping operator (Reinhard or ACES) preserving hue
- **7.3** Implement `color/encoding.py` — XYZ to sRGB:
  - XYZ → linear sRGB matrix transform
  - Gamut mapping for out-of-gamut colors
  - sRGB gamma encoding (linear → sRGB transfer function)
- **7.4** Implement 16-bit output support *(SPEC §3)*:
  - Linear light internal representation
  - 16-bit per channel PNG output (preferred)
  - 8-bit dithering fallback if needed
- **7.5** Tests *(SPEC §5.5 Atomic Tests)*:
  - `test_energy_luminance` — More energy → higher luminance
  - `test_hue_ordering` — Wavelength dominance preserves hue
  - `test_zero_radiance_black` — Zero radiance → black
  - `test_tonemapping_hue_preservation` — Tone mapping preserves hue ordering

---

### 8. Software Renderer
*Source: SPEC §5.6, §3, §7*

- **8.1** Implement `renderer/engine.py` — CPU-only rendering engine:
  - **Constraint**: No GPU, OpenGL, Vulkan, compute shaders *(SPEC §5.6)*
  - Deterministic floating-point behavior
- **8.2** Implement horizon geometry:
  - Curved horizon based on body radius and observer altitude
  - Horizon line calculation (geometric horizon angle)
- **8.3** Implement sun rendering:
  - Solar disc with correct angular diameter from Geometry Resolver
  - Sun position intersecting/within ±1.5° of horizon *(SPEC §3)*
  - Limb darkening (optional enhancement)
- **8.4** Implement sky rendering:
  - For each pixel above horizon, compute viewing ray direction
  - Call Optical Path Integrator for spectral radiance
  - Convert to sRGB via Spectral → Color Pipeline
- **8.5** Implement ground/horizon rendering:
  - Simple dark ground below horizon (color derived from atmosphere)
  - No hardcoded colors *(SPEC §7)*
- **8.6** Ensure resolution ≥ 1024 × 512 *(SPEC §3)*
- **8.7** Tests *(SPEC §5.6 Atomic Tests)*:
  - `test_horizon_present` — Horizon visible and horizontal
  - `test_sun_intersects_horizon` — Sun within tolerance of horizon
  - `test_no_color_banding` — No visible banding at 16-bit depth

---

### 9. Caption Generator
*Source: SPEC §5.7*

- **9.1** Implement `caption/generator.py`:
  - **Input**: Scene metadata dataclass
  - **Output**: Human-readable caption string
- **9.2** Required content from metadata only:
  - Body name
  - Location (lat, lon) and altitude
  - UTC time
  - Atmospheric explanation (or explicit absence for airless bodies)
- **9.3** No speculative or fictional claims *(SPEC §5.7)*
- **9.4** Tests *(SPEC §5.7 Atomic Tests)*:
  - `test_altitude_mentioned` — Mentions altitude if not surface
  - `test_no_atmosphere_mentioned` — Mentions lack of atmosphere if applicable
  - `test_no_speculation` — No speculative/fictional claims

---

### 10. PNG Metadata Embedding
*Source: SPEC §4, §2*

- **10.1** Implement `metadata/embedder.py`:
  - Write PNG text chunks with required fields
  - Validate metadata consistency with image content
- **10.2** Required metadata keys:
  - `body`: Solar System body name
  - `utc_time`: ISO-8601 UTC timestamp
  - `latitude`: degrees
  - `longitude`: degrees
  - `altitude_m`: meters above reference surface
  - `random_seed`: integer or hash
  - `renderer_id`: version or commit hash
- **10.3** Implement 16-bit PNG writing with embedded metadata
- **10.4** Tests:
  - `test_metadata_present` — All required fields present in PNG
  - `test_metadata_parseable` — Metadata values are parseable
  - `test_metadata_consistent` — Metadata matches rendered scene

---

### 11. CLI Entry Point & Integration
*Source: SPEC §2, src/sunset/main.py*

- **11.1** Implement CLI in `main.py`:
  - Arguments: `--utc-time`, `--body` (optional), `--seed` (optional), `--output`
  - Default behavior: auto-select body with visible sunset
- **11.2** Orchestrate full pipeline:
  1. Parse inputs and initialize random seed
  2. Call Celestial Geometry Resolver to find sunset location
  3. Call Atmospheric Profile Provider
  4. Call Visibility Elevation Resolver
  5. Call Optical Path Integrator
  6. Call Software Renderer to produce image
  7. Call Caption Generator
  8. Call PNG Metadata Embedder
  9. Write output file
- **11.3** Implement body iteration:
  - If a body is rejected (no visible sunset), try next body
  - Seed-controlled ordering of body selection
- **11.4** Error handling and validation:
  - Validate output meets contract *(SPEC §2)*
  - Return non-zero exit code on failure

---

### 12. Determinism & Reproducibility Verification
*Source: SPEC §6*

- **12.1** Implement reproducibility tests:
  - Same inputs → identical PNG output
  - ByteCompare verification
- **12.2** Ensure all randomness is seed-controlled:
  - Audit all random number usage
  - Numpy random state management
- **12.3** Document numeric tolerance if any floating-point variance is allowed
- **12.4** Tests:
  - `test_deterministic_output` — Run twice with same inputs, compare bytes
  - `test_seed_affects_location` — Different seeds produce different locations
  - `test_seed_preserved_metadata` — Seed recorded in PNG metadata

---

### 13. End-to-End Testing & Validation
*Source: SPEC §8, §2, §7*

- **13.1** Create integration test suite:
  - Full pipeline execution for each supported body
  - Verify output PNG meets all contracts
- **13.2** Contract validation:
  - PNG format, lossless *(SPEC §3)*
  - Resolution ≥ 1024 × 512 *(SPEC §3)*
  - Sun within ±1.5° of horizon *(SPEC §3)*
  - All required metadata present *(SPEC §4)*
- **13.3** Visual validation (manual + automated):
  - No hardcoded color gradients *(SPEC §7)*
  - Rayleigh scattering not approximated *(SPEC §7)*
  - Physics-derived sky colors
- **13.4** Cross-body validation:
  - Earth sunsets match "golden hour" expectations
  - Mars sunsets have blue hues near sun (known phenomenon)
  - Airless bodies have black sky with sun disc only
- **13.5** Performance benchmarking:
  - Ensure reasonable render time on CPU
  - Memory usage within practical limits

---

## Dependencies Summary

| Package | Purpose |
|---------|---------|
| `numpy` | Numerical arrays, integration |
| `Pillow` or `pypng` | PNG I/O with 16-bit support |
| `astropy` or `skyfield` | Ephemeris, celestial mechanics |
| `colour-science` | Spectral→XYZ→sRGB, CIE CMFs |
| `pytest` | Testing framework |

---

## Notes

- All components are independently replaceable and testable *(SPEC §5)*
- No GPU acceleration — pure CPU rendering *(SPEC §5.6)*
- No hardcoded colors — all derived from physics *(SPEC §7)*
- Determinism is mandatory *(SPEC §6)*

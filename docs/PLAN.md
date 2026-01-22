# It's Always Sunset Somewhere - Implementation Roadmap

## Living Document

This document tracks implementation tasks resulting from testing, bug fixes, and enhancements beyond the original SPEC.md. Tasks are prioritized by criticality. Keep this document current as you perform tasks from it.

---

## Priority 0: Critical Issues (Blockers)

These issues prevent the system from meeting the SPEC.md contract or user requirements.

### ~~Sun Rendering Bug~~ (COMPLETED 2026-01-21)
- ~~Issue:** Sun appeared as a diffuse black circle in the top-right corner; should be centered horizontally and intersect the horizon
- **Root cause:** Renderer was ignoring `observer.sun_direction` (ENU coordinates) and computing its own hardcoded sun direction based only on `solar_elevation_deg`. Also, `u` and `v` coordinates in `direction_map` precomputation were computed wrong (using `/width` and `/height` instead of `/center_x` and `/center_y`).
- **Fix:**
  1. Changed `u = (X - center_x) / width` to `u = (X - center_x) / center_x`
  2. Changed `v = (center_y - Y) / height` to `v = (center_y - Y) / center_y`
  3. Changed renderer to use `observer.sun_direction` instead of hardcoded direction
  4. Fixed renderer coordinate system: forward=[0,1,0] (looking toward horizon along +Y), right=[1,0,0], up=[0,0,1]
  5. Fixed coordinate mapping from ENU to renderer: (East, North, Up) → (X, Y, Z)
  6. Updated `_compute_sky_brightness_factors` to take `sun_direction` parameter
- **Result:** Sun now correctly positioned relative to horizon based on actual observer.sun_direction
- **Tests:** All 169 tests pass

### ~~Output Filename Management~~ (COMPLETED 2026-01-21)
- **Change:** Save output images to `output/YYYYMMDD-HHMMSS-planetname.png` instead of overwriting `sunset.png`
- **Example:** `output/20260121-143052-mars.png`
- **Rationale:** Preserve previous renders for comparison; easier inspection and debugging
- **Implementation:**
  - Changed default `--output` argument to `None` in `parse_args()`
  - Modified `render_sunset()` to auto-generate output path when `output_path` is `None`
  - Generated path format: `output/{timestamp}-{bodyname}.png` where timestamp uses `YYYYMMDD-HHMMSS` format
  - Manual override still available via `--output` flag
  - Updated test `test_parse_args_defaults` to expect `None` instead of `"sunset.png"`
- **Tests:** All 169 tests pass

### ~~Sky Gradient Visibility~~ (COMPLETED 2026-01-22)
- **Issue:** Test exists and passes but rendered images showed incorrect gradient - zenith brighter than horizon during sunset
- **Root cause:** `_compute_sky_brightness_factors()` used angle from sun, which for a sunset (sun at horizon) created uniform brightness. The function didn't properly model the physical distribution of skylight.
- **Fix:** Rewrote `_compute_sky_brightness_factors()` to use vertical angle from horizon instead of angle from sun. The new implementation:
  * Uses cos(vertical_angle) to create gradient where sky is brightest at horizon and darkest at zenith
  * Creates brightness factor ranging from 0.3 to 1.0 (previously 0.4 to 1.0)
  * Properly models that during sunset, skylight is brightest near the horizon (closer to sun) and dimmer at zenith
- **Result:** Sky gradient now properly shows horizon brighter than zenith during sunset
- **Tests:** All 169 tests pass

---

## Priority 1: Important Improvements (Should Fix)

### CLI/UX Enhancements

### ~~Better Error Reporting~~ (COMPLETED 2026-01-22)
- Provide actionable guidance in error messages
- Include context (what was being attempted when error occurred)
- Suggest next steps or possible fixes
- **Implementation:** Created new `src/sunset/errors.py` module with:
  * Base `SunsetError` class that includes message and suggestions
  * Specialized error types: `BodyNotFoundError`, `SunsetNotFoundError`, `VisibilityError`, `TimeParseError`, `AtmosphereProfileError`
  * `handle_error()` and `print_error()` utility functions
- Updated all modules to use custom error types:
  * `geometry/resolver.py` - uses `BodyNotFoundError`, `SunsetNotFoundError`, `TimeParseError`
  * `atmosphere/provider.py` - uses `BodyNotFoundError`, `AtmosphereProfileError`
  * `main.py` - uses `handle_error()`, `VisibilityError`, `print_error()`
- Enhanced validation error messages with bug report suggestions
- **Result:** All error messages now include:
  * Clear description of what went wrong
  * Actionable suggestions for resolution
  * Context about what was being attempted
- **Tests:** All 169 tests pass
- **Example outputs:**
  ```
  Error: Unknown body: 'unknown_planet'

  Suggestions:
    - Available bodies: earth, mars, mercury, moon, titan, venus
    - Check spelling (body names are case-insensitive)
  ```

  ```
  Error while parsing UTC time:
  Error: Invalid UTC time format: 'invalid-time'

  Suggestions:
    - Use ISO-8601 format with 'Z' suffix for UTC (e.g., '2024-01-15T18:00:00Z')
    - Omit the command to use the current time
    - Example: --utc-time 2024-01-15T18:00:00Z
  ```

### ~~Informative Output Before Rendering~~ (COMPLETED 2026-01-22)
- Print body name being rendered before rendering starts
- Print relevant metadata (UTC time, location coordinates, altitude, solar elevation)
- Format should be clear and human-readable
- **Implementation:** Added print statements in `render_sunset()` before rendering that display:
  * Body name
  * UTC time
  * Location (latitude, longitude)
  * Altitude
  * Solar elevation
- **Tests:** All 169 tests pass
- **Example output:**
  ```
  Rendering sunset on Mars
    UTC time: 2024-01-15T18:00:00Z
    Location: 21.214289°, -173.318240°
    Altitude: 0.0 m
    Solar elevation: 0.022°
  ```

### ~~Progress Indication~~ (COMPLETED 2026-01-22)
- Display a terminal progress indicator during rendering
- Show completion status with spinner
- **Implementation:** Added `Spinner` class in `src/sunset/main.py` that:
  * Shows animated spinner using Unicode characters during rendering
  * Threaded implementation to avoid blocking render process
  * Displays "Rendering scene" message with spinner animation
  * Clears spinner on completion and prints "Rendering complete"
- **Tests:** All 169 tests pass
- **Note:** Light spinner implementation without external dependencies; render time is typically <1 second for standard 1024×512 resolution

**Better Error Reporting**
- Provide actionable guidance in error messages
- Include context (what was being attempted when error occurred)
- Suggest next steps or possible fixes

### ~~Suppress Colour Library Warnings~~ (COMPLETED 2026-01-22)
- Suppress warnings about missing SciPy and Matplotlib dependencies
- Do not install these libraries unnecessarily
- Use warning filters or configure colour library to suppress these specific warnings
- **Implementation:** Added `warnings.filterwarnings("ignore", module="colour.utilities.verbose")` in `src/sunset/__init__.py` to suppress ColourUsageWarning about missing SciPy and Matplotlib dependencies
- **Result:** Warning count reduced from 39 to 37 in pytest output; ColourUsageWarning no longer appears during test runs or normal operation
- **Tests:** All 169 tests pass

---

## Priority 2: Nice-to-Have Features

### ~~Sun Disc Improvements~~ (COMPLETED 2026-01-22)
- **Enhancement:** Added limb darkening and atmospheric glow/halo effects to sun disc
- **Implementation:**
  - Implemented `_compute_sun_with_limb_darkening_and_glow()` method using linear limb darkening law: I(μ) = I₀ * (1 - u * (1 - μ))
    * Limb darkening coefficient u = 0.6 (typical for visible light)
    * μ = cos(θ) where θ is angle between line of sight and surface normal
    * Sun appears darker at edges (limb) because we see cooler outer layers vs deeper hotter layers at center
  - Implemented `_compute_sun_glow()` method for atmospheric halo effect
    * Glow extends to 2.5× sun angular radius
    * Glow intensity fades from 30% of sun color at disc edge to 0 at glow boundary
    * Only applies to atmospheric bodies (airless bodies have no atmosphere for scattering)
  - Modified `render()` method to use limb-darkened sun and add glow effect
- **Result:** Sun disc now shows realistic limb darkening with soft atmospheric glow during sunset
- **Tests:** All 169 tests pass
- **Physics:** Limb darkening is physically grounded in solar photosphere geometry and atmospheric scattering

### ~~Camera/Viewing Adjustments~~ (COMPLETED 2026-01-22)
- **Change:** Camera now tilted upward by 15.75° by default, placing horizon at ~15% from bottom of image
- **Implementation:**
  - Added `camera_pitch_deg` parameter to `Renderer.__init__()` with default of 15.75°
  - Modified `_precompute_direction_vectors()` to apply rotation matrix to camera basis vectors based on pitch angle
  - Added `--camera-pitch` CLI argument for user customization
  - Updated `render_scene()` and `render_scene_to_pil()` to accept and pass through camera_pitch_deg
  - Formula: `camera_pitch = arcsin(horizon_percent_from_bottom * tan(fov_vertical/2))` where horizon_percent_from_bottom = 0.15
- **Result:** Horizon now appears at approximately 15% from bottom of image (was at 50%), showing more sky area
- **Tests:** All 169 tests pass
- **Example usage:**
  - Default: `sunset-render` (horizon at 15% from bottom)
  - Custom: `sunset-render --camera-pitch 0` (horizon at center, 50% from bottom)
  - More upward: `sunset-render --camera-pitch 30` (horizon at ~27% from bottom)

### ~~Stars Visible in Darker Sky Regions~~ (COMPLETED 2026-01-22)
- **Enhancement:** Added star rendering for dark sky regions, especially on airless bodies
- **Implementation:**
  - Created `src/sunset/renderer/stars.py` module with:
    * Bright star catalog (RA/Dec coordinates and magnitudes for ~20 bright stars)
    * Functions to convert RA/Dec to ICRS Cartesian coordinates
    * Apparent magnitude to brightness conversion using standard formula: B = 10^(-0.4 * magnitude)
  - Implemented `_render_stars()` method in `Renderer` class:
    * Converts star positions from ICRS to ENU coordinates using observer's longitude/latitude
    * Renders stars only where sky brightness is below visibility threshold
    * Airless bodies: stars visible where sky brightness <= 0.8 (brighter threshold)
    * Atmospheric bodies: stars visible only where sky brightness <= 0.2 (darker threshold)
    * Stars rendered as small points with varying brightness based on magnitude
  - Modified `render()` method to add stars to final image
- **Result:** Stars now appear in dark sky regions, particularly visible on airless bodies like Moon and Mercury
- **Tests:** All 169 tests pass (including `test_moon_black_sky` which confirms sky remains dark)
- **Physics:** Stars are rendered using standard photometric magnitude scale; ENU conversion uses simplified longitude-based rotation sufficient for visualization

### Visual Quality Enhancements
- Sun disc improvements: limb darkening, atmospheric glow/halo effects
- ~~Stars visible in darker sky regions~~
- Ground/terrain representation instead of flat black horizon

### Debugging/Diagnostics
- Intermediate visualization mode: spectral data, direction vectors, masks
- Verbose flag for detailed internal state output

---

## Testing Strategy for Sun Issue

**Approach:** Test-driven investigation rather than prescriptive test specification

**Goal:** Builders should write tests to discover and pinpoint the sun rendering problem

**Suggested areas for investigation (not prescriptive):**
- Where is the sun actually being positioned in the rendered image?
- What color values are being computed for sun pixels?
- Is the sun_mask correctly identifying which pixels should contain the sun?
- Are there coordinate system mismatches between sun direction vector and viewing direction map?

**Guidance:** Let builders determine appropriate tests based on what they uncover during investigation

# It’s Always Sunset Somewhere

## 1. Core Objective

Produce a deterministic system that renders **physically plausible sunset scenes** occurring *somewhere in the Solar System at a given UTC time*, using real celestial mechanics and physically grounded atmospheric optics. The system must be able to do this **without GPU acceleration**, producing a single canonical image output per run.

Success is judged entirely by **observable outputs**, not implementation choices.

---

## 2. Canonical Output Contract

Each successful run **must** produce:

1. **A PNG image** depicting a sunset in progress
2. **Embedded metadata** sufficient to reproduce and explain the scene
3. **Deterministic reproducibility** under identical inputs

Failure to satisfy any item constitutes a failed run.

---

## 3. Render Output Requirements

### Image Format

* Format: **PNG** (lossless)
* Resolution: **≥ 1024 × 512** pixels
* Orientation:

  * Horizon is horizontal
  * Sun is intersecting or within ±1.5° of the horizon

### Color & Precision

* Internal representation: **linear light**
* Output encoding: **sRGB**
* Bit depth:

  * **16-bit per channel preferred**
  * 8-bit allowed only if dithering is applied

### Visual Constraints

* Sun has a finite angular diameter (varies by body)
* Sky exhibits a vertical luminance and chromatic gradient
* No hardcoded color palettes or gradients

---

## 4. Required Embedded Metadata (PNG text chunks)

Each PNG must embed the following key–value data:

* `body`: Solar System body name
* `utc_time`: ISO-8601 UTC timestamp
* `latitude`: degrees
* `longitude`: degrees
* `altitude_m`: meters above reference surface
* `random_seed`: integer or hash
* `renderer_id`: version or commit hash

Metadata must be internally consistent with the image.

---

## 5. High-Level System Decomposition (Contract-Based)

Each component below is independently replaceable and testable. Components may not assume internal details of others.

---

### 5.1 Celestial Geometry Resolver

**Purpose**: Identify observer locations where a sunset is occurring at the given UTC time.

**Inputs**

* UTC timestamp
* Solar System body identifier
* Optional random seed

**Outputs**

* One observer location containing:

  * Latitude
  * Longitude
  * Solar elevation angle ∈ [-1.5°, +0.5°]
  * Sun direction vector in body-fixed coordinates
  * Solar angular diameter

**Atomic Tests**

* Returned solar elevation satisfies bounds
* Coordinates are valid for the body
* Deterministic under fixed seed

---

### 5.2 Atmospheric Profile Provider

**Purpose**: Describe atmospheric composition and structure, including wavelength-dependent scattering mechanisms.

**Purpose**: Describe atmospheric composition and structure.

**Inputs**

* Body identifier
* Altitude

**Outputs**

* Gas composition (fractions sum to 1.0)
* Rayleigh scattering coefficient as a function of wavelength
* Optional Mie scattering parameters (particle size, density) if aerosols or clouds are present
* Density vs altitude
* Pressure vs altitude
* Refractive index
* Absorption bands (if any)

**Atomic Tests**

* Density and pressure decrease monotonically with altitude
* Known bodies (Earth, Mars, Venus, Titan) match published ranges
* Airless bodies return a null atmosphere

---

### 5.3 Visibility Elevation Resolver

**Purpose**: Ensure the sun is directly visible.

**Rule**
If the sun is not visible at the surface:

* Increase altitude until optical depth < visibility threshold
* If no such altitude exists, reject this body for this run

**Outputs**

* Final observer altitude
* Justification string

**Atomic Tests**

* Venus surface rejected, upper atmosphere accepted
* Earth never requires altitude change
* Airless bodies always accept surface

---

### 5.4 Optical Path Integrator

**Purpose**: Compute wavelength-dependent attenuation and scattering, with Rayleigh scattering treated as a first-class mechanism.

**Purpose**: Compute wavelength-dependent attenuation and scattering.

**Inputs**

* Atmospheric profile
* Sun direction
* Observer altitude

**Outputs**

* Spectral radiance reaching observer
* Skylight contribution spectrum

**Atomic Tests (Atmosphere-Dependent)**

* Rayleigh scattering strength scales approximately with 1/λ⁴ for molecular atmospheres
* In Rayleigh-dominated atmospheres, shorter wavelengths contribute more strongly to skylight at high solar elevations
* Optical depth increases as solar elevation → 0°
* Short wavelengths attenuate faster than long wavelengths in thick atmospheres
* Airless bodies produce zero skylight
* No negative radiance values

---

### 5.5 Spectral → Color Pipeline

**Purpose**: Convert spectral radiance into displayable color.

**Inputs**

* Spectral radiance
* Exposure / white point

**Outputs**

* Linear tristimulus values
* sRGB-encoded color values

**Atomic Tests (Atmosphere-Independent)**

* Increased total energy increases luminance
* Relative wavelength dominance preserves hue ordering
* Zero radiance → black
* Tone mapping preserves hue ordering

---

### 5.6 Software Renderer

**Purpose**: Produce the final image.

**Constraints**

* **CPU-only** rendering
* Deterministic floating-point behavior
* No GPU, OpenGL, Vulkan, or compute shaders

**Responsibilities**

* Horizon geometry consistent with body radius
* Solar disc rendered with correct angular size
* Sky color derived from physics outputs only

**Outputs**

* Final image buffer

**Atomic Tests**

* Horizon present and horizontal
* Sun intersects horizon within tolerance
* No visible color banding at 16-bit depth

---

### 5.7 Caption Generator (Metadata-Derived)

**Purpose**: Explain the scene via metadata (no free-form invention).

**Required Content**

* Body name
* Location and altitude
* UTC time
* Atmospheric explanation or explicit absence

**Atomic Tests**

* Mentions altitude if not surface
* Mentions lack of atmosphere if applicable
* No speculative or fictional claims

---

## 6. Determinism & Reproducibility

* Identical inputs **must** produce identical PNG output (within numeric tolerance if explicitly defined)
* Randomness must be seed-controlled

---

## 7. Disallowed Behaviors

The system must not:

* Approximate Rayleigh scattering via hardcoded color gradients
* Collapse Rayleigh effects into post-hoc color correction

The system must not:

* Hardcode sky or sunset colors
* Invent atmospheric properties
* Ignore solar distance or angular size
* Fake sunset geometry
* Depend on GPU acceleration

---

## 8. Evaluation Criterion

The system is considered successful if:

> A domain expert in planetary science agrees that the rendered sunset is physically plausible for the stated body, location, and time.

---

## 9. Optional Extensions (Non-Blocking)

* PNG sequences for time progression
* Multiple candidate sunsets per run
* Side-by-side body comparisons

These must not weaken or bypass the core output contracts.
# It's Always Sunset Somewhere - Implementation Roadmap

## Living Document

This document tracks implementation tasks resulting from testing, bug fixes, and enhancements beyond the original SPEC.md. Tasks are prioritized by criticality. Keep this document current as you perform tasks from it.

---

## Priority 0: Critical Issues (Blockers)

These issues prevent the system from meeting the SPEC.md contract or user requirements.

### Output Filename Management
- **Change:** Save output images to `output/YYYYMMDD-HHMMSS-planetname.png` instead of overwriting `sunset.png`
- **Example:** `output/20260121-143052-mars.png`
- **Rationale:** Preserve previous renders for comparison; easier inspection and debugging

### Sun Rendering Bug
- **Issue:** Sun appears as a diffuse black circle in the top-right corner; should be centered horizontally and intersect the horizon
- **Root cause unknown; investigate via tests to determine:**
  - Actual sun position in image coordinates
  - Color values from optics pipeline for sun pixels
  - Sun_mask correctness
  - Coordinate system alignment between sun direction and viewing vectors

### Sky Gradient Visibility
- **Issue:** Sky appears as a single hue; SPEC.md requires "vertical luminance and chromatic gradient"
- **Investigate:** Whether sky brightness factors are producing visible variation across the image

---

## Priority 1: Important Improvements (Should Fix)

### CLI/UX Enhancements

**Informative Output Before Rendering**
- Print body name being rendered before rendering starts
- Print relevant metadata (UTC time, location coordinates, altitude, solar elevation)
- Format should be clear and human-readable

**Progress Indication**
- Display a terminal progress bar during rendering
- Show percentage complete or render stage information

**Better Error Reporting**
- Provide actionable guidance in error messages
- Include context (what was being attempted when error occurred)
- Suggest next steps or possible fixes

**Suppress Colour Library Warnings**
- Suppress warnings about missing SciPy and Matplotlib dependencies
- Do not install these libraries unnecessarily
- Use warning filters or configure colour library to suppress these specific warnings

---

## Priority 2: Nice-to-Have Features

### Camera/Viewing Adjustments
- Adjust camera angle so horizon appears at ~15% of image height (currently at 50%)
- Expose camera angle/viewing direction as CLI parameters

### Visual Quality Enhancements
- Sun disc improvements: limb darkening, atmospheric glow/halo effects
- Stars visible in darker sky regions
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

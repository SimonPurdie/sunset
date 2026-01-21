import argparse
import sys
from datetime import datetime
from pathlib import Path

from .geometry.resolver import resolve_sunset_location
from .atmosphere.provider import get_body_profile
from .visibility.resolver import resolve_visibility_elevation
from .renderer.engine import render_scene
from .caption.generator import generate_caption
from .metadata.embedder import embed_metadata
from .models import Observer, SOLAR_SYSTEM_BODIES, Scene
from . import __version__


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Render a physically plausible sunset scene for a Solar System body."
    )
    parser.add_argument(
        "--utc-time",
        type=str,
        default=datetime.utcnow().isoformat(timespec="seconds") + "Z",
        help="ISO-8601 UTC timestamp (default: current time)",
    )
    parser.add_argument(
        "--body",
        type=str,
        choices=list(SOLAR_SYSTEM_BODIES.keys()),
        help="Solar System body to render (default: auto-select)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for deterministic output",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="sunset.png",
        help="Output PNG file path (default: sunset.png)",
    )
    parser.add_argument(
        "--caption",
        action="store_true",
        help="Print caption to stdout",
    )
    return parser.parse_args()


def get_body_order(seed: int | None) -> list[str]:
    """Get ordered list of bodies to try based on seed.

    Args:
        seed: Random seed for deterministic ordering

    Returns:
        List of body IDs in order to attempt
    """
    import numpy as np

    bodies = list(SOLAR_SYSTEM_BODIES.keys())
    rng = np.random.default_rng(seed)
    return list(rng.permutation(bodies))


def create_scene(
    observer: Observer, utc_time: str, random_seed: int, renderer_id: str
) -> Scene:
    """Create a Scene object from an Observer.

    Args:
        observer: Observer location and geometry
        utc_time: UTC timestamp string
        random_seed: Random seed used
        renderer_id: Renderer version identifier

    Returns:
        Scene dataclass with all required metadata
    """
    return Scene(
        body=observer.body,
        utc_time=utc_time,
        latitude=observer.latitude,
        longitude=observer.longitude,
        altitude_m=observer.altitude_m,
        random_seed=random_seed,
        renderer_id=renderer_id,
    )


def validate_scene(scene: Scene, observer: Observer) -> bool:
    """Validate that the scene meets the output contract.

    Args:
        scene: Scene metadata
        observer: Observer geometry

    Returns:
        True if validation passes
    """
    if -1.5 > observer.solar_elevation_deg or observer.solar_elevation_deg > 0.5:
        print(
            f"Error: Solar elevation {observer.solar_elevation_deg:.3f}° outside "
            f"required range [-1.5°, +0.5°]",
            file=sys.stderr,
        )
        return False

    if not (-90 <= scene.latitude <= 90):
        print(
            f"Error: Latitude {scene.latitude}° outside valid range [-90, 90]",
            file=sys.stderr,
        )
        return False

    if not (-180 <= scene.longitude <= 180):
        print(
            f"Error: Longitude {scene.longitude}° outside valid range [-180, 180]",
            file=sys.stderr,
        )
        return False

    if scene.altitude_m < 0:
        print(
            f"Error: Altitude {scene.altitude_m}m cannot be negative",
            file=sys.stderr,
        )
        return False

    return True


def render_sunset(
    utc_time: str,
    body_id: str | None = None,
    random_seed: int | None = None,
    output_path: str = "sunset.png",
    print_caption: bool = False,
) -> int:
    """Render a sunset scene and save to PNG with embedded metadata.

    Args:
        utc_time: ISO-8601 UTC timestamp
        body_id: Body identifier (None for auto-selection)
        random_seed: Random seed for deterministic output
        output_path: Path to save output PNG
        print_caption: If True, print caption to stdout

    Returns:
        Exit code (0 for success, non-zero for failure)
    """
    try:
        import numpy as np
        import hashlib

        seed_value: int
        if random_seed is None:
            hash_digest = hashlib.sha256(utc_time.encode()).digest()
            seed_value = (
                int.from_bytes(hash_digest[:4], byteorder="little") & 0x7FFFFFFF
            )
        else:
            seed_value = random_seed

        renderer_id = f"sunset-{__version__}"

        bodies_to_try = [body_id] if body_id else get_body_order(seed_value)

        for body_to_try in bodies_to_try:
            if body_to_try is None:
                continue

            try:
                observer = resolve_sunset_location(utc_time, body_to_try, seed_value)
            except ValueError as e:
                if body_id:
                    print(f"Error: {e}", file=sys.stderr)
                    return 1
                continue

            atmospheric_profile = get_body_profile(
                observer.body.name, observer.altitude_m
            )

            final_altitude, justification, is_visible = resolve_visibility_elevation(
                atmospheric_profile, observer
            )

            if not is_visible:
                if body_id:
                    print(
                        f"Error: Sun not visible on {observer.body.name} at this time",
                        file=sys.stderr,
                    )
                    return 1
                continue

            observer = Observer(
                latitude=observer.latitude,
                longitude=observer.longitude,
                altitude_m=final_altitude,
                body=observer.body,
                solar_elevation_deg=observer.solar_elevation_deg,
                sun_direction=observer.sun_direction,
                solar_angular_diameter_deg=observer.solar_angular_diameter_deg,
            )

            scene = create_scene(observer, utc_time, seed_value, renderer_id)

            if not validate_scene(scene, observer):
                return 1

            image_array = render_scene(observer)

            output_file = Path(output_path)
            output_file.parent.mkdir(parents=True, exist_ok=True)

            embed_metadata(image_array, scene, str(output_file))

            if print_caption:
                caption = generate_caption(scene, atmospheric_profile)
                print(caption)

            return 0

        print(
            "Error: Could not find a body with a visible sunset at this time",
            file=sys.stderr,
        )
        return 1

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback

        traceback.print_exc()
        return 1


def main():
    """CLI entry point."""
    args = parse_args()

    exit_code = render_sunset(
        utc_time=args.utc_time,
        body_id=args.body,
        random_seed=args.seed,
        output_path=args.output,
        print_caption=args.caption,
    )

    sys.exit(exit_code)


if __name__ == "__main__":
    main()

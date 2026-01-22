import argparse
import sys
import threading
import time
from datetime import datetime
from pathlib import Path
from itertools import cycle

from .geometry.resolver import resolve_sunset_location
from .atmosphere.provider import get_body_profile
from .visibility.resolver import resolve_visibility_elevation
from .renderer.engine import render_scene
from .caption.generator import generate_caption
from .metadata.embedder import embed_metadata
from .models import Observer, SOLAR_SYSTEM_BODIES, Scene
from . import __version__
from .errors import (
    SunsetError,
    VisibilityError,
    handle_error,
    print_error,
)


class Spinner:
    """Simple terminal spinner for long-running operations."""

    def __init__(self, message: str):
        self.message = message
        self._stop_event = threading.Event()
        self._spinner_thread = None
        self._chars = cycle(["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"])

    def _spin(self):
        while not self._stop_event.is_set():
            char = next(self._chars)
            sys.stdout.write(f"\r{self.message} {char} ")
            sys.stdout.flush()
            time.sleep(0.1)

    def start(self):
        self._stop_event.clear()
        self._spinner_thread = threading.Thread(target=self._spin, daemon=True)
        self._spinner_thread.start()

    def stop(self):
        if self._spinner_thread:
            self._stop_event.set()
            self._spinner_thread.join()
            sys.stdout.write("\r" + " " * (len(self.message) + 3) + "\r")
            sys.stdout.flush()


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
        default=None,
        help="Output PNG file path (default: output/YYYYMMDD-HHMMSS-planetname.png)",
    )
    parser.add_argument(
        "--caption",
        action="store_true",
        help="Print caption to stdout",
    )
    parser.add_argument(
        "--camera-pitch",
        type=float,
        default=15.75,
        help="Camera pitch angle in degrees (positive = looking up, default 15.75° for horizon at 15%% from bottom)",
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
        error = SunsetError(
            f"Solar elevation {observer.solar_elevation_deg:.3f}° outside "
            f"required range [-1.5°, +0.5°]",
            suggestions=[
                "This indicates a bug in the geometry resolver",
                "The resolver should always return a valid sunset location",
                "Please report this issue with the command-line arguments used",
            ],
        )
        print_error(error)
        return False

    if not (-90 <= scene.latitude <= 90):
        error = SunsetError(
            f"Latitude {scene.latitude}° outside valid range [-90, 90]",
            suggestions=[
                "This indicates a bug in the geometry resolver",
                "The resolver should always return a valid latitude",
                "Please report this issue with the command-line arguments used",
            ],
        )
        print_error(error)
        return False

    if not (-180 <= scene.longitude <= 180):
        error = SunsetError(
            f"Longitude {scene.longitude}° outside valid range [-180, 180]",
            suggestions=[
                "This indicates a bug in the geometry resolver",
                "The resolver should always return a valid longitude",
                "Please report this issue with the command-line arguments used",
            ],
        )
        print_error(error)
        return False

    if scene.altitude_m < 0:
        error = SunsetError(
            f"Altitude {scene.altitude_m}m cannot be negative",
            suggestions=[
                "This indicates a bug in the visibility resolver",
                "Altitude should never be negative",
                "Please report this issue with the command-line arguments used",
            ],
        )
        print_error(error)
        return False

    return True


def render_sunset(
    utc_time: str,
    body_id: str | None = None,
    random_seed: int | None = None,
    output_path: str | None = None,
    print_caption: bool = False,
    camera_pitch_deg: float = 15.75,
) -> int:
    """Render a sunset scene and save to PNG with embedded metadata.

    Args:
        utc_time: ISO-8601 UTC timestamp
        body_id: Body identifier (None for auto-selection)
        random_seed: Random seed for deterministic output
        output_path: Path to save output PNG (None for auto-generated)
        print_caption: If True, print caption to stdout
        camera_pitch_deg: Camera pitch angle in degrees (positive = looking up, default 15.75° for horizon at 15% from bottom)

    Returns:
        Exit code (0 for success, non-zero for failure)
    """
    try:
        import numpy as np
        import hashlib
        from .errors import SunsetError, TimeParseError

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
            except TimeParseError as e:
                return handle_error(e, "parsing UTC time")
            except SunsetError as e:
                if body_id:
                    return handle_error(
                        e, f"resolving sunset location for {body_to_try}"
                    )
                continue

            atmospheric_profile = get_body_profile(
                observer.body.name, observer.altitude_m
            )

            final_altitude, justification, is_visible = resolve_visibility_elevation(
                atmospheric_profile, observer
            )

            if not is_visible:
                if body_id:
                    return handle_error(
                        VisibilityError(observer.body.name, justification),
                        f"checking visibility on {observer.body.name}",
                    )
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

            print(f"Rendering sunset on {observer.body.name}")
            print(f"  UTC time: {utc_time}")
            print(f"  Location: {observer.latitude:.6f}°, {observer.longitude:.6f}°")
            print(f"  Altitude: {final_altitude:.1f} m")
            print(f"  Solar elevation: {observer.solar_elevation_deg:.3f}°")
            print()

            spinner = Spinner("Rendering scene")
            spinner.start()
            image_array = render_scene(observer, camera_pitch_deg=camera_pitch_deg)
            spinner.stop()
            print("Rendering complete")

            if output_path is None:
                timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
                body_name = observer.body.name.lower()
                output_path = f"output/{timestamp}-{body_name}.png"

            output_file = Path(output_path)
            output_file.parent.mkdir(parents=True, exist_ok=True)

            embed_metadata(image_array, scene, str(output_file))

            if print_caption:
                caption = generate_caption(scene, atmospheric_profile)
                print(caption)

            return 0

        error = SunsetError(
            "Could not find a body with a visible sunset at this time",
            suggestions=[
                "Try a different UTC time when a sunset might be occurring somewhere",
                "Specify a particular body with --body to see detailed error messages",
                "Check that the time is not during a period where the body is in continuous daylight or darkness",
            ],
        )
        return handle_error(error)

    except Exception as e:
        return handle_error(e, "rendering sunset scene")


def main():
    """CLI entry point."""
    args = parse_args()

    exit_code = render_sunset(
        utc_time=args.utc_time,
        body_id=args.body,
        random_seed=args.seed,
        output_path=args.output,
        print_caption=args.caption,
        camera_pitch_deg=args.camera_pitch,
    )

    sys.exit(exit_code)


if __name__ == "__main__":
    main()

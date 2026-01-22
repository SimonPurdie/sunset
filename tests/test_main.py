import pytest
from unittest.mock import patch, MagicMock
from pathlib import Path
import numpy as np
import tempfile
from PIL import Image

from sunset.main import (
    parse_args,
    get_body_order,
    create_scene,
    validate_scene,
    render_sunset,
    main,
)
from sunset.models import Observer, Scene
from sunset import __version__


def test_get_body_order_with_seed():
    """Test that get_body_order returns deterministic order with seed."""
    seed = 12345
    order1 = get_body_order(seed)
    order2 = get_body_order(seed)

    assert order1 == order2
    assert len(order1) == 6
    assert set(order1) == {"mercury", "venus", "earth", "moon", "mars", "titan"}


def test_get_body_order_without_seed():
    """Test that get_body_order works with None seed."""
    order = get_body_order(None)

    assert len(order) == 6
    assert set(order) == {"mercury", "venus", "earth", "moon", "mars", "titan"}


def test_get_body_order_different_seeds():
    """Test that different seeds produce different orders."""
    order1 = get_body_order(123)
    order2 = get_body_order(456)

    assert order1 != order2


def test_create_scene():
    """Test create_scene creates correct Scene object."""
    from sunset.models import SOLAR_SYSTEM_BODIES

    body = SOLAR_SYSTEM_BODIES["earth"]
    observer = Observer(
        latitude=45.0,
        longitude=-75.0,
        altitude_m=0.0,
        body=body,
        solar_elevation_deg=-0.3,
        sun_direction=(0.1, 0.2, 0.9),
        solar_angular_diameter_deg=0.5,
    )

    utc_time = "2026-01-20T12:00:00Z"
    random_seed = 42
    renderer_id = f"sunset-{__version__}"

    scene = create_scene(observer, utc_time, random_seed, renderer_id)

    assert scene.body == body
    assert scene.utc_time == utc_time
    assert scene.latitude == 45.0
    assert scene.longitude == -75.0
    assert scene.altitude_m == 0.0
    assert scene.random_seed == 42
    assert scene.renderer_id == renderer_id


def test_validate_scene_valid():
    """Test validate_scene with valid scene."""
    from sunset.models import SOLAR_SYSTEM_BODIES

    body = SOLAR_SYSTEM_BODIES["earth"]
    observer = Observer(
        latitude=45.0,
        longitude=-75.0,
        altitude_m=0.0,
        body=body,
        solar_elevation_deg=-0.3,
        sun_direction=(0.1, 0.2, 0.9),
        solar_angular_diameter_deg=0.5,
    )

    scene = Scene(
        body=body,
        utc_time="2026-01-20T12:00:00Z",
        latitude=45.0,
        longitude=-75.0,
        altitude_m=0.0,
        random_seed=42,
        renderer_id="test",
    )

    assert validate_scene(scene, observer) is True


def test_validate_scene_solar_elevation_too_high():
    """Test validate_scene with solar elevation too high."""
    from sunset.models import SOLAR_SYSTEM_BODIES

    body = SOLAR_SYSTEM_BODIES["earth"]
    observer = Observer(
        latitude=45.0,
        longitude=-75.0,
        altitude_m=0.0,
        body=body,
        solar_elevation_deg=1.0,
        sun_direction=(0.1, 0.2, 0.9),
        solar_angular_diameter_deg=0.5,
    )

    scene = Scene(
        body=body,
        utc_time="2026-01-20T12:00:00Z",
        latitude=45.0,
        longitude=-75.0,
        altitude_m=0.0,
        random_seed=42,
        renderer_id="test",
    )

    assert validate_scene(scene, observer) is False


def test_validate_scene_solar_elevation_too_low():
    """Test validate_scene with solar elevation too low."""
    from sunset.models import SOLAR_SYSTEM_BODIES

    body = SOLAR_SYSTEM_BODIES["earth"]
    observer = Observer(
        latitude=45.0,
        longitude=-75.0,
        altitude_m=0.0,
        body=body,
        solar_elevation_deg=-2.0,
        sun_direction=(0.1, 0.2, 0.9),
        solar_angular_diameter_deg=0.5,
    )

    scene = Scene(
        body=body,
        utc_time="2026-01-20T12:00:00Z",
        latitude=45.0,
        longitude=-75.0,
        altitude_m=0.0,
        random_seed=42,
        renderer_id="test",
    )

    assert validate_scene(scene, observer) is False


def test_validate_scene_latitude_invalid():
    """Test validate_scene with invalid latitude."""
    from sunset.models import SOLAR_SYSTEM_BODIES

    body = SOLAR_SYSTEM_BODIES["earth"]
    observer = Observer(
        latitude=95.0,
        longitude=-75.0,
        altitude_m=0.0,
        body=body,
        solar_elevation_deg=-0.3,
        sun_direction=(0.1, 0.2, 0.9),
        solar_angular_diameter_deg=0.5,
    )

    scene = Scene(
        body=body,
        utc_time="2026-01-20T12:00:00Z",
        latitude=95.0,
        longitude=-75.0,
        altitude_m=0.0,
        random_seed=42,
        renderer_id="test",
    )

    assert validate_scene(scene, observer) is False


def test_validate_scene_longitude_invalid():
    """Test validate_scene with invalid longitude."""
    from sunset.models import SOLAR_SYSTEM_BODIES

    body = SOLAR_SYSTEM_BODIES["earth"]
    observer = Observer(
        latitude=45.0,
        longitude=-200.0,
        altitude_m=0.0,
        body=body,
        solar_elevation_deg=-0.3,
        sun_direction=(0.1, 0.2, 0.9),
        solar_angular_diameter_deg=0.5,
    )

    scene = Scene(
        body=body,
        utc_time="2026-01-20T12:00:00Z",
        latitude=45.0,
        longitude=-200.0,
        altitude_m=0.0,
        random_seed=42,
        renderer_id="test",
    )

    assert validate_scene(scene, observer) is False


def test_validate_scene_altitude_negative():
    """Test validate_scene with negative altitude."""
    from sunset.models import SOLAR_SYSTEM_BODIES

    body = SOLAR_SYSTEM_BODIES["earth"]
    observer = Observer(
        latitude=45.0,
        longitude=-75.0,
        altitude_m=-100.0,
        body=body,
        solar_elevation_deg=-0.3,
        sun_direction=(0.1, 0.2, 0.9),
        solar_angular_diameter_deg=0.5,
    )

    scene = Scene(
        body=body,
        utc_time="2026-01-20T12:00:00Z",
        latitude=45.0,
        longitude=-75.0,
        altitude_m=-100.0,
        random_seed=42,
        renderer_id="test",
    )

    assert validate_scene(scene, observer) is False


@patch("sunset.main.render_scene")
@patch("sunset.main.embed_metadata")
@patch("sunset.main.generate_caption")
@patch("sunset.main.resolve_visibility_elevation")
@patch("sunset.main.get_body_profile")
@patch("sunset.main.resolve_sunset_location")
def test_render_sunset_success_with_body(
    mock_resolve, mock_profile, mock_visibility, mock_caption, mock_embed, mock_render
):
    """Test render_sunset successfully renders with specified body."""
    from sunset.models import SOLAR_SYSTEM_BODIES
    from sunset.models.atmosphere import NullAtmosphere

    body = SOLAR_SYSTEM_BODIES["earth"]
    observer = Observer(
        latitude=45.0,
        longitude=-75.0,
        altitude_m=0.0,
        body=body,
        solar_elevation_deg=-0.3,
        sun_direction=(0.1, 0.2, 0.9),
        solar_angular_diameter_deg=0.5,
    )

    mock_resolve.return_value = observer
    mock_profile.return_value = NullAtmosphere(body_name="Earth")
    mock_visibility.return_value = (0.0, "Justification", True)
    mock_render.return_value = np.zeros((512, 1024, 3), dtype=np.uint16)

    result = render_sunset(
        utc_time="2026-01-20T18:00:00Z",
        body_id="earth",
        random_seed=42,
        output_path="/tmp/test.png",
        print_caption=False,
    )

    assert result == 0
    mock_resolve.assert_called_once()
    mock_embed.assert_called_once()


@patch("sunset.main.render_scene")
@patch("sunset.main.embed_metadata")
@patch("sunset.main.generate_caption")
@patch("sunset.main.resolve_visibility_elevation")
@patch("sunset.main.get_body_profile")
@patch("sunset.main.resolve_sunset_location")
def test_render_sunset_print_caption(
    mock_resolve,
    mock_profile,
    mock_visibility,
    mock_caption,
    mock_embed,
    mock_render,
    capsys,
):
    """Test render_sunset with caption printing."""
    from sunset.models import SOLAR_SYSTEM_BODIES
    from sunset.models.atmosphere import NullAtmosphere

    body = SOLAR_SYSTEM_BODIES["earth"]
    observer = Observer(
        latitude=45.0,
        longitude=-75.0,
        altitude_m=0.0,
        body=body,
        solar_elevation_deg=-0.3,
        sun_direction=(0.1, 0.2, 0.9),
        solar_angular_diameter_deg=0.5,
    )

    mock_resolve.return_value = observer
    mock_profile.return_value = NullAtmosphere(body_name="Earth")
    mock_visibility.return_value = (0.0, "Justification", True)
    mock_render.return_value = np.zeros((512, 1024, 3), dtype=np.uint16)
    mock_caption.return_value = "Sunset on Earth at 45.0000°N, 75.0000°W at 2026-01-20T18:00:00Z. No atmosphere."

    result = render_sunset(
        utc_time="2026-01-20T18:00:00Z",
        body_id="earth",
        random_seed=42,
        output_path="/tmp/test.png",
        print_caption=True,
    )

    assert result == 0
    captured = capsys.readouterr()
    assert "Sunset on Earth" in captured.out


@patch("sunset.main.render_scene")
@patch("sunset.main.embed_metadata")
@patch("sunset.main.generate_caption")
@patch("sunset.main.resolve_visibility_elevation")
@patch("sunset.main.get_body_profile")
@patch("sunset.main.resolve_sunset_location")
def test_render_sunset_verbose(
    mock_resolve,
    mock_profile,
    mock_visibility,
    mock_caption,
    mock_embed,
    mock_render,
    capsys,
):
    """Test render_sunset with verbose flag."""
    from sunset.models import SOLAR_SYSTEM_BODIES
    from sunset.atmosphere.profiles import earth_profile

    body = SOLAR_SYSTEM_BODIES["earth"]
    observer = Observer(
        latitude=45.0,
        longitude=-75.0,
        altitude_m=0.0,
        body=body,
        solar_elevation_deg=-0.3,
        sun_direction=(0.1, 0.2, 0.9),
        solar_angular_diameter_deg=0.5,
    )

    mock_resolve.return_value = observer
    mock_profile.return_value = earth_profile
    mock_visibility.return_value = (0.0, "Justification", True)
    mock_render.return_value = np.zeros((512, 1024, 3), dtype=np.uint16)

    result = render_sunset(
        utc_time="2026-01-20T18:00:00Z",
        body_id="earth",
        random_seed=42,
        output_path="/tmp/test.png",
        print_caption=False,
        verbose=True,
    )

    assert result == 0
    captured = capsys.readouterr()
    assert "=== VERBOSE: Internal State ===" in captured.out
    assert "Body Information:" in captured.out
    assert "Atmospheric Profile:" in captured.out
    assert "Observer Geometry:" in captured.out
    assert "Scene Metadata:" in captured.out


@patch("sunset.main.resolve_sunset_location")
def test_render_sunset_unknown_body(mock_resolve):
    """Test render_sunset with unknown body."""
    mock_resolve.side_effect = ValueError("Unknown body: unknown_planet")

    result = render_sunset(
        utc_time="2026-01-20T18:00:00Z",
        body_id="unknown_planet",
        random_seed=42,
    )

    assert result == 1


@patch("sunset.main.resolve_sunset_location")
@patch("sunset.main.get_body_profile")
@patch("sunset.main.resolve_visibility_elevation")
def test_render_sunset_body_not_visible(mock_resolve, mock_profile, mock_visibility):
    """Test render_sunset when sun is not visible on specified body."""
    from sunset.models import SOLAR_SYSTEM_BODIES
    from sunset.models.atmosphere import NullAtmosphere

    body = SOLAR_SYSTEM_BODIES["earth"]
    observer = Observer(
        latitude=45.0,
        longitude=-75.0,
        altitude_m=0.0,
        body=body,
        solar_elevation_deg=-0.3,
        sun_direction=(0.1, 0.2, 0.9),
        solar_angular_diameter_deg=0.5,
    )

    mock_resolve.return_value = observer
    mock_profile.return_value = NullAtmosphere(body_name="Earth")
    mock_visibility.return_value = (0.0, "Not visible", False)

    result = render_sunset(
        utc_time="2026-01-20T18:00:00Z",
        body_id="earth",
        random_seed=42,
    )

    assert result == 1


def test_parse_args_defaults():
    """Test parse_args with default arguments."""
    with patch("sys.argv", ["sunset"]):
        args = parse_args()

        assert args.utc_time.endswith("Z")
        assert args.body is None
        assert args.seed is None
        assert args.output is None
        assert args.caption is False


def test_parse_args_with_all_options():
    """Test parse_args with all options specified."""
    with patch(
        "sys.argv",
        [
            "sunset",
            "--utc-time",
            "2026-01-20T12:00:00Z",
            "--body",
            "mars",
            "--seed",
            "123",
            "--output",
            "output.png",
            "--caption",
        ],
    ):
        args = parse_args()

        assert args.utc_time == "2026-01-20T12:00:00Z"
        assert args.body == "mars"
        assert args.seed == 123
        assert args.output == "output.png"
        assert args.caption is True
        assert args.verbose is False


def test_parse_args_verbose():
    """Test parse_args with verbose flag."""
    with patch(
        "sys.argv",
        ["sunset", "--verbose"],
    ):
        args = parse_args()

        assert args.verbose is True


def test_parse_args_verbose_short():
    """Test parse_args with verbose short flag."""
    with patch(
        "sys.argv",
        ["sunset", "-v"],
    ):
        args = parse_args()

        assert args.verbose is True


def test_parse_args_invalid_body():
    """Test parse_args with invalid body choice."""
    with patch("sys.argv", ["sunset", "--body", "invalid"]):
        with pytest.raises(SystemExit):
            parse_args()


def test_deterministic_output():
    """Test that same inputs produce identical PNG output."""
    utc_time = "2026-01-20T18:00:00Z"
    seed = 42

    with tempfile.TemporaryDirectory() as tmpdir:
        output1 = Path(tmpdir) / "sunset1.png"
        output2 = Path(tmpdir) / "sunset2.png"

        result1 = render_sunset(
            utc_time=utc_time,
            body_id="earth",
            random_seed=seed,
            output_path=str(output1),
            print_caption=False,
        )
        assert result1 == 0, "First render failed"

        result2 = render_sunset(
            utc_time=utc_time,
            body_id="earth",
            random_seed=seed,
            output_path=str(output2),
            print_caption=False,
        )
        assert result2 == 0, "Second render failed"

        img1 = Image.open(output1)
        img2 = Image.open(output2)

        assert np.array_equal(np.array(img1), np.array(img2)), (
            "Images are not identical"
        )

        exif1 = img1.info
        exif2 = img2.info

        assert exif1 == exif2, "PNG metadata differs"


def test_deterministic_output_without_explicit_seed():
    """Test that same inputs produce identical output when seed is not provided."""
    utc_time = "2026-01-20T18:00:00Z"

    with tempfile.TemporaryDirectory() as tmpdir:
        output1 = Path(tmpdir) / "sunset1.png"
        output2 = Path(tmpdir) / "sunset2.png"

        result1 = render_sunset(
            utc_time=utc_time,
            body_id="earth",
            random_seed=None,
            output_path=str(output1),
            print_caption=False,
        )
        assert result1 == 0, "First render failed"

        result2 = render_sunset(
            utc_time=utc_time,
            body_id="earth",
            random_seed=None,
            output_path=str(output2),
            print_caption=False,
        )
        assert result2 == 0, "Second render failed"

        img1 = Image.open(output1)
        img2 = Image.open(output2)

        assert np.array_equal(np.array(img1), np.array(img2)), (
            "Images are not identical when seed not provided"
        )


def test_seed_affects_location():
    """Test that different seeds produce different outputs."""
    utc_time = "2026-01-20T18:00:00Z"

    with tempfile.TemporaryDirectory() as tmpdir:
        output1 = Path(tmpdir) / "sunset1.png"
        output2 = Path(tmpdir) / "sunset2.png"

        result1 = render_sunset(
            utc_time=utc_time,
            body_id="earth",
            random_seed=123,
            output_path=str(output1),
            print_caption=False,
        )
        assert result1 == 0, "First render failed"

        result2 = render_sunset(
            utc_time=utc_time,
            body_id="earth",
            random_seed=456,
            output_path=str(output2),
            print_caption=False,
        )
        assert result2 == 0, "Second render failed"

        img1 = Image.open(output1)
        img2 = Image.open(output2)

        assert not np.array_equal(np.array(img1), np.array(img2)), (
            "Images should differ with different seeds"
        )


def test_seed_preserved_metadata():
    """Test that seed is recorded in PNG metadata."""
    utc_time = "2026-01-20T18:00:00Z"
    seed = 789

    with tempfile.TemporaryDirectory() as tmpdir:
        output = Path(tmpdir) / "sunset.png"

        result = render_sunset(
            utc_time=utc_time,
            body_id="earth",
            random_seed=seed,
            output_path=str(output),
            print_caption=False,
        )
        assert result == 0, "Render failed"

        img = Image.open(output)
        metadata = img.info

        assert "random_seed" in metadata, "random_seed not in metadata"
        assert int(metadata["random_seed"]) == seed, (
            "Seed in metadata does not match input seed"
        )

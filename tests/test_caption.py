import pytest

from sunset.caption.generator import generate_caption
from sunset.models.atmosphere import AtmosphericProfile, NullAtmosphere
from sunset.models.bodies import Body
from sunset.models.scene import Scene


class TestAltitudeMentioned:
    """Test that altitude is mentioned when not at surface."""

    def test_surface_altitude_not_mentioned(self):
        scene = Scene(
            body=Body(name="Earth", radius_m=6_371_000, has_atmosphere=True),
            utc_time="2026-01-21T12:00:00Z",
            latitude=37.7749,
            longitude=-122.4194,
            altitude_m=0.0,
            random_seed=12345,
            renderer_id="test",
        )

        caption = generate_caption(scene)
        assert "meters" not in caption.lower()
        assert "above surface" not in caption.lower()

    def test_positive_altitude_mentioned(self):
        scene = Scene(
            body=Body(name="Earth", radius_m=6_371_000, has_atmosphere=True),
            utc_time="2026-01-21T12:00:00Z",
            latitude=37.7749,
            longitude=-122.4194,
            altitude_m=5000.0,
            random_seed=12345,
            renderer_id="test",
        )

        caption = generate_caption(scene)
        assert "5000 meters" in caption
        assert "above surface" in caption

    def test_fractional_altitude_rounded(self):
        scene = Scene(
            body=Body(name="Earth", radius_m=6_371_000, has_atmosphere=True),
            utc_time="2026-01-21T12:00:00Z",
            latitude=37.7749,
            longitude=-122.4194,
            altitude_m=5123.7,
            random_seed=12345,
            renderer_id="test",
        )

        caption = generate_caption(scene)
        assert "5124 meters" in caption


class TestNoAtmosphereMentioned:
    """Test that lack of atmosphere is mentioned for airless bodies."""

    def test_airless_body_with_null_atmosphere(self):
        scene = Scene(
            body=Body(name="Moon", radius_m=1_737_400, has_atmosphere=False),
            utc_time="2026-01-21T12:00:00Z",
            latitude=0.0,
            longitude=0.0,
            altitude_m=0.0,
            random_seed=12345,
            renderer_id="test",
        )
        atmosphere = NullAtmosphere(body_name="Moon")

        caption = generate_caption(scene, atmosphere)
        assert "No atmosphere" in caption

    def test_airless_body_without_atmosphere_arg(self):
        scene = Scene(
            body=Body(name="Mercury", radius_m=2_439_700, has_atmosphere=False),
            utc_time="2026-01-21T12:00:00Z",
            latitude=0.0,
            longitude=0.0,
            altitude_m=0.0,
            random_seed=12345,
            renderer_id="test",
        )

        caption = generate_caption(scene, atmosphere=None)
        assert "No atmosphere" in caption

    def test_body_with_atmosphere_shows_composition(self):
        scene = Scene(
            body=Body(name="Mars", radius_m=3_389_500, has_atmosphere=True),
            utc_time="2026-01-21T12:00:00Z",
            latitude=0.0,
            longitude=0.0,
            altitude_m=0.0,
            random_seed=12345,
            renderer_id="test",
        )
        atmosphere = AtmosphericProfile(
            body_name="Mars",
            gas_composition={"CO2": 0.953, "N2": 0.027, "Ar": 0.016, "O2": 0.004},
            rayleigh_coefficient=lambda w: w**-4,
            density_profile=lambda a: 0.02 * 2.718 ** (-a / 11100),
            pressure_profile=lambda a: 610 * 2.718 ** (-a / 11100),
            refractive_index=1.00029,
        )

        caption = generate_caption(scene, atmosphere)
        assert "Atmosphere primarily composed of CO2" in caption


class TestNoSpeculation:
    """Test that no speculative or fictional claims are made."""

    def test_caption_uses_only_metadata(self):
        scene = Scene(
            body=Body(name="Earth", radius_m=6_371_000, has_atmosphere=True),
            utc_time="2026-01-21T12:00:00Z",
            latitude=37.7749,
            longitude=-122.4194,
            altitude_m=0.0,
            random_seed=12345,
            renderer_id="test",
        )

        caption = generate_caption(scene)

        speculative_words = [
            "beautiful",
            "stunning",
            "amazing",
            "breathtaking",
            "magnificent",
            "perhaps",
            "might",
            "possibly",
            "probably",
            "likely",
            "wonderful",
        ]

        for word in speculative_words:
            assert word not in caption.lower()

    def test_caption_ends_with_period(self):
        scene = Scene(
            body=Body(name="Earth", radius_m=6_371_000, has_atmosphere=True),
            utc_time="2026-01-21T12:00:00Z",
            latitude=0.0,
            longitude=0.0,
            altitude_m=0.0,
            random_seed=12345,
            renderer_id="test",
        )

        caption = generate_caption(scene)
        assert caption.endswith(".")

    def test_caption_includes_body_name(self):
        scene = Scene(
            body=Body(name="Venus", radius_m=6_051_800, has_atmosphere=True),
            utc_time="2026-01-21T12:00:00Z",
            latitude=0.0,
            longitude=0.0,
            altitude_m=0.0,
            random_seed=12345,
            renderer_id="test",
        )

        caption = generate_caption(scene)
        assert "Venus" in caption

    def test_caption_includes_utc_time(self):
        scene = Scene(
            body=Body(name="Earth", radius_m=6_371_000, has_atmosphere=True),
            utc_time="2026-01-21T18:30:45Z",
            latitude=0.0,
            longitude=0.0,
            altitude_m=0.0,
            random_seed=12345,
            renderer_id="test",
        )

        caption = generate_caption(scene)
        assert "2026-01-21T18:30:45Z" in caption

    def test_caption_includes_location(self):
        scene = Scene(
            body=Body(name="Earth", radius_m=6_371_000, has_atmosphere=True),
            utc_time="2026-01-21T12:00:00Z",
            latitude=40.7128,
            longitude=-74.0060,
            altitude_m=0.0,
            random_seed=12345,
            renderer_id="test",
        )

        caption = generate_caption(scene)
        assert "40.7128" in caption
        assert "N" in caption
        assert "74.0060" in caption
        assert "W" in caption

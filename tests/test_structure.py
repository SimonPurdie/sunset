"""Basic test to verify test infrastructure is working."""


def test_project_structure():
    """Verify that the project structure is set up correctly."""
    import sunset

    assert hasattr(sunset, "__version__")
    assert sunset.__version__ == "0.1.0"

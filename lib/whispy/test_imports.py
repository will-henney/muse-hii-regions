def test_relative_imports():
    """
    I have restuctured these files as a package, so we have to use relative imports.  This is to check that they work.
    """
    from . import extract
    from . import moments
    from . import sky
    from . import masktools

    # This will need changing if I ever give this package a proper name
    assert __name__ == "lib.test_imports"
    # Check that a random function is where it should be
    assert "find_moments" in dir(moments)

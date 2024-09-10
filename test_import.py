# test_import.py

def test_import_package():
    try:
        import xraybinaryorbit  # Replace with the actual name of your package
    except ImportError as e:
        assert False, f"Failed to import the package: {e}"

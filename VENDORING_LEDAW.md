# LEDAW Vendoring Documentation

## What Was Done

LEDAW (LED Analysis Workflow) has been vendored (bundled) directly into the xbpy package to eliminate installation complexity and path configuration issues.

## Vendored Package Location

```
xbpy/
├── _vendor/
│   ├── __init__.py
│   ├── VENDOR_README.md
│   └── ledaw_package/        # Complete LEDAW package copied here
│       ├── ...
│       └── (all LEDAW files)
└── orcautil/
    └── led_extract.py         # Imports from xbpy._vendor.ledaw_package
```

## Why Vendoring?

### Problems Solved

1. **Segmentation Faults**: Dynamic import path searching was causing segfaults with LEDAW's compiled extensions
2. **Installation Complexity**: Users needed to manually configure `sys.path` before importing
3. **Reproducibility**: Different LEDAW installations/versions caused inconsistent behavior
4. **User Experience**: Required complex setup instructions before using LED functionality

### Benefits

- ✅ **Zero configuration**: Just `from xbpy.orcautil import compute_led_interaction_matrix`
- ✅ **No segfaults**: Clean import from known location
- ✅ **Reproducible**: Same LEDAW version bundled with xbpy
- ✅ **Portable**: Works anywhere xbpy works

## Attribution

**IMPORTANT**: LEDAW is developed and maintained by its original authors. All credit for LEDAW functionality belongs to the LEDAW development team. xbpy simply bundles a copy for convenience.

The vendored LEDAW remains under its original license (see `xbpy/_vendor/ledaw_package/` for details). xbpy's license does NOT apply to the vendored LEDAW code.

## Usage - Before vs After

### Before (Complex)

```python
import sys
# User had to add LEDAW to path before importing xbpy
sys.path.insert(0, "/path/to/LEDAW")

from xbpy.orcautil import compute_led_interaction_matrix

# This often caused segfaults during import path searching
```

### After (Simple)

```python
from xbpy.orcautil import compute_led_interaction_matrix

# That's it! LEDAW is bundled and imported automatically
```

## Implementation Details

### Import Strategy

```python
# In led_extract.py
try:
    from xbpy._vendor.ledaw_package import engine_LED_N_body as _LEDAW_ENGINE
    _LEDAW_AVAILABLE = True
except ImportError as e:
    _LEDAW_ENGINE = None
    _LEDAW_AVAILABLE = False
```

- Imports LEDAW at module level from vendored location
- Clean import from known path avoids segfaults
- Graceful degradation if import fails (with helpful error message)

### Function Usage

```python
def _run_ledaw(...):
    if not _LEDAW_AVAILABLE:
        raise ImportError("LEDAW not available - see error details")
    
    # Use vendored LEDAW
    _LEDAW_ENGINE(
        main_filenames=[...],
        ...
    )
```

## Updating Vendored LEDAW

To update the bundled LEDAW version:

1. **Obtain new LEDAW version**:
   ```bash
   # Get updated LEDAW from original source
   ```

2. **Replace vendored copy**:
   ```bash
   cd /path/to/xbpy/module/xbpy/_vendor
   rm -rf ledaw_package
   cp -r /path/to/new/LEDAW/ledaw_package .
   ```

3. **Test thoroughly**:
   ```bash
   # Test import
   python3 -c "from xbpy._vendor.ledaw_package import engine_LED_N_body"
   
   # Test full LED workflow
   python3 -c "from xbpy.orcautil import compute_led_interaction_matrix"
   ```

4. **Update documentation**:
   - Update version info in `_vendor/VENDOR_README.md`
   - Document any API changes in `orcautil/LED_README.md`
   - Update this file with any new considerations

5. **Commit and release**:
   ```bash
   git add xbpy/_vendor/ledaw_package
   git commit -m "Update vendored LEDAW to version X.Y.Z"
   ```

## Compatibility Notes

### Python Version

LEDAW contains compiled C/Fortran extensions. The bundled version must be compatible with:
- Python 3.11+ (current xbpy target)
- Linux/macOS/Windows (depending on compiled binaries)

### System Dependencies

LEDAW may require system libraries:
- BLAS (Basic Linear Algebra Subprograms)
- LAPACK (Linear Algebra PACKage)
- MKL (Math Kernel Library) - optional

If vendored LEDAW fails to import, these dependencies may be missing.

### Recompilation

If bundled LEDAW is incompatible with a user's system:

1. **User rebuilds LEDAW**:
   ```bash
   cd /path/to/LEDAW/source
   python3 setup.py build_ext --inplace
   ```

2. **Replace vendored copy**:
   ```bash
   cp -r /path/to/rebuilt/ledaw_package /path/to/xbpy/module/xbpy/_vendor/
   ```

## Legal Considerations

### License Compliance

- LEDAW remains under its original license
- xbpy's license applies only to xbpy code, NOT vendored dependencies
- Users must comply with LEDAW's license terms
- See `_vendor/VENDOR_README.md` for full details

### Distribution

When distributing xbpy:
- Include full `_vendor/` directory with LEDAW
- Include `_vendor/VENDOR_README.md` with attribution
- Ensure LEDAW license file is preserved
- Document that LEDAW is bundled, not developed by xbpy team

## Future Considerations

### Alternative Approaches

If vendoring becomes problematic:

1. **Optional dependency**: Make LEDAW optional, fall back to external
2. **Dynamic compilation**: Compile LEDAW on first use for user's system
3. **Package separately**: Create `xbpy-led` extension package
4. **Binary wheels**: Distribute pre-compiled wheels for common platforms

### Maintenance Burden

Vendoring adds maintenance responsibility:
- Keep LEDAW version up to date
- Test compatibility with new Python versions
- Handle platform-specific issues
- Respond to LEDAW bugs/updates

Balance against user experience improvements.

## Testing Vendored LEDAW

### Unit Tests

```python
def test_vendored_ledaw_import():
    """Test that vendored LEDAW can be imported."""
    from xbpy._vendor.ledaw_package import engine_LED_N_body
    assert callable(engine_LED_N_body)

def test_led_extraction_with_vendored():
    """Test LED extraction uses vendored LEDAW."""
    from xbpy.orcautil import compute_led_interaction_matrix
    # ... test with actual ORCA files ...
```

### Integration Tests

```bash
# Test on clean Python environment
python3 -m venv test_env
source test_env/bin/activate
pip install /path/to/xbpy
python3 -c "from xbpy.orcautil import compute_led_interaction_matrix; print('OK')"
```

## References

- LEDAW Original Source: [Add URL if known]
- CovaLED Methodology: https://www.faccts.de/docs/opi/nightly/docs/contents/notebooks/covaled_ethane.html
- Python Vendoring Best Practices: https://caremad.io/posts/2013/07/setup-vs-requirement/





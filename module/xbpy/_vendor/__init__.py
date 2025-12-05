"""
Vendored third-party dependencies.

This directory contains bundled third-party packages to ensure compatibility
and avoid installation/path issues.

IMPORTANT: These are NOT maintained by xbpy developers. See VENDOR_README.md
for information about each vendored package.
"""

# Make vendored packages available
from . import ledaw_package

__all__ = ['ledaw_package']





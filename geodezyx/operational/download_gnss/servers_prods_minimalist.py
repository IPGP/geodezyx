#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GNSS Products Server Dispatcher - Minimalist Config-Driven Version

Replaces the original servers_prods.py (86 lines with 9 hardcoded servers)
with a unified, config-driven approach using products_config.yaml.

All archive centers are now 100% config-driven, loaded from YAML.
No special-case functions remain.

This version maintains 100% backward compatibility with the original
_server_select_prods interface used by download_prods.py.

Version: 2.0 - 2024
"""

import logging

from .products_config_loader import get_config, ProductsConfigError

log = logging.getLogger("geodezyx")


def _server_select_prods(archive_center, mgex=False, repro=0):
    """
    Resolve archive center name to FTP host, base directory, and protocol.

    This is the main dispatcher function called by download_prods.py.
    It has been refactored to use the config-driven approach via YAML.

    Parameters
    ----------
    archive_center : str
        Name of the IGS archive/data center (e.g., "igs", "cddis", "esa").
        Supports aliases for backward compatibility.
    mgex : bool, optional
        Get MGEX (multi-constellation) products. Default is False.
    repro : int, optional
        Reprocessing campaign number. Default is 0 (operational).
        Valid values: 0 (operational), 1-3 (reprocessing campaigns).

    Returns
    -------
    host : str
        FTP server hostname (e.g., "igs.ign.fr")
    basedir : str
        Base directory on the server (e.g., "/pub/igs/products/")
    protocol : str
        Protocol type: "ftp" or "sftp"
    secure_ftp : bool
        Whether SFTP (secure) connection is required

    Returns
    -------
    tuple of (None, None, None, None) if archive center not found.

    Examples
    --------
    >>> host, basedir, protocol, secure = _server_select_prods("igs")
    >>> host, basedir, protocol, secure = _server_select_prods("cddis", mgex=True, repro=2)

    Notes
    -----
    All archive centers are 100% configuration-driven via products_config.yaml.
    No special-case functions are needed—the YAML approach is sufficiently
    flexible to handle all 9 archive centers plus MGEX and reprocessing variants.

    Backward Compatibility
    ----------------------
    This function maintains the exact same interface and behavior as the
    original servers_prods.py, so it can be used as a drop-in replacement
    in download_prods.py without any code changes.
    """
    try:
        cfg = get_config()
        host, basedir, protocol, secure_ftp = cfg.get_archive_center(
            archive_center, mgex=mgex, repro=repro
        )
        
        if host is None:
            log.error(
                "Unknown archive center: %s (available: %s)",
                archive_center,
                ", ".join(cfg.list_archive_centers())
            )
            return None, None, None, None

        return host, basedir, protocol, secure_ftp

    except ProductsConfigError as e:
        log.error("Configuration error: %s", e)
        return None, None, None, None

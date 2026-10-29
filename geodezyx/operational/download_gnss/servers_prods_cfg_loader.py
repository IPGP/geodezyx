#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Products Configuration Loader

Provides a flexible API for loading and querying GNSS product archive
center configurations from YAML, replacing the hardcoded
servers_prods.py functions.

This module is analogous to servers_rinex_cfg_loader.py but specifically
for GNSS product archive centers (sp3, clk, erp, etc.).
"""

import logging
import os
from pathlib import Path

try:
    import yaml
except ImportError:
    yaml = None

log = logging.getLogger("geodezyx")


class ProductsConfigError(Exception):
    """Exception raised for products configuration errors."""
    pass


class ProductsConfig:
    """
    Load and query GNSS product archive center configurations from YAML.

    This class provides a unified interface for accessing product server
    configurations, replacing the hardcoded if/elif chain in
    _server_select_prods().

    Example:
        >>> cfg = ProductsConfig()
        >>> host, basedir, protocol, secure = cfg.get_archive_center("igs")
        >>> archive_centers = cfg.list_archive_centers()
        >>> all_product_types = cfg.list_product_types()
    """

    _instance = None  # Singleton instance
    _config = None

    def __init__(self, config_path=None):
        """
        Initialize the products configuration loader.

        Parameters
        ----------
        config_path : str, optional
            Path to servers_prods_cfg.yaml file.
            If None, searches relative to this module.
        """
        if yaml is None:
            raise ProductsConfigError(
                "PyYAML is required. Install with: pip install pyyaml"
            )

        if config_path is None:
            # Default to file in same directory as this module
            module_dir = Path(__file__).parent
            config_path = module_dir / "servers_prods_cfg.yaml"
        else:
            config_path = Path(config_path)

        if not config_path.exists():
            raise ProductsConfigError(
                f"Products configuration not found: {config_path}"
            )

        self.config_path = config_path
        self._load_config()

    def _load_config(self):
        """Load YAML configuration file."""
        with open(self.config_path, 'r') as f:
            self._config = yaml.safe_load(f)

        if not self._config:
            raise ProductsConfigError(
                f"Failed to load configuration from {self.config_path}"
            )

        log.debug(
            "Loaded products configuration: %s (%d archive centers)",
            self.config_path,
            len(self.get_archive_centers())
        )

    def get_config(self):
        """
        Get the entire configuration dictionary.

        Returns
        -------
        dict
            Complete configuration object
        """
        return self._config

    def get_archive_center(self, archive_center_name, mgex=False, repro=0):
        """
        Get archive center configuration (replaces _server_select_prods).

        Parameters
        ----------
        archive_center_name : str
            Name/alias of the archive center.
        mgex : bool, optional
            Get MGEX (multi-constellation) products. Default is False.
        repro : int, optional
            Reprocessing campaign (0=operational, 1/2/3=repro campaigns).
            Default is 0.

        Returns
        -------
        tuple
            (host, basedir, protocol, secure_ftp) matching _server_select_prods
            output format for backward compatibility.
            Returns (None, None, None, None) if not found.
        """
        # Resolve alias if needed
        name = self._resolve_alias(archive_center_name)

        if name not in self.get_archive_centers():
            log.error("Unknown archive center: %s", archive_center_name)
            return None, None, None, None

        center = self._config['products'][name]

        host = center.get('host')
        protocol = center.get('protocol', 'ftp')
        secure_ftp = center.get('secure', False)

        # Build basedir with MGEX and repro suffixes
        basedir = center.get('basedir', '')

        if mgex and 'mgex_suffix' in center:
            basedir = os.path.join(basedir, center['mgex_suffix'])

        if repro and repro in self._config.get('reprocessing', {}):
            repro_subdir = self._config['reprocessing'][f'repro{repro}'].get('subdir', '')
            basedir = os.path.join(basedir, repro_subdir)
        elif repro > 0:
            log.warning("Unknown reprocessing campaign: %d, using operational", repro)

        # Ensure consistent trailing slashes
        basedir = basedir.rstrip('/') + '/'

        return host, basedir, protocol, secure_ftp

    def _resolve_alias(self, name):
        """Resolve archive center alias to actual name."""
        aliases = self._config.get('aliases', {})
        return aliases.get(name, name)

    def get_archive_centers(self):
        """
        List all available archive center names.

        Returns
        -------
        dict
            Archive center configurations keyed by name
        """
        return self._config.get('products', {})

    def list_archive_centers(self):
        """
        List all available archive center names.

        Returns
        -------
        list
            Names of all archive centers
        """
        return list(self.get_archive_centers().keys())

    def get_archive_center_info(self, archive_center_name):
        """
        Get full configuration dict for an archive center.

        Parameters
        ----------
        archive_center_name : str
            Name or alias of archive center

        Returns
        -------
        dict
            Full configuration dict for the center, or None if not found
        """
        name = self._resolve_alias(archive_center_name)
        return self._config['products'].get(name)

    def get_product_types(self):
        """
        Get all product type configurations.

        Returns
        -------
        dict
            Product type information keyed by product code (e.g., "sp3", "clk")
        """
        return self._config.get('product_types', {})

    def list_product_types(self):
        """
        List all available product type codes.

        Returns
        -------
        list
            Product type codes (e.g., ["sp3", "clk", "erp"])
        """
        return list(self.get_product_types().keys())

    def get_product_info(self, product_code):
        """
        Get information about a product type.

        Parameters
        ----------
        product_code : str
            Product code (e.g., "sp3", "clk")

        Returns
        -------
        dict
            Product configuration dict, or None if not found
        """
        return self._config['product_types'].get(product_code)

    def get_analysis_centers(self):
        """
        Get all analysis center information.

        Returns
        -------
        dict
            Analysis center info keyed by abbreviation (e.g., "cod", "grg")
        """
        return self._config.get('analysis_centers', {})

    def list_analysis_centers(self):
        """
        List all available analysis center abbreviations.

        Returns
        -------
        list
            Analysis center abbreviations (e.g., ["cod", "grg", "wum"])
        """
        return list(self.get_analysis_centers().keys())

    def get_analysis_center_info(self, ac_abbr):
        """
        Get information about an analysis center.

        Parameters
        ----------
        ac_abbr : str
            Analysis center abbreviation (e.g., "cod")

        Returns
        -------
        dict
            Analysis center configuration dict, or None if not found
        """
        return self._config['analysis_centers'].get(ac_abbr)

    def get_aliases(self):
        """
        Get all alias mappings.

        Returns
        -------
        dict
            Alias → name mappings
        """
        return self._config.get('aliases', {})

    def uses_generic_builder(self, archive_center_name):
        """
        Check if an archive center is config-driven (uses_generic: true).

        In this implementation, all archive centers are config-driven.

        Parameters
        ----------
        archive_center_name : str
            Name or alias of archive center

        Returns
        -------
        bool
            True if uses generic builder, False if requires special-case code
        """
        name = self._resolve_alias(archive_center_name)
        center = self._config['products'].get(name, {})
        return center.get('uses_generic', False)


# Global singleton instance for convenience
_products_config_instance = None


def get_config(config_path=None):
    """
    Get or create the global ProductsConfig singleton instance.

    Parameters
    ----------
    config_path : str, optional
        Path to servers_prods_cfg.yaml. Only used on first call.

    Returns
    -------
    ProductsConfig
        Singleton configuration instance
    """
    global _products_config_instance

    if _products_config_instance is None:
        _products_config_instance = ProductsConfig(config_path)

    return _products_config_instance


# For backward compatibility with original _server_select_prods interface
def server_select_prods(archive_center, mgex=False, repro=0):
    """
    Drop-in replacement for _server_select_prods function.

    Returns the same tuple format as the original function for
    100% backward compatibility.

    Parameters
    ----------
    archive_center : str
        Name of archive center
    mgex : bool, optional
        Get MGEX products
    repro : int, optional
        Reprocessing campaign

    Returns
    -------
    tuple
        (host, basedir, protocol, secure_ftp)
    """
    cfg = get_config()
    return cfg.get_archive_center(archive_center, mgex, repro)

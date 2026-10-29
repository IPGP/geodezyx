# -*- coding: utf-8 -*-
"""
Server configuration loader for GNSS RINEX downloads.

This module loads and manages server configurations from a YAML file,
providing a centralized way to manage server parameters without needing
to modify Python code.
"""

import logging
import os
from pathlib import Path
from typing import Dict, Any, Optional

try:
    import yaml
except ImportError:
    yaml = None

log = logging.getLogger("geodezyx")


class ServerCfgRinexError(Exception):
    """Raised when there's an error loading or accessing server configuration."""
    pass


class ServerCfgRinex:
    """
    Manages GNSS server configurations loaded from YAML file.

    This class provides centralized access to server parameters, supporting
    both direct server lookups and aliases.

    Parameters
    ----------
    config_file : str or Path, optional
        Path to the YAML configuration file. If None, uses default location.

    Raises
    ------
    ServerConfigError
        If YAML is not installed or config file cannot be loaded.
    """

    def __init__(self, config_file: Optional[str] = None):
        """Initialize the ServerConfig with optional custom config file path."""
        if yaml is None:
            raise ServerCfgRinexError(
                "PyYAML is required for server configuration. "
                "Install it with: pip install pyyaml"
            )

        self.config_file = config_file or self._get_default_config_path()
        self._config = None
        self._load_config()

    @staticmethod
    def _get_default_config_path() -> Path:
        """Get the default configuration file path."""
        current_dir = Path(__file__).parent
        config_path = current_dir / "servers_rinex_cfg.yaml"
        return config_path

    def _load_config(self) -> None:
        """Load configuration from YAML file."""
        config_path = Path(self.config_file)

        if not config_path.exists():
            raise ServerCfgRinexError(
                f"Server configuration file not found: {config_path}"
            )

        try:
            with open(config_path, 'r', encoding='utf-8') as f:
                self._config = yaml.safe_load(f)
        except yaml.YAMLError as e:
            raise ServerCfgRinexError(
                f"Error parsing YAML configuration: {e}"
            ) from e
        except IOError as e:
            raise ServerCfgRinexError(
                f"Error reading configuration file: {e}"
            ) from e

        if not self._config or 'rinex_servers' not in self._config:
            raise ServerCfgRinexError(
                "Invalid configuration: 'rinex_servers' section not found"
            )

        log.info(
            "Loaded server configuration from %s with %d servers",
            config_path,
            len(self._config.get('rinex_servers', {}))
        )

    def get_server(self, datacenter: str) -> Dict[str, Any]:
        """
        Get server configuration by name.

        Supports both direct server names and aliases.

        Parameters
        ----------
        datacenter : str
            Server name or alias (e.g., 'igs_cddis', 'igs', 'brdc')

        Returns
        -------
        dict
            Server configuration dictionary

        Raises
        ------
        ServerConfigError
            If server/alias is not found in configuration
        """
        # Resolve aliases
        actual_server = self._resolve_alias(datacenter)

        if actual_server not in self._config['rinex_servers']:
            raise ServerCfgRinexError(
                f"Unknown server: '{datacenter}'. "
                f"Available servers: {list(self._config['rinex_servers'].keys())}"
            )

        return self._config['rinex_servers'][actual_server]

    def _resolve_alias(self, name: str) -> str:
        """
        Resolve a server name or alias to the actual server name.

        Parameters
        ----------
        name : str
            Server name or alias

        Returns
        -------
        str
            Resolved server name
        """
        aliases = self._config.get('aliases', {})
        return aliases.get(name, name)

    def list_servers(self) -> list:
        """
        Get list of all available servers.

        Returns
        -------
        list
            List of server names
        """
        return list(self._config['rinex_servers'].keys())

    def list_aliases(self) -> Dict[str, str]:
        """
        Get all configured aliases.

        Returns
        -------
        dict
            Mapping of aliases to server names
        """
        return self._config.get('aliases', {})

    def get_observation_servers(self) -> list:
        """
        Get list of servers providing observation data.

        Returns
        -------
        list
            List of server names offering observation files
        """
        return [
            name for name, config in self._config['rinex_servers'].items()
            if config.get('file_type') == 'observation'
        ]

    def get_navigation_servers(self) -> list:
        """
        Get list of servers providing navigation data.

        Returns
        -------
        list
            List of server names offering navigation files
        """
        return [
            name for name, config in self._config['rinex_servers'].items()
            if config.get('file_type') == 'navigation'
        ]

    def uses_generic_builder(self, datacenter: str) -> bool:
        """
        Check if a server uses the generic URL builder.

        Parameters
        ----------
        datacenter : str
            Server name or alias

        Returns
        -------
        bool
            True if server uses _generic_server, False for custom builders
        """
        config = self.get_server(datacenter)
        return config.get('uses_generic', False)


# Global configuration instance (lazy-loaded)
_config_instance: Optional[ServerCfgRinex] = None


def get_cfg_rinex(config_file: Optional[str] = None) -> ServerCfgRinex:
    """
    Get the global server configuration instance (lazy-loaded).

    Parameters
    ----------
    config_file : str or Path, optional
        Path to YAML config file. Only used on first call.

    Returns
    -------
    ServerCfgRinex
        Global configuration instance
    """
    global _config_instance
    if _config_instance is None:
        _config_instance = ServerCfgRinex(config_file)
    return _config_instance


def reload_config(config_file: Optional[str] = None) -> ServerCfgRinex:
    """
    Reload the server configuration (useful for testing or config updates).

    Parameters
    ----------
    config_file : str or Path, optional
        Path to YAML config file

    Returns
    -------
    ServerCfgRinex
        Reloaded configuration instance
    """
    global _config_instance
    _config_instance = ServerCfgRinex(config_file)
    return _config_instance

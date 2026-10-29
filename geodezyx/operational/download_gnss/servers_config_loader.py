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


class ServerConfigError(Exception):
    """Raised when there's an error loading or accessing server configuration."""
    pass


class ServerConfig:
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
            raise ServerConfigError(
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
        config_path = current_dir / "servers_config.yaml"
        return config_path

    def _load_config(self) -> None:
        """Load configuration from YAML file."""
        config_path = Path(self.config_file)

        if not config_path.exists():
            raise ServerConfigError(
                f"Server configuration file not found: {config_path}"
            )

        try:
            with open(config_path, 'r', encoding='utf-8') as f:
                self._config = yaml.safe_load(f)
        except yaml.YAMLError as e:
            raise ServerConfigError(
                f"Error parsing YAML configuration: {e}"
            ) from e
        except IOError as e:
            raise ServerConfigError(
                f"Error reading configuration file: {e}"
            ) from e

        if not self._config or 'servers' not in self._config:
            raise ServerConfigError(
                "Invalid configuration: 'servers' section not found"
            )

        log.info(
            "Loaded server configuration from %s with %d servers",
            config_path,
            len(self._config.get('servers', {}))
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

        if actual_server not in self._config['servers']:
            raise ServerConfigError(
                f"Unknown server: '{datacenter}'. "
                f"Available servers: {list(self._config['servers'].keys())}"
            )

        return self._config['servers'][actual_server]

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
        return list(self._config['servers'].keys())

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
            name for name, config in self._config['servers'].items()
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
            name for name, config in self._config['servers'].items()
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
_config_instance: Optional[ServerConfig] = None


def get_config(config_file: Optional[str] = None) -> ServerConfig:
    """
    Get the global server configuration instance (lazy-loaded).

    Parameters
    ----------
    config_file : str or Path, optional
        Path to YAML config file. Only used on first call.

    Returns
    -------
    ServerConfig
        Global configuration instance
    """
    global _config_instance
    if _config_instance is None:
        _config_instance = ServerConfig(config_file)
    return _config_instance


def reload_config(config_file: Optional[str] = None) -> ServerConfig:
    """
    Reload the server configuration (useful for testing or config updates).

    Parameters
    ----------
    config_file : str or Path, optional
        Path to YAML config file

    Returns
    -------
    ServerConfig
        Reloaded configuration instance
    """
    global _config_instance
    _config_instance = ServerConfig(config_file)
    return _config_instance

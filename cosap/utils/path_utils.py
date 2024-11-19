
import os

class PathUtils:
    @staticmethod
    def join_paths(path: str, *paths) -> str:
        """Joins and normalizes paths according to os standards"""
        return os.path.normpath(os.path.join(path, *paths))

    @staticmethod
    def convert_to_absolute_path(path: str) -> str:
        """Converts path to absolute path"""
        return os.path.abspath(path)

    @staticmethod
    def is_valid_path(path: str) -> bool:
        """Returns True if path exists and is not empty."""
        return os.path.exists(path) or os.path.isabs(path)

    @staticmethod
    def get_commonpath_from_config(config: dict) -> str:
        """Returns the commonpath of paths that are in the config."""
        paths = []
        for key, value in config.items():
            if isinstance(value, str):
                if PathUtils.is_valid_path(value):
                    paths.append(value)
            elif isinstance(value, list):
                for i in value:
                    if PathUtils.is_valid_path(i):
                        paths.append(i)
            elif isinstance(value, dict):
                for v in value.values():
                    if PathUtils.is_valid_path(v):
                        paths.append(v)
            else:
                raise ValueError("Config value is not a string or list of strings.")
        return os.path.commonpath(paths)
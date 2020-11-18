import os


def join_paths(path: str, *paths) -> str:
    """Joins and normalizes paths according to os standards"""
    return os.path.normpath(os.path.join(path, *paths))

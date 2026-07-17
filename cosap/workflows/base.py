from typing import Protocol

from ..pipeline_builder import Pipeline


class Workflow(Protocol):
    """Common contract for preconstructed COSAP workflows."""

    config: dict

    def build(self) -> Pipeline:
        ...

    def run_pipeline(self) -> str:
        ...

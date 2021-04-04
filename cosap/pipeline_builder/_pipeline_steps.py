from uuid import uuid4


class _PipelineStep:
    def _get_name(self) -> str:
        return uuid4().hex[: 4].upper()

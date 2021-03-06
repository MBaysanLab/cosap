from ._indexer import BamIndexer, _Indexer


class IndexerFactory:
    def create(cls) -> _Indexer:
        return BamIndexer

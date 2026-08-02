# Source Builds

COSAP depends on `cosap-formats` for generic VCF and annotated CNV VCF parsing.
Released packages can resolve `cosap-formats>=0.1,<0.2` normally. Development and
source builds should install both repositories from sibling checkouts and do not
require `cosap-formats` to have been published to PyPI.

Expected layout:

```text
workspace/
|-- cosap/
`-- cosap_formats/
```

Install both editable packages from `cosap/`:

```bash
python -m pip install -e ../cosap_formats -e ".[dev]"
```

The Dockerfile intentionally requires a BuildKit named context so the image is built
from the same local `cosap_formats` source:

```bash
docker buildx build \
  --build-context cosap_formats_source=../cosap_formats \
  -t cosap:local .
```

Callers building through Compose or CI must provide the same required context under
the name `cosap_formats_source`. This avoids silently falling back to a potentially
different package from PyPI.

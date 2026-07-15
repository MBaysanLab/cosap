# Mitochondrial nf-core Fixture

These integration-test fixtures are downloaded on demand from the public
`scorreard/Variant_catalogue_pipeline` test data used by
`nf-core/variantcatalogue`'s `test` profile.

Source repository: https://github.com/scorreard/Variant_catalogue_pipeline

Pinned source commit: `c3d772f51d883a2c858c3a8c017e08774383da81`

License: Apache-2.0

The COSAP repository does not vendor the binary FASTQ/reference files. The gated
integration test downloads them into a temporary directory, validates file sizes,
stages a temporary COSAP hg38 resource bundle, and runs the mitochondrial
workflow against that bundle.

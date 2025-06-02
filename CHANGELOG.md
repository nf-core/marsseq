# nf-core/marsseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0 - [02-06-2025]

### Bug Fixes

- linting issues - ([bedc4a8](https://github.com/nf-core/marsseq/commit/bedc4a887ed804a0f2a5168a4f59090b083d40bd))
- test config - ([0b9410b](https://github.com/nf-core/marsseq/commit/0b9410bcb19ac0e3630ebd500e5891851e0c2c2b))
- check that provided xlsx files are not empty after subset - ([e88075b](https://github.com/nf-core/marsseq/commit/e88075be0f56d9791229bfe590dba594a10122ce))
- missing openpyxl in conda specification - ([41d617b](https://github.com/nf-core/marsseq/commit/41d617b3aa0deece111487ea29d208489662948c))
- missing openpyxl conda dependency - ([d47f7f0](https://github.com/nf-core/marsseq/commit/d47f7f02a9b3367d870bd9b460df6811b61e95aa))
- replace ubuntu image with conda packages - ([b4ef5a5](https://github.com/nf-core/marsseq/commit/b4ef5a520ff401a75451242df34e02aebccd3068))
- perl not anymore in bioconda channel - ([7c8cf13](https://github.com/nf-core/marsseq/commit/7c8cf1320f0ff01710f48d1f08000c71351f1ab1))
- move uncompressed genome and annotation to reference folder - ([4dea965](https://github.com/nf-core/marsseq/commit/4dea96557937785970ca12755a212817113e0753))
- small and full AWS test - ([b72d716](https://github.com/nf-core/marsseq/commit/b72d716416f2c37499c708f2f9f4d1e415e0dfd5))
- linting - ([4f6feee](https://github.com/nf-core/marsseq/commit/4f6feeefe85bfd5547bef81b38e052f6cf2f79e1))

### Documentation

- linting for CITATIONS - ([98b5159](https://github.com/nf-core/marsseq/commit/98b5159c997ad964166c95de01c0673adc964695))
- updated CHANGELOG - ([5a4e7f2](https://github.com/nf-core/marsseq/commit/5a4e7f2cdd5757cd1094a5534cddd702c173449e))
- updated workflow and outputs - ([cddbe3b](https://github.com/nf-core/marsseq/commit/cddbe3b44b7fd0aa189b78eb3f79379ba74091b0))

### Features

- merge with template and full refactor - ([7057892](https://github.com/nf-core/marsseq/commit/7057892182f2eb1e9b846e61e2fc30ba41d8726b))
- support for compressed genome references - ([13b7cb6](https://github.com/nf-core/marsseq/commit/13b7cb67f5e7ed0f184411ba9a42e9f2e8f46b13))

### Miscellaneous Chores

- merge pull request #22 from nf-core/nf-core-template-merge-3.0.2 - ([64a3ee5](https://github.com/nf-core/marsseq/commit/64a3ee5a7012a9f2fe7f4ee2bc116c92de24ce66))

### Wip

- merge template 3.0.2 to dev - ([0b72b26](https://github.com/nf-core/marsseq/commit/0b72b2639d7f480cbf460533750343d5ad54d531))

## v1.0.3 - [11-07-2023]

### `Fixed`

- AWS testing

## v1.0.2 - [11-07-2023]

### `Dependencies`

- [#6](https://github.com/nf-core/marsseq/pull/6) sync with template 2.8 -> 2.9
- [#4](https://github.com/nf-core/marsseq/issues/4) - Fix AWS testing with s3 bucket

## v1.0.1 - [26-06-2023]

### `Fixed`

- [#4](https://github.com/nf-core/marsseq/issues/4) - Fix AWS testing with s3 bucket
- missing zenodo in `lib/`
- added credit and provided short MARS-seq description

## v1.0.0 - [21-06-2023]

### `Dependencies`

- Bump minimal Nextflow version to 23.04.0
- sync with template 2.8

## v1.0dev - [date]

Initial release of nf-core/marsseq, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`

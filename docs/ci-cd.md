# CI/CD

## Workflows

- `CI` runs Ruff, the fast test suite, and a non-publishing container build for pull requests and pushes to `main`.
- `Publish Docker Image` publishes `ghcr.io/marifergun/cosap-base` on every push to `main` and for release tags matching `v*`.
- After publishing, it dispatches the immutable COSAP revision, formats revision, and image digest to `marifergun/cosap_celery` so a compatible compute worker is built automatically.
- `Integration Tests` runs the mitochondrial fixture every Monday and exposes the NA12878 benchmark as a manual action.

## GitHub Container Registry Setup

Configure these GitHub Actions secrets:

- `COSAP_BUILD_TOKEN`: read-only access to `marifergun/cosap_formats`.
- `COSAP_GHCR_TOKEN`: `marifergun` package token with `write:packages` permission.
- `COSAP_DISPATCH_TOKEN`: fine-grained token with Contents write access to `marifergun/cosap_celery`, used only to create repository dispatch events.

Every `main` publication produces an immutable `sha-<commit>` tag and updates the moving `main` tag. Release tags matching `v*` produce semver tags. The workflow never publishes `latest`; deployments should ultimately pin the published manifest digest.

## Repository Settings

Protect `main` with these rules:

- Require a pull request before merging.
- Require one approval.
- Require the `Lint and fast tests` and `Build container image` checks.
- Require branches to be up to date before merging.
- Block force pushes and branch deletion.

Use GitHub Environments for deployment credentials and require approval for production deployments.

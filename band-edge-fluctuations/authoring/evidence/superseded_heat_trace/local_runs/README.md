# Local oracle / nop runs (2026-10-08)

The images were built from the shipped `environment/Dockerfile` and `tests/Dockerfile` with `docker build --network host`. In the authoring sandbox, outbound HTTPS from containers goes through a TLS-intercepting proxy. For these local runs only, the two `python:*-slim-bookworm` base images were re-tagged with that proxy's CA added. The shipped Dockerfiles are unchanged. `harbor run` itself could not build images in this sandbox because compose builds have no network egress. Each step below therefore reproduces what harbor does.

* **oracle**: `docker run --network none --cpus 2 -m 4g -v solution:/solution:ro -v out:/app/output htf-env bash /solution/solve.sh`, then `bash /tests/test.sh` in the verifier image with `out` mounted read-only and `--network none`.
* **nop**: the same verifier run on an empty output directory.

| run | solve wall time | reward |
|---|---|---|
| oracle 1 | 62 s | 1 (5/5 tests passed) |
| oracle 2 | 61 s | 1 (5/5 tests passed) |
| nop | n/a | 0 (`answer.json` missing, so the fixture errors) |

The answers from both oracle runs are stored next to this file. They are bit-identical, because the solver is deterministic.

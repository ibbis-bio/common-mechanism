# Containerfile for the commec DNA sequence screening tool.
#
# Works with both `podman build` and `docker build` (Docker reads this file with
# `-f Containerfile`).
#
# Reference databases are NOT included in this image. They are downloaded at
# runtime by `commec setup` from https://databases.commec.io and are tens of
# gigabytes; bind-mount a host directory instead. See README.md.
#
# Build targets:
#   runtime  (default) the distributable image: commec + its binary dependencies
#   test                runtime + pytest + the source tree, for CI
#
# Architecture: linux/amd64 only. bioconda does not publish linux-aarch64 builds
# for `blast` (and `infernal` coverage is incomplete), so an arm64 image cannot
# be solved today. Apple Silicon / arm64 hosts run this image under emulation.
# TODO(arm64): revisit once bioconda ships linux-aarch64 for blast and infernal
# (https://bioconda.github.io/recipes/blast/README.html,
#  https://bioconda.github.io/recipes/infernal/README.html), then switch the
# CI workflow to buildx with platforms: linux/amd64,linux/arm64.

# ---------------------------------------------------------------------------
# Stage 1: build the conda environment.
#
# The environment is created at its FINAL path (/opt/commec) rather than being
# relocated later: conda packages bake the prefix into their contents (notably
# OpenSSL's CA bundle path), so the prefix must be identical in both stages.
# ---------------------------------------------------------------------------
FROM --platform=linux/amd64 mambaorg/micromamba:2.3.2-debian12-slim AS build

ENV COMMEC_PREFIX=/opt/commec

USER root

COPY environment.yaml /tmp/environment.yaml

# environment.yaml is the dev environment. Everything from its
# "# Development dependencies" marker onwards (pytest, matplotlib, pip, `-e .`)
# is stripped here so the runtime image carries only what screening needs.
# Keeping one dependency list avoids drift between the conda env and the image;
# the marker comment in environment.yaml is load-bearing for this build.
# The `name:` key is dropped too, since the env is created by prefix (-p).
# `pip` is added back because it installs the commec package itself below.
RUN set -eux; \
    sed -e '/# Development dependencies/,$d' -e '/^name:/d' \
        /tmp/environment.yaml > /tmp/environment.runtime.yaml; \
    printf '  - pip\n' >> /tmp/environment.runtime.yaml; \
    cat /tmp/environment.runtime.yaml

RUN micromamba create -y -p "${COMMEC_PREFIX}" -f /tmp/environment.runtime.yaml \
    && micromamba clean --all --yes

COPY . /src

# --no-deps is required, not merely an optimisation: pyproject.toml declares no
# `dependencies`, so every Python dependency comes from the conda environment
# above. Letting pip resolve would either no-op or fight the conda env.
RUN "${COMMEC_PREFIX}/bin/pip" install --no-deps --no-cache-dir /src

# Trim build-time-only artefacts. Static archives, headers, man/info pages and
# bytecode caches are not needed to run commec. Deliberately preserved: ssl/
# (CA bundle), share/ncbi-data and the whole of bin/.
RUN set -eux; \
    rm -rf "${MAMBA_ROOT_PREFIX:-/opt/conda}/pkgs" /root/.cache /tmp/*.yaml /src; \
    find "${COMMEC_PREFIX}" -depth -type d -name __pycache__ -exec rm -rf {} +; \
    find "${COMMEC_PREFIX}" -type f -name '*.pyc' -delete; \
    find "${COMMEC_PREFIX}" -type f -name '*.a' -delete; \
    rm -rf "${COMMEC_PREFIX}/include" \
           "${COMMEC_PREFIX}/share/doc" \
           "${COMMEC_PREFIX}/share/man" \
           "${COMMEC_PREFIX}/share/info" \
           "${COMMEC_PREFIX}/share/terminfo"

# ---------------------------------------------------------------------------
# Stage 2: runtime image.
# ---------------------------------------------------------------------------
FROM --platform=linux/amd64 debian:bookworm-slim AS runtime

ARG COMMEC_VERSION=dev

LABEL org.opencontainers.image.title="commec" \
      org.opencontainers.image.description="Free, open-source, globally available tool for DNA sequence screening. Reference databases are not included; fetch them with 'commec setup'." \
      org.opencontainers.image.source="https://github.com/ibbis-bio/common-mechanism" \
      org.opencontainers.image.url="https://ibbis.bio/common-mechanism" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.version="${COMMEC_VERSION}" \
      org.opencontainers.image.vendor="International Biosecurity and Biosafety Initiative for Science (IBBIS)" \
      io.commec.database-compatibility=">=1.0,<2.0"

# The conda environment ships its own CA bundle; system certificates are added
# so anything falling back to the OS trust store also works.
RUN set -eux; \
    apt-get update; \
    apt-get install -y --no-install-recommends ca-certificates; \
    rm -rf /var/lib/apt/lists/*

COPY --from=build /opt/commec /opt/commec

ENV PATH=/opt/commec/bin:$PATH \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    PYTHONDONTWRITEBYTECODE=1

# Non-root by default. /databases is the conventional mount point for reference
# data and /data is the working directory where input is read and output written.
RUN set -eux; \
    useradd --create-home --uid 1000 --shell /bin/bash commec; \
    mkdir -p /data /databases; \
    chown commec:commec /data /databases

WORKDIR /data
USER commec

ENTRYPOINT ["commec"]
CMD ["--help"]

# ---------------------------------------------------------------------------
# Stage 3: test image (CI only).
#
# The tests need a source checkout, not just the installed package:
# commec/tests/test_data and commec/tests/test_dbs are not declared in
# package-data, so they are absent from site-packages. test_dbs.py actually
# executes blastn/blastx/hmmscan/cmscan against those fixtures and writes into
# the working directory, hence the chown.
# ---------------------------------------------------------------------------
FROM runtime AS test

USER root
RUN /opt/commec/bin/pip install --no-cache-dir pytest

COPY --chown=commec:commec . /src
WORKDIR /src
USER commec

ENTRYPOINT []
CMD ["pytest", "-vv"]

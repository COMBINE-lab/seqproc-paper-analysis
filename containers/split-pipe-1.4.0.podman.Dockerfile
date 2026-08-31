# Rootless-Podman-compatible adaptation of ParseNIP commit
# f1ba1291aee95b4c6fa85bf0ce6678ffd078f6a1, dockerfiles/split-pipe.Dockerfile.
#
# The split-pipe archive is not redistributed. It is fetched from the same
# public FTP path used by ParseNIP and verified before installation. Use of
# split-pipe is governed by the Parse Biosciences Software License Agreement
# bundled in that archive.

FROM ubuntu:noble-20241011

LABEL org.opencontainers.image.title="split-pipe 1.4.0 calibration environment"
LABEL org.opencontainers.image.source="https://github.com/StevenWingett/ParseNIP"
LABEL org.opencontainers.image.version="1.4.0"
LABEL seqproc.parsenip.commit="f1ba1291aee95b4c6fa85bf0ce6678ffd078f6a1"

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# This account has a single-ID rootless user namespace. Keeping apt's download
# method as container root avoids its otherwise failing transition to `_apt`.
RUN printf 'APT::Sandbox::User "root";\n' > /etc/apt/apt.conf.d/00sandbox-root \
    && apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        ca-certificates curl unzip \
    && rm -rf /var/lib/apt/lists/*

ARG PBP_ARCHIVE="PBP.1.4.0.zip"
ARG PBP_SHA256="6a8d54452f28489585cc7e3923825278aa1deb0b96fab9ebb082d6b594a122d8"
ARG PBP_URL="ftp://ftp.mrc-lmb.cam.ac.uk/pub/swingett/software_download/PBP.1.4.0.zip"
ARG MINICONDA_INSTALLER="Miniconda3-py311_23.5.2-0-Linux-x86_64.sh"
ARG MINICONDA_SHA256="634d76df5e489c44ade4085552b97bebc786d49245ed1a830022b0b406de5817"
ARG MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-py311_23.5.2-0-Linux-x86_64.sh"

WORKDIR /opt/myfiles
RUN curl --fail --location --retry 5 --output "${PBP_ARCHIVE}" "${PBP_URL}" \
    && echo "${PBP_SHA256}  ${PBP_ARCHIVE}" | sha256sum --check - \
    && unzip "${PBP_ARCHIVE}" \
    && rm "${PBP_ARCHIVE}" \
    && curl --fail --location --retry 5 --output "${MINICONDA_INSTALLER}" "${MINICONDA_URL}" \
    && echo "${MINICONDA_SHA256}  ${MINICONDA_INSTALLER}" | sha256sum --check - \
    && bash "${MINICONDA_INSTALLER}" -b -p /opt/miniconda3 \
    && rm "${MINICONDA_INSTALLER}"

ENV CONDA_DIR=/opt/miniconda3
ENV PATH=/opt/miniconda3/bin:${PATH}

RUN conda create --name spipe python=3.10 --yes

WORKDIR /opt/myfiles/ParseBiosciences-Pipeline.1.4.0
RUN eval "$(conda shell.bash hook)" \
    && conda activate spipe \
    && bash ./install_dependencies_conda.sh --install --yes \
    && pip install --no-cache-dir ./ \
    && split-pipe --version \
    && python --version \
    && conda list --explicit > /opt/split-pipe-1.4.0-conda-explicit.txt \
    && python -m pip freeze > /opt/split-pipe-1.4.0-pip-freeze.txt

ENV PATH=/opt/miniconda3/envs/spipe/bin:/opt/miniconda3/bin:${PATH}

ENTRYPOINT ["split-pipe"]
CMD ["--help"]

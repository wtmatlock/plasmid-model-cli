# syntax=docker/dockerfile:1.7

ARG RUST_VERSION=1.85
ARG CMDSTAN_VERSION=2.36.0
ARG STAN_MODEL=stan/model.stan

FROM debian:bookworm AS cmdstan-builder
ARG CMDSTAN_VERSION
ARG STAN_MODEL

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential g++ make gfortran wget ca-certificates tar python3 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt
RUN wget -q https://github.com/stan-dev/cmdstan/releases/download/v${CMDSTAN_VERSION}/cmdstan-${CMDSTAN_VERSION}.tar.gz \
    && tar -xzf cmdstan-${CMDSTAN_VERSION}.tar.gz \
    && rm cmdstan-${CMDSTAN_VERSION}.tar.gz

WORKDIR /opt/cmdstan-${CMDSTAN_VERSION}
RUN make build -j"$(nproc)"

WORKDIR /work
COPY . .

RUN test -f "${STAN_MODEL}" && /opt/cmdstan-${CMDSTAN_VERSION}/bin/stanc "${STAN_MODEL}" >/dev/null
RUN test -f "${STAN_MODEL}" \
    && rm -f "/work/${STAN_MODEL%.stan}" \
    && make -C /opt/cmdstan-${CMDSTAN_VERSION} "/work/${STAN_MODEL%.stan}" -j"$(nproc)" \
    && test -x "/work/${STAN_MODEL%.stan}"

FROM rust:${RUST_VERSION}-bookworm AS rust-builder

WORKDIR /work
COPY Cargo.toml Cargo.lock ./
COPY src ./src

RUN cargo build --release

FROM debian:bookworm-slim AS runtime
ARG CMDSTAN_VERSION
ARG STAN_MODEL

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates libstdc++6 libgomp1 libtbb12 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
ENV CMDSTAN=/opt/cmdstan
ENV STAN_MODEL_PATH=/app/${STAN_MODEL%.stan}
ENV LD_LIBRARY_PATH=/opt/cmdstan/stan/lib/stan_math/lib/tbb:${LD_LIBRARY_PATH}

COPY --from=rust-builder /work/target/release/plasmid-model-cli /app/plasmid-model-cli
COPY --from=cmdstan-builder /opt/cmdstan-${CMDSTAN_VERSION} /opt/cmdstan
COPY --from=cmdstan-builder /work/${STAN_MODEL%.stan} /app/${STAN_MODEL%.stan}
COPY --from=cmdstan-builder /work/${STAN_MODEL} /app/${STAN_MODEL}

ENTRYPOINT ["/app/plasmid-model-cli"]
CMD ["--help"]
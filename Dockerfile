FROM rust:1.85.1 AS builder

RUN apt-get update && \
    apt-get install -y cmake

FROM builder AS build

COPY . /opt/filter-clipped/
WORKDIR /opt/filter-clipped
RUN cargo install --path .

FROM debian:buster-slim as exec
RUN apt-get update && rm -rf /var/lib/apt/lists/*
COPY --from=build /usr/local/cargo/bin/filter-clipped /usr/local/bin/filter-clipped
ENV RUST_LOG=info
ENTRYPOINT ["/usr/local/bin/filter-clipped"]

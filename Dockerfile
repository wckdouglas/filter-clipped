FROM rust:1.62.1 AS builder

RUN apt-get update && \
    apt-get install -y cmake

FROM builder AS build

COPY . /opt/filter-clipped/
WORKDIR /opt/filter-clipped
RUN cargo install --path .

FROM build AS final
ENV RUST_LOG=info
ENTRYPOINT ["/usr/local/cargo/bin/filter-clipped"]

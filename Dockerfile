FROM rust:1.62.1 as builder

COPY . /opt/filter-clipped/
WORKDIR /opt/filter-clipped
RUN cargo install --path .

FROM builder as final
ENV RUST_LOG=info
ENTRYPOINT ["/usr/local/cargo/bin/filter-clipped"]

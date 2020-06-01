FROM ubuntu:20.04 AS builder

RUN apt update && apt install -y --no-install-recommends \
  ca-certificates \
  libgsl-dev \
  libhts-dev \
  libbigwig-dev \
  libcurl4-openssl-dev \
  gcc \
  python \
  make

WORKDIR /WiggleTools

COPY . .

RUN make LIBS='-lwiggletools -lBigWig -lcurl -lhts -lgsl  -lgslcblas -lz -lpthread -lm -llzma' \
  && make test

FROM ubuntu:20.04

RUN apt update && apt install -y --no-install-recommends \
  ca-certificates \
  libbigwig0 \
  libcurl4 \
  libgsl23 \
  libhts3

COPY --from=builder /WiggleTools/bin/wiggletools /usr/local/bin

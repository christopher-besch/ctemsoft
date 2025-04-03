FROM debian AS builder

RUN apt-get update && \
    apt-get -y install make gfortran liblapack3 libblas3 liblapack-dev libblas-dev

RUN mkdir -p /code/exe
COPY ./src/ /code/src

WORKDIR /code/src
RUN make all

FROM debian

RUN apt-get update && \
    apt-get -y install libgfortran5 python3 python3-numpy ghostscript ffmpeg

COPY --from=builder /code/exe/lens /usr/bin/lens
COPY ./python_src/ /python_src
WORKDIR /python_src
ENTRYPOINT ["/bin/python3", "/python_src/main.py"]

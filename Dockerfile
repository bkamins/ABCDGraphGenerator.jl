FROM julia:1.7.3-alpine

WORKDIR /app

ADD utils/ .

RUN julia install.jl





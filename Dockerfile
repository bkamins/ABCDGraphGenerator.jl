FROM julia:1.6-alpine

WORKDIR /app

ADD utils/ .

RUN julia install.jl





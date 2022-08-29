FROM julia:1.7.3

WORKDIR /app

ADD utils/ .

RUN julia install.jl





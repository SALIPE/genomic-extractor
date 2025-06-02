FROM julia:1.11.5
WORKDIR /app

COPY ./GREAC/ /app


RUN julia -e 'using Pkg; Pkg.activate(".");Pkg.instantiate()'

CMD [ "julia" ]

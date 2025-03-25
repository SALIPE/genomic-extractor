FROM julia:1.11.3
WORKDIR /app

ADD /modules /app/modules/
ADD Main.jl /app/
ADD *.toml /app/


RUN julia -e 'using Pkg; Pkg.add("AbstractFFTs")'
RUN julia -e 'using Pkg; Pkg.add("ArgParse")'
RUN julia -e 'using Pkg; Pkg.add("BenchmarkTools")'
RUN julia -e 'using Pkg; Pkg.add("BioSequences")'
RUN julia -e 'using Pkg; Pkg.add("DSP")'
RUN julia -e 'using Pkg; Pkg.add("Distributed")'
RUN julia -e 'using Pkg; Pkg.add("FASTX")'
RUN julia -e 'using Pkg; Pkg.add("FLoops")'
RUN julia -e 'using Pkg; Pkg.add("Normalization")'
RUN julia -e 'using Pkg; Pkg.add("Pickle")'
RUN julia -e 'using Pkg; Pkg.add("Polynomials")'
# RUN julia -e 'using Pkg; Pkg.add("Plots")'
RUN julia -e 'using Pkg; Pkg.add("StatsBase")'

# CMD export JULIA_DEPOT_PATH="/root/.julia/" &&
CMD julia --project 
# Main.jl -w 0.004 --classify -o ./teste

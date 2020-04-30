# Start from Ubuntu (e.g., bionic - 18.04.4)
FROM ubuntu

# Set default language to UTF-8 US English
RUN apt-get update -y \
    && apt-get install --no-install-recommends -y \
        locales \
    && rm -rf /var/lib/apt/lists/* \
    && locale-gen en_US.UTF-8
ENV LC_ALL=en_US.UTF-8

# Install Python (e.g., 3.6.9), pip, and setuptools
RUN apt-get update -y \
    && apt-get install --no-install-recommends -y \
        python3 \
        python3-pip \
    && pip3 install -U pip setuptools \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*
    
# Install Julia v1.0.5
RUN apt update
RUN apt upgrade -y
RUN apt install curl -y
RUN curl https://julialang-s3.julialang.org/bin/linux/x64/1.0/julia-1.0.5-linux-x86_64.tar.gz | tar -xz -C /opt/
RUN ln -s /opt/julia-1.0.5/bin/julia /usr/local/bin/julia

# Install Julia packages
RUN julia -e 'using Pkg; Pkg.add(PackageSpec(name="JuMP", version="0.21.1"))'
RUN julia -e 'using Pkg; Pkg.add(PackageSpec(name="Ipopt", version="0.6.1"))'
RUN julia -e 'using Pkg; Pkg.add(PackageSpec(name="CSV", version="0.6.1"))'
RUN julia -e 'using Pkg; Pkg.add(PackageSpec(name="PyCall", version="1.91.4"))'

# Copy DisFit and associated files to Docker image
#COPY DisFit /DisFit/DisFit
#COPY docs /DisFit/docs
#COPY examples /DisFit/examples
#COPY tests /DisFit/tests
COPY . /DisFit 

# Install DisFit
RUN pip3 install /DisFit
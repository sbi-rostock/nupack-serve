<div id="top"></div>

[![Docker Repository on Quay](https://quay.io/repository/bagnacan/nupack-serve/status "Docker Repository on Quay")](https://quay.io/repository/bagnacan/nupack-serve)

# Nupack-serve

Nupack-serve is a wrapper for [NUPACK](http://www.nupack.org/): a software
suite for the analysis and design of nucleic acid structures, devices, and
systems ().  
This wrapper modifies the default behavior of some of NUPACK's functions in
order to:
- Incoporate NUPACK within computational pipelines
- Make results machine readabile (JSON)
- Restrict computation ot CPU-bound operations
- Provide provenance medatada
- Enhance reproducibility

Nupack-serve was initially created to port NUPACK's *mfe*, *complexes* and
*concentrations* within the [Triplexer](https://github.com/sbi-rostock/triplexer).
However, the C functions written to remove I/O bound operation, provide
provenance medatada and JSON output, can be adopted by the remaining NUPACK
utilities.

- [Installation requirements](#installation-requirements)
- [Run nupack-serve](#run-nupack-serve)


## Installation requirements

The only requirement is [Docker](https://www.docker.com/), which can be
installed in different ways depending on the underlying operative system:
- Unix users should follow the [Docker installation for Linux](https://docs.docker.com/compose/install/#install-compose-on-linux-systems#install-compose-on-linux-systems),
and install both Docker and Docker compose
- MacOS 10.13+ users should follow the [Docker installation for Mac](https://docs.docker.com/docker-for-mac/install/)
- Windows 10+ users, should follow the [Docker installation for Windows](https://docs.docker.com/docker-for-windows/install/)
- For legacy systems, users can rely on the [Docker Toolbox](https://docs.docker.com/toolbox/overview/).
<p align="right"><a href="#top">&#x25B2; back to top</a></p>



## Run nupack-serve

To run nupack-serve, you need to run its docker container. Type:
```
docker run -it quay.io/bagnacan/nupack-serve
```

You can now access the wrapped NUPACK functions.
<p align="right"><a href="#top">&#x25B2; back to top</a></p>

#!/usr/bin/env python3

# This module implements the nupack-serve app, which handles calls to:
# - serve-mfe
# - serve-complexes
# - serve-concentrations

from starlette.applications import Starlette
from starlette.responses import JSONResponse
from starlette.requests import Request
from common import *
import json
import uvicorn
import serve_mfe
import serve_complexes
import serve_concentrations


USAGE    = "usage"
HOMEPAGE = "homepage"
LICENSE  = "license"
STATUS   = "status"
RESULT   = "result"

# values
NUPACK_LICENSE_LOCATION = "/tmp/nupack3.2.2/LICENSE"
USAGE_DESCRIPTION = "Nupack-serve is a wrapper for NUPACK: a software suite for the analysis and design of nucleic acid structures, devices, and systems (Wolfe et al. 2017)."
HOMEPAGE_URL = "https://github.com/sbi-rostock/nupack-serve"



# return NUPACK's license
#
def nupack_license():
    """
    Returns NUPACK's license.
    """

    terms = None

    with open(NUPACK_LICENSE_LOCATION) as nupack_license:
        terms = nupack_license.read()

    return terms



# return nupack-serve response:
# - NUPACK license
# - Subprocess status code
# - Result
#
def serve_result(status, result):
    """
    Assembles the nupack-serve response and returns it.
    """

    return JSONResponse({
        LICENSE: nupack_license(),
        STATUS: status,
        RESULT: result
    })



#
# nupack-serve
#
app = Starlette(debug=True)



#
# about
#
@app.route("/", methods=["GET"])
def about(request):
    """
    Returns the about page of nupack-serve.
    """

    return JSONResponse({
        USAGE: USAGE_DESCRIPTION,
        HOMEPAGE: HOMEPAGE_URL,
        LICENSE: nupack_license()
    })



#
# mfe
#
@app.route("/mfe", methods=["GET"])
async def mfe(request):
    """
    Returns nupack-serve mfe.
    """

    parameters = await request.json()

    status, result = serve_mfe.mfe(parameters)
    return serve_result(status, result)



#
# mfe example:
#
# Compute the MFE of the complex formed by aligning Human target gene E2F1's
# binding site sequence and the putatively cooperating miRNA pair hsa-miR-205
# and hsa-miR-342-3p
#
@app.route("/example/mfe", methods=["GET"])
def example_mfe(request):
    """
    A nupack-serve mfe example:
    Compute the MFE of the complex formed by aligning Human target gene E2F1's
    binding site sequence and the putatively cooperating miRNA pair
    hsa-miR-205 and hsa-miR-342-3p.
    """

    parameters = {
        SEQ_NUM: "3",
        SEQ_TARGET: "ccgggggugaaugugugugagcaugugugugugcauguaccggggaaugaaggu",
        SEQ_MIR1: "uccuucauuccaccggagucug",
        SEQ_MIR2: "ucucacacagaaaucgcacccgu",
        PERMUTATIONS: ["1 2 3"]
    }
    status, result = serve_mfe.mfe(parameters)
    return serve_result(status, result)



#
# complexes
#
@app.route("/complexes", methods=["GET"])
async def complexes(request):
    """
    Returns nupack-serve complexes.
    """

    parameters = await request.json()

    status, result = serve_complexes.complexes(parameters)
    return serve_result(status, result)



#
# complexes example:
#
# Compute the partition function and equilibrium base-pairing properties for
# each complex formed by aligning Human target gene E2F1's  binding site
# sequence and the putatively cooperating miRNA pair hsa-miR-205 and
# hsa-miR-342-3p
#
@app.route("/example/complexes", methods=["GET"])
def example_complexes(request):
    """
    A nupack-serve complexes example:
    Compute the partition function and equilibrium base-pairing properties for
    each complex formed by aligning Human target gene E2F1's  binding site
    sequence and the putatively cooperating miRNA pair hsa-miR-205 and
    hsa-miR-342-3p.
    """

    parameters = {
        SEQ_NUM: "3",
        SEQ_TARGET: "ccgggggugaaugugugugagcaugugugugugcauguaccggggaaugaaggu",
        SEQ_MIR1: "uccuucauuccaccggagucug",
        SEQ_MIR2: "ucucacacagaaaucgcacccgu",
        MAX_COMPLEX_SIZE: "1",
        PERMUTATIONS: ["1 2 3", "1 2", "1 3"]
    }
    status, result = serve_complexes.complexes(parameters)
    return serve_result(status, result)



#
# concentrations
#
@app.route("/concentrations", methods=["POST"])
async def concentrations(request):
    """
    Returns nupack-serve concentrations.
    """

    parameters = await request.json()

    status, result = serve_concentrations.concentrations(parameters)
    return serve_result(status, result)



#
# nupac-serve concentrations example:
#
# Compute the equilibrium concentration for each complex that is formed by
# aligning Human target gene E2F1's with the putatively cooperating miRNA
# pair hsa-miR-205 and hsa-miR-342-3p
#
@app.route("/example/concentrations", methods=["GET"])
def example_concentrations(request):
    """
    A nupack-serve concentrations example:
    Compute the equilibrium concentration for each complex that is formed by
    aligning Human target gene E2F1's with the putatively cooperating miRNA
    pair hsa-miR-205 and hsa-miR-342-3p.
    """

    parameters = {
        NUM_COMPLEXES: "6",
        LIST_CONCENTRATIONS: ["1e-7", "1e-7", "1e-7"],
        TEMP: "37.0",
        OCX: [
            "1,1,1,0,0,-7.92078773e+00",
            "2,1,0,1,0,-9.79502400e+00",
            "3,1,0,0,1,-9.79502400e+00",
            "4,1,1,1,0,-4.84277745e+01",
            "5,1,1,0,1,-4.84277745e+01",
            "6,1,1,1,1,-6.36285141e+01"
        ]
    }
    status, result = serve_concentrations.concentrations(parameters)
    return serve_result(status, result)



#
# main
#
if __name__ == "__main__":
    uvicorn.run(app, host='0.0.0.0', port=8000)


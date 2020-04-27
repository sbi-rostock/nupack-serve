#!/usr/bin/env python3

# This module implements the nupack-serve app, which handles calls to:
# - serve-mfe
# - serve-complexes
# - serve-concentrations

from starlette.applications import Starlette
from starlette.responses import JSONResponse
from starlette.routing import Route
from starlette.staticfiles import StaticFiles
from common import *
import json
import uvicorn
import serve_mfe
import serve_complexes
import serve_concentrations


NUPACK_LICENSE_LOCATION = "/tmp/nupack3.2.2/LICENSE"
NUPACK_LICENSE_TERMS = None
LICENSE = "license"
STATUS  = "status"
RESULT  = "result"



# response assembler
def response(status, result):
  r = {}

  r[LICENSE] = NUPACK_LICENSE_TERMS
  r[STATUS]  = status
  r[RESULT]  = result
  return JSONResponse(r)



#
# app
#
app = Starlette(debug=True)



#
# nupac-serve mfe
#
@app.route("/mfe", methods=["POST"])
def mfe(payload):

  status, result = serve_mfe.mfe(payload)

  # send response
  return response(status, result)

#
# nupac-serve mfe example:
#
# Compute the MFE of the complex formed by aligning Human target gene E2F1's
# binding site sequence and the putatively cooperating miRNA pair hsa-miR-205
# and hsa-miR-342-3p
#
@app.route("/example/mfe")
def example_mfe(request):

  payload = {
    SEQ_NUM: "3",
    SEQ_TARGET: "ccgggggugaaugugugugagcaugugugugugcauguaccggggaaugaaggu",
    SEQ_MIR1: "uccuucauuccaccggagucug",
    SEQ_MIR2: "ucucacacagaaaucgcacccgu",
    PERMUTATIONS: ["1 2 3"]
  }

  return mfe(json.dumps(payload))



#
# nupac-serve complexes
#
@app.route("/complexes", methods=["POST"])
def complexes(payload):

  status, result = serve_complexes.complexes(payload)

  # send response
  return response(status, result)

#
# nupac-serve complexes example:
#
# Compute the partition function and equilibrium base-pairing properties for
# each complex formed by aligning Human target gene E2F1's  binding site
# sequence and the putatively cooperating miRNA pair hsa-miR-205 and
# hsa-miR-342-3p
#
@app.route("/example/complexes")
def example_complexes(request):

  payload = {
    SEQ_NUM: "3",
    SEQ_TARGET: "ccgggggugaaugugugugagcaugugugugugcauguaccggggaaugaaggu",
    SEQ_MIR1: "uccuucauuccaccggagucug",
    SEQ_MIR2: "ucucacacagaaaucgcacccgu",
    MAX_COMPLEX_SIZE: "1",
    PERMUTATIONS: ["1 2 3", "1 2", "1 3"]
  }

  return complexes(json.dumps(payload))



#
# nupac-serve concentrations
#
@app.route("/concentrations", methods=["POST"])
def concentrations(payload):

  status, result = serve_concentrations.concentrations(payload)

  # send response
  return response(status, result)

#
# nupac-serve concentrations example:
#
# Compute the equilibrium concentration for each complex that is formed by
# aligning Human target gene E2F1's with the putatively cooperating miRNA
# pair hsa-miR-205 and hsa-miR-342-3p
#
@app.route("/example/concentrations")
def example_concentrations(request):

  payload = {
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

  return concentrations(json.dumps(payload))



#
# main
#
if __name__ == "__main__":
  with open(NUPACK_LICENSE_LOCATION) as nupack_license:
    NUPACK_LICENSE_TERMS = nupack_license.read()

  uvicorn.run(app, host='0.0.0.0', port=8000)

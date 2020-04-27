# nupack-serve app

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



#
# app
#
app = Starlette(debug=True)



#
# nupac-serve mfe
#
@app.route("/mfe/{target}/{mir1}/{mir2}")
def mfe(request):

  global NUPACK_LICENSE_TERMS
  target = request.path_params["target"]
  mir1 = request.path_params["mir1"]
  mir2 = request.path_params["mir2"]
  response = {}

  status, result = serve_mfe.mfe(target, mir1, mir2)

  # send response
  response[LICENSE] = NUPACK_LICENSE_TERMS
  response[STATUS]  = status
  response[RESULT]  = result
  return JSONResponse(response)

#
# nupac-serve mfe example:
#
# Compute the MFE of the complex formed by aligning Human target gene E2F1's
# binding site sequence and the putatively cooperating miRNA pair hsa-miR-205
# and hsa-miR-342-3p
#
@app.route("/example/mfe")
def example_mfe(request):

  global NUPACK_LICENSE_TERMS
  target = "ccgggggugaaugugugugagcaugugugugugcauguaccggggaaugaaggu"
  mir1 = "uccuucauuccaccggagucug"
  mir2 = "ucucacacagaaaucgcacccgu"
  response = {}

  status, result = serve_mfe.mfe(target, mir1, mir2)

  # send response
  response[LICENSE] = NUPACK_LICENSE_TERMS
  response[STATUS]  = status
  response[RESULT]  = result
  return JSONResponse(response)



#
# nupac-serve complexes
#
@app.route("/complexes", methods=["POST"])
def complexes(payload):

  global NUPACK_LICENSE_TERMS
  response = {}

  status, result = serve_complexes.complexes(payload)

  # send response
  response[LICENSE] = NUPACK_LICENSE_TERMS
  response[STATUS]  = status
  response[RESULT]  = result
  return JSONResponse(response)

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

  response = complexes(json.dumps(payload))
  return response



#
# nupac-serve concentrations
#
@app.route("/concentrations", methods=["POST"])
def concentrations(payload):

  global NUPACK_LICENSE_TERMS
  response = {}

  status, result = serve_concentrations.concentrations(payload)

  # send response
  response[LICENSE] = NUPACK_LICENSE_TERMS
  response[STATUS]  = status
  response[RESULT]  = result
  return JSONResponse(response)

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

  response = concentrations(json.dumps(payload))
  return response



#
# main
#
if __name__ == "__main__":
  with open(NUPACK_LICENSE_LOCATION) as nupack_license:
    NUPACK_LICENSE_TERMS = nupack_license.read()

  uvicorn.run(app, host='0.0.0.0', port=8000)

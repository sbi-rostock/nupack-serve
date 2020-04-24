# nupack-serve app

from starlette.applications import Starlette
from starlette.responses import JSONResponse
from starlette.routing import Route
from starlette.staticfiles import StaticFiles
import uvicorn
import serve_mfe


NUPACK_LICENSE_LOCATION = "nupack3.2.2/LICENSE"
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
# main
#
if __name__ == "__main__":
  with open(NUPACK_LICENSE_LOCATION) as nupack_license:
    NUPACK_LICENSE_TERMS = nupack_license.read()

  uvicorn.run(app, host='0.0.0.0', port=8000)

from starlette.applications import Starlette
from starlette.responses import JSONResponse
from starlette.routing import Route
from starlette.staticfiles import StaticFiles
import ast
import io
import pexpect
import json
import uvicorn


NUPACK_LICENSE_LOCATION = "nupack3.2.2/LICENSE"
NUPACK_LICENSE_TERMS = None
LICENSE = "license"
STATUS  = "status"
RESULT  = "result"



# app
app = Starlette(debug=True)



# nupac-serve mfe
@app.route("/mfe/{target}/{mir1}/{mir2}")
def mfe(request):
  target = request.path_params["target"]
  mir1 = request.path_params["mir1"]
  mir2 = request.path_params["mir2"]

  global NUPACK_LICENSE_TERMS
  result = None
  response = {}

  # spawn a subprocess for the provided call
  p = pexpect.spawn("mfe -multi", encoding="ascii")

  # suppress echo
  p.setecho(False)

  # provide number of sequences
  p.expect("Enter number of strands")
  p.sendline("3")

  # provide target sequence
  p.expect("Enter sequence for strand type")
  p.sendline(target)

  # provide mir1 sequence
  p.expect("Enter sequence for strand type")
  p.sendline(mir1)

  # provide mir2 sequence
  p.expect("Enter sequence for strand type")
  p.sendline(mir2)

  # provide order
  p.expect("Enter strand permutation")
  p.sendline("1 2 3")

  # start log the subprocess' result
  log = io.StringIO()
  p.logfile = log
  p.expect("{") # starting of JSON string
  log.seek(0)

  # store the result
  result = log.read()
  result = result.replace("\x00", "").replace("\n", "").replace("\r", "")

  log.close()
  p.close()

  # send response
  response[LICENSE] = NUPACK_LICENSE_TERMS
  response[STATUS]  = p.exitstatus
  response[RESULT]  = json.loads(result)
  return JSONResponse(response)



# routes
#routes = [
#  Route('/mfe/{target}/{mir1}/{mir2}', mfe),
#]

# app
#app = Starlette(debug=True, routes=routes)

if __name__ == "__main__":
  with open(NUPACK_LICENSE_LOCATION) as nupack_license:
    NUPACK_LICENSE_TERMS = nupack_license.read()

  uvicorn.run(app, host='0.0.0.0', port=8000)

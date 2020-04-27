#!/usr/bin/env python3

from common import *
import io
import json
import pexpect


def concentrations(payload):

    payload = json.loads(payload)
    result = None

    # spawn a subprocess
    p = pexpect.spawn("concentrations", encoding="ascii")

    # suppress echo
    p.setecho(False)

    # provide number of complexes
    p.expect("Enter number of complex IDs")
    p.sendline(str(payload[NUM_COMPLEXES]))

    # provide number of concentrations
    p.expect("Enter number of different concentrations")
    p.sendline(str(len(payload[LIST_CONCENTRATIONS])))

    # provide entries concentrations
    for entry in range(len(payload[LIST_CONCENTRATIONS])):
      p.expect("Enter concentration")
      p.sendline(payload[LIST_CONCENTRATIONS][entry])

    # provide temperature
    p.expect("Enter temperature")
    p.sendline(payload[TEMP])

    # provide entries ocx
    for entry in range(len(payload[OCX])):
      p.expect("Enter complex ocx")
      p.sendline(payload[OCX][entry])

    # start log the subprocess' result
    log = io.StringIO()
    p.logfile = log
    p.expect("{") # starting of JSON string
    log.seek(0)

    # store the result
    result = log.read()
    result = result.replace("\x00", "").replace("\n", "").replace("\r", "")
    result = json.loads(result)

    log.close()
    p.close()

    return (p.exitstatus, result)


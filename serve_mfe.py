#!/usr/bin/env python3

import io
import json
import pexpect


def mfe(target, mir1, mir2):

    result = None

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
    result = json.loads(result)

    log.close()
    p.close()

    return (p.exitstatus, result)


#!/usr/bin/env python3

import io
import json
import pexpect


def complexes(target, mir1, mir2):

    result = None

    # spawn a subprocess for the provided call
    p = pexpect.spawn(call, encoding="ascii")

    # suppress echo
    p.setecho(False)

    # interact to provide its input
    p.expect("Enter number of different sequences")
    p.sendline("3")

    # provide target sequence
    p.expect("Enter sequence")
    p.sendline(target)

    # provide mir1 sequence
    p.expect("Enter sequence")
    p.sendline(mir1)

    # provide mir2 sequence
    p.expect("Enter sequence")
    p.sendline(mir2)

    p.expect("Enter max complex size to completely enumerate")
    p.sendline("1")

    # premutation 1
    p.expect("Enter permutation")
    p.sendline("1 2 3")

    # premutation 2
    p.expect("Enter permutation")
    p.sendline("1 2")

    # premutation 2
    p.expect("Enter permutation")
    p.sendline("1 3")

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


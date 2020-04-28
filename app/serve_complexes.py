#!/usr/bin/env python3

# This module handles the execution of the modified NUPACK complexes function.

from common import *
import io
import json
import pexpect


def complexes(parameters):

    result = None

    # spawn a subprocess
    p = pexpect.spawn("complexes", encoding="ascii")

    # suppress echo
    p.setecho(False)

    # provide number of sequences
    p.expect("Enter number of different sequences")
    p.sendline(parameters[SEQ_NUM])

    # provide all sequences
    p.expect("Enter sequence")
    p.sendline(parameters[SEQ_TARGET])

    # provide mir1 sequence
    p.expect("Enter sequence")
    p.sendline(parameters[SEQ_MIR1])

    # provide mir2 sequence
    p.expect("Enter sequence")
    p.sendline(parameters[SEQ_MIR2])

    # provide max complex size
    p.expect("Enter max complex size to completely enumerate")
    p.sendline(parameters[MAX_COMPLEX_SIZE])

    # provide permutations
    for entry in parameters[PERMUTATIONS]:
        p.expect("Enter permutation")
        p.sendline(entry)

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


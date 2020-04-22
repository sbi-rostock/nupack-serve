#!/usr/bin/env python3

LICENSE_PATH = "/tmp/nupack3.2.2/LICENSE"

def print_license():
    with open(LICENSE_PATH, 'r') as license:
        print(license.read(), end="")



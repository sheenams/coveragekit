#!/usr/bin/env python

import sys
import coveragekit.covdb as covdb
import coveragekit.covbam as covbam

from coveragekit.version import __version__

usage = '''
coveragekit - v{}
coveragekit.py <command> [options]

bam     import bam data
db      work with coverage database
'''.format(__version__)

try:
    if sys.argv[1] == "bam":
        covbam.run(sys.argv[2:])
    elif sys.argv[1] == "db":
        covdb.run(sys.argv[2:])
    else:
        print usage
        sys.exit(1)
except IndexError as e:
    print usage
    raise
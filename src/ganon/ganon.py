#!/usr/bin/env python3
import sys, time

from ganon.build_update import build, update
from ganon.classify_report import classify, report
from ganon.config import Config
from ganon.util import print_log

def main(which: str=None, **kwargs):

    # sys.argv by default or get args from call
    cfg = Config(which, **kwargs)

    # validate arguments
    if not cfg.validate(): sys.exit(1)

    tx_total = time.time()
    
    print_log("- - - - - - - - - -")
    print_log("   _  _  _  _  _   ")
    print_log("  (_|(_|| |(_)| |  ")
    print_log("   _|   v. "+ str(cfg.version))
    print_log("- - - - - - - - - -")

    if cfg.which=='build':
        ret=build(cfg)
    elif cfg.which=='update': 
        ret=update(cfg) 
    elif cfg.which=='classify':
        ret=classify(cfg)
    elif cfg.which=='report':
        ret=report(cfg)

    print_log("Total elapsed time: " + str("%.2f" % (time.time() - tx_total)) + " seconds.")
    return ret

if __name__ == '__main__':
    main()

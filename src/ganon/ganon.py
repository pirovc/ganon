#!/usr/bin/env python3
import sys, time

from src.ganon.build_update import build, update
from src.ganon.classify_report import classify, report
from src.ganon.config import Config
from src.ganon.util import print_log

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
        build(cfg)
    elif cfg.which=='update': 
        update(cfg) 
    elif cfg.which=='classify':
        classify(cfg)
    elif cfg.which=='report':
        report(cfg)

    print_log("Total elapsed time: " + str("%.2f" % (time.time() - tx_total)) + " seconds.")

if __name__ == '__main__':
    main()

#!/usr/bin/env python3
import sys, time

from ganon.build_update import build, update
from ganon.classify_report import classify, report
from ganon.config import Config
from ganon.util import print_log

def main(which: str=None, cfg=None, **kwargs):
    # 3 entry points: 
    # main() without args, cfg is parsed from sys.argv -> call from CLI
    # main(cfg) run directly with Config() -> used for tests
    # main(which, **kwargs) e.g. main("build", db_prefix="test", ...) -> generate config and run
    
    # flag to exit proper code if direct call from cli
    cli = False

    if cfg is None: 
        if which is None: cli = True
        cfg = Config(which, **kwargs)

    # Validate
    if not cfg.validate(): sys.exit(1)

    # Set paths
    if not cfg.set_paths(): sys.exit(1)

    tx_total = time.time()
    
    print_log("- - - - - - - - - -", cfg.quiet)
    print_log("   _  _  _  _  _   ", cfg.quiet)
    print_log("  (_|(_|| |(_)| |  ", cfg.quiet)
    print_log("   _|   v. "+ str(cfg.version), cfg.quiet)
    print_log("- - - - - - - - - -", cfg.quiet)

    if cfg.which=='build':
        ret=build(cfg)
    elif cfg.which=='update': 
        ret=update(cfg) 
    elif cfg.which=='classify':
        ret=classify(cfg)
    elif cfg.which=='report':
        ret=report(cfg)

    print_log("Total elapsed time: " + str("%.2f" % (time.time() - tx_total)) + " seconds.", cfg.quiet)
    
    if cli:
        sys.exit(0 if ret else 1)
    else:
        return ret

if __name__ == '__main__':
    main()

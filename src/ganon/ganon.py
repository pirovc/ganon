#!/usr/bin/env python3
import sys, time

from ganon.build_update import build, build_custom, update, update_custom
from ganon.classify import classify
from ganon.report import report
from ganon.table import table
from ganon.config import Config
from ganon.util import print_log

def main(which: str=None, cfg=None, **kwargs):
    # 3 entry points: 
    # main() without args, cfg is parsed from sys.argv
    # main(which, **kwargs) -> main("build", db_prefix="test", ...) generate config and run
    # main(cfg) run directly with Config()

    if cfg is None: cfg = Config(which, **kwargs)

    # Validate
    if not cfg.validate(): return False

    # Set paths
    if not cfg.set_paths(): return False

    tx_total = time.time()
    
    print_log("- - - - - - - - - -", cfg.quiet)
    print_log("   _  _  _  _  _   ", cfg.quiet)
    print_log("  (_|(_|| |(_)| |  ", cfg.quiet)
    print_log("   _|   v. "+ str(cfg.version), cfg.quiet)
    print_log("- - - - - - - - - -", cfg.quiet)

    if cfg.which=='build':
        ret=build(cfg)
    if cfg.which=='build-custom':
        ret=build_custom(cfg)
    elif cfg.which=='update': 
        ret=update(cfg) 
    elif cfg.which=='update-custom': 
        ret=update_custom(cfg) 
    elif cfg.which=='classify':
        ret=classify(cfg)
    elif cfg.which=='report':
        ret=report(cfg)
    elif cfg.which=='table':
        ret=table(cfg)

    print_log("Total elapsed time: " + str("%.2f" % (time.time() - tx_total)) + " seconds.", cfg.quiet)
    return ret

def main_cli():
    sys.exit(0 if main() else 1)

if __name__ == '__main__':
    sys.exit(0 if main() else 1)
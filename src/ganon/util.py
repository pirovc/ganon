import sys, subprocess, shlex, shutil, os, time

def run(cmd, print_stderr: bool=False, shell: bool=False, exit_on_error: bool=True):
    errcode=0
    stdout=""
    stderr=""
    try:
        process = subprocess.Popen(shlex.split(cmd) if not shell else cmd, 
                                    shell=shell, 
                                    universal_newlines=True, 
                                    stdout=subprocess.PIPE, 
                                    stderr=subprocess.PIPE)   
        stdout, stderr = process.communicate() # wait for the process to terminate
        errcode = process.returncode
        if errcode!=0: raise Exception()
        if print_stderr and stderr: print_log(stderr)
 
    #except OSError as e: # The most common exception raised is OSError. This occurs, for example, when trying to execute a non-existent file. Applications should prepare for OSError exceptions.
    #except ValueError as e: #A ValueError will be raised if Popen is called with invalid arguments.
    except Exception as e:
        print_log('The following command failed to run:\n'+cmd)
        print_log(str(e))
        print_log("Error code: "+str(errcode))
        print_log("Out: ")
        if stdout: print_log(stdout)
        print_log("Error: ")
        if stderr: print_log(stderr)
        if exit_on_error: sys.exit(errcode)

    return stdout, stderr

def print_log(text, quiet: bool=False):
    if not quiet:
        sys.stderr.write(text+"\n")
        sys.stderr.flush()

def set_tmp_folder(fld):
    # Create temporary working directory
    if os.path.exists(fld): 
        print_log("ERROR: temp folder already exists " + os.path.abspath(fld))
        return False
    else:
        os.makedirs(fld)
        return True

def rm_tmp_folder(fld):
    shutil.rmtree(fld)

def check_files(files):
    checked_files = [file for file in files if os.path.isfile(file) and os.path.getsize(file) > 0]
    if len(checked_files)<len(files):
        print_log(str(len(files)-len(checked_files)) + " input file(s) could not be found/empty")
    return checked_files

def check_db(prefix):
    for db_file_type in [".ibf", ".map", ".tax", ".gnn"]:
        if not os.path.isfile(prefix+db_file_type):
            print_log("Incomplete database [" + prefix  + "] (.ibf, .map, .tax and .gnn)")
            return False
    return True

def set_taxdump_files(taxdump_file, tmp_output_folder, quiet):
    if not taxdump_file:
        ncbi_nodes_file, ncbi_names_file, ncbi_merged_file = unpack_taxdump(get_taxdump(tmp_output_folder, quiet), tmp_output_folder, quiet)
    elif taxdump_file[0].endswith(".tar.gz"):
        ncbi_nodes_file, ncbi_names_file, ncbi_merged_file = unpack_taxdump(taxdump_file[0], tmp_output_folder, quiet)
    else:
        ncbi_nodes_file = taxdump_file[0]
        ncbi_names_file = taxdump_file[1]
        ncbi_merged_file =  taxdump_file[2] if len(taxdump_file)==3 else ""

    return ncbi_nodes_file, ncbi_merged_file, ncbi_names_file

def get_taxdump(tmp_output_folder, quiet):
    tx = time.time()
    print_log("Downloading taxdump", quiet)
    taxdump_file = tmp_output_folder+'taxdump.tar.gz'
    run_wget_taxdump_cmd = 'wget -qO {0} "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"'.format(taxdump_file)
    stdout, stderr = run(run_wget_taxdump_cmd, print_stderr=True)
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", quiet)
    return taxdump_file

def unpack_taxdump(taxdump_file, tmp_output_folder, quiet):
    tx = time.time()
    print_log("Unpacking taxdump", quiet)
    unpack_taxdump_cmd = 'tar xf {0} -C "{1}" nodes.dmp merged.dmp names.dmp'.format(taxdump_file, tmp_output_folder)
    stdout, stderr = run(unpack_taxdump_cmd, print_stderr=True)
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", quiet)
    return tmp_output_folder+'nodes.dmp', tmp_output_folder+'names.dmp', tmp_output_folder+'merged.dmp'

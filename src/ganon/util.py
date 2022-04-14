import sys, subprocess, shlex, shutil, os, time
import urllib.request

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
        if exit_on_error and errcode!=0: 
            raise Exception()
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
        sys.exit(errcode)

    return stdout, stderr

def print_log(text, quiet: bool=False):
    if not quiet:
        sys.stderr.write(text+"\n")
        sys.stderr.flush()

def set_tmp_folder(fld, restart):
    # Create temporary working directory
    if os.path.exists(fld): 
        if restart:
            rm_tmp_folder(fld)
        else:
            print_log("ERROR: temp folder already exists " + os.path.abspath(fld))
            return False
    
    os.makedirs(fld)
    return True

def rm_tmp_folder(fld):
    shutil.rmtree(fld)

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

def validate_input_files(input_files_folder, input_extension, quiet):
    """
    given a list of input files and/or folders and an file extension
    check for valid files and return them in a set
    """
    valid_input_files = set()
    for i in input_files_folder:
        if check_file(i):
            valid_input_files.add(i)
        elif os.path.isdir(i):
            if not input_extension:
                print_log("--input-extension is required when using folders in the --input. Skipping: " + i, quiet)
                continue
            files_in_dir=0
            for file in os.listdir(i):
                if file.endswith(input_extension):
                    f = os.path.join(i, file)
                    if check_file(f):
                        files_in_dir += 1
                        valid_input_files.add(f)
            print_log(str(files_in_dir) + " valid file(s) [--input-extension " + input_extension + "] found in " + i, quiet)
        else:
            print_log("Skipping invalid file/folder: " + i, quiet)

    print_log("Total valid files: " + str(len(valid_input_files)), quiet)
    return valid_input_files

def check_file(file):
    if os.path.isfile(file) and os.path.getsize(file) > 0:
        return True
    else:
        return False

def download(urls: list, output_prefix: str):
    """
    Parameters:
    * **urls** *[list]*: List of urls to download
    * **output_prefix** *[str]*: Output directory to save files

    Returns:
    * list of files saved
    """
    files = []
    for url in urls:
        outfile = output_prefix + os.path.basename(url)
        urlstream = urllib.request.urlopen(url)
        with open(outfile, 'b+w') as f:
            f.write(urlstream.read())
        urlstream.close()
        files.append(outfile)
    return files
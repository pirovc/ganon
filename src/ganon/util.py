import sys
import subprocess
import shlex
import shutil
import os
import time
import urllib.request


def run(cmd, ret_stdout: bool=False, shell: bool=False, quiet: bool=False):
    errcode = 1
    stdout = None

    try:
        # print stdout to stderr if not captured (default log ganon)
        process = subprocess.Popen(shlex.split(cmd) if not shell else cmd,
                                   shell=shell,
                                   text=True,
                                   stdout=subprocess.PIPE if ret_stdout else sys.stderr,
                                   stderr=None if quiet else sys.stderr)
        if ret_stdout:
            stdout, stderr = process.communicate()
            if stderr and not quiet:
                print_log(stderr)
        else:
            process.wait()

        errcode = process.returncode
        if errcode != 0:
            raise Exception()

    except Exception as e:
        print_log("The following command failed to run:\n" + cmd)
        print_log(str(e))
        print_log("Error code: " + str(errcode))
        sys.exit(errcode)

    return stdout


def print_log(text, quiet: bool=False):
    if not quiet:
        sys.stderr.write(text+"\n")
        sys.stderr.flush()


def set_out_folder(fld, restart):
    # Create working directory
    if os.path.exists(fld):
        if restart:
            rm_folder(fld)
        else:
            print_log("ERROR: temp folder already exists " + os.path.abspath(fld))
            return False

    os.makedirs(fld)
    return True


def set_out_files(prefix, ext, restart):
    for e in ext:
        file = prefix + "." + e
        if os.path.exists(file):
            if restart:
                os.remove(file)
            else:
                print_log("ERROR: output file already exists " + os.path.abspath(file))
                return False
    return True


def rm_folder(fld):
    shutil.rmtree(fld)


def set_taxdump_files(taxdump_file, tmp_output_folder, quiet):
    if not taxdump_file:
        ncbi_nodes_file, ncbi_names_file, ncbi_merged_file = unpack_taxdump(get_taxdump(tmp_output_folder, quiet),
                                                                            tmp_output_folder,
                                                                            quiet)
    elif taxdump_file[0].endswith(".tar.gz"):
        ncbi_nodes_file, ncbi_names_file, ncbi_merged_file = unpack_taxdump(taxdump_file[0], tmp_output_folder, quiet)
    else:
        ncbi_nodes_file = taxdump_file[0]
        ncbi_names_file = taxdump_file[1]
        ncbi_merged_file = taxdump_file[2] if len(taxdump_file) == 3 else ""

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
            files_in_dir = 0
            for file in os.listdir(i):
                if file.endswith(input_extension):
                    f = os.path.join(i, file)
                    if check_file(f):
                        files_in_dir += 1
                        valid_input_files.add(f)
            print_log(str(files_in_dir) + " valid file(s) [--input-extension " + input_extension +
                      "] found in " + i, quiet)
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

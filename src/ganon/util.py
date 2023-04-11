import sys
import subprocess
import shlex
import os
import urllib.request
from pathlib import Path


def run(cmd, ret_stdout: bool=False, shell: bool=False, quiet: bool=False):
    errcode = 1
    stdout = None

    try:
        # print stdout to stderr if not captured (default log ganon)
        process = subprocess.Popen(shlex.split(cmd) if not shell else cmd,
                                   shell=shell,
                                   universal_newlines=True,  # text= (from py3.7)
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


def logo(version):
    logo = ""
    logo += "- - - - - - - - - -\n"
    logo += "   _  _  _  _  _   \n"
    logo += "  (_|(_|| |(_)| |  \n"
    logo += "   _|   v. " + str(version) + "\n"
    logo += "- - - - - - - - - -"
    return logo


def print_log(text, quiet: bool=False):
    if not quiet:
        sys.stderr.write(text+"\n")
        sys.stderr.flush()


def rm_files(files):
    if isinstance(files, str):
        files = [files]
    for f in files:
        if os.path.isfile(f):
            os.remove(f)


def validate_input_files(input_files_folder, input_extension, quiet, input_recursive: bool = False):
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

            if input_recursive:
                for path in Path(i).rglob('*' + input_extension):
                    f = str(path.joinpath())
                    if check_file(f):
                        files_in_dir += 1
                        valid_input_files.add(f)
            else:
                for file in os.listdir(i):
                    if file.endswith(input_extension):
                        f = os.path.join(i, file)
                        if check_file(f):
                            files_in_dir += 1
                            valid_input_files.add(f)

            print_log(str(files_in_dir) + " valid file(s) [--input-extension " + input_extension +
                      (", --input-recursive" if input_recursive else "") + 
                      "] found in " + i, quiet)
        else:
            print_log("Skipping invalid file/folder: " + i, quiet)

    print_log("Total valid files: " + str(len(valid_input_files)) + "\n", quiet)
    return valid_input_files


def check_file(file):
    """
    Check if file exists and it's not empty
    """
    if os.path.isfile(file) and os.path.getsize(file) > 0:
        return True
    else:
        return False


def check_folder(folder):
    """
    Check if folder exists and it's not empty
    """
    if os.path.isdir(folder) and len(os.listdir(folder)) > 0:
        return True
    else:
        return False


def save_state(state, folder):
    Path(folder + state).touch()


def load_state(state, folder):
    return os.path.isfile(folder + state)


def set_output_folder(db_prefix):
    """
    set general working directory for downloads and temporary files
    """
    return db_prefix + "_files/"


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

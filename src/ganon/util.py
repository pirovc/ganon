import sys, subprocess, shlex

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
        print_log('The following command failed to execute:\n'+cmd)
        print_log(str(e))
        print_log("Error code: "+str(errcode))
        print_log("Out: ")
        if stdout: print_log(stdout)
        print_log("Error: ")
        if stderr: print_log(stderr)
        if exit_on_error: sys.exit(errcode)

    return stdout, stderr

def print_log(text):
    sys.stderr.write(text+"\n")
    sys.stderr.flush()
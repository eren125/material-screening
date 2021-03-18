import os

def bash(command):
    return os.system(command)

def copy_dir(path_to_src,dirname_src,path_to_dst,dirname_dst,delete=False):
    """
    Using rsync to copy a directory to a given location with a custom name
    --args--
        path_to_src (str): location of the source directory
        dirname_src (str): name of the source directory
        path_to_dst (str): location of the destination directory
        dirname_dst (str): name of the destination directory
    --option--
        delete = whether you want to delete non existing files in the destination file
    """
    option = ''
    if delete==True:
        option = ' --delete-during'
    print("Copying %s to the destination directory %s"%(dirname_src,path_to_dst))
    return bash("rsync -rtv%s %s/ %s/"%(option,os.path.join(path_to_src,dirname_src),os.path.join(path_to_dst,dirname_dst)))

def copy_file(path_to_src,filename_src,path_to_dst,filename_dst):
    """
    Using rsync to copy a directory to a given location with a custom name
    --args--
        path_to_src (str): location of the source directory
        filename_src (str): name of the source file
        path_to_dst (str): location of the destination directory
        filename_dst (str): name of the destination file
    """
    return bash("cp %s %s"%(os.path.join(path_to_src,filename_src),os.path.join(path_to_dst,filename_dst)))

def sed(pattern, replace, source):
    """
    Using sed to replace a pattern with another one
    --args--
        pattern (str): pattern to match
        replace (str): replacement str
        source  (str): input filename
    """
    return bash('sed -i "s/%s/%s/g" %s'%(pattern,replace,source))

def run_shell_script(shell_script,verbose=1):
    """Runs a bash script and print out stdout if verbose == 1
    Make sure to put <cd "$(dirname "$0")"> in the bash script if the script is reading files using relative paths (e.g. ../filename)

    Args:
        shell_script (str): Path to shell script
        verbose = 0 or 1 depending on whether ou want to print out messages from the bash script      
    """
    return bash("bash %s"%shell_script)

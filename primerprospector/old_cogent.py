from tempfile import gettempdir

class FilePath(str):
    """ Hold paths for proper handling
        Paths in this sense are filenames, directory paths, or filepaths. 
        Some examples include:
         file.txt
         ./path/to/file.txt
         ./path/to/dir/
         /path/to/file.txt
         .
         /
        The purpose of this class is to allow all paths to be handled the
         same since they sometimes need to be treated differently than 
         simple strings. For example, if a path has a space in it, and it
         is being passed to system, it needs to be wrapped in quotes. But,
         you wouldn't want it as a string wrapped in quotes b/c, e.g.,
         isabs('"/absolute/path"') == False, b/c the first char is a ", not
         a /.
        * This would make more sense to call Path, but that conflicts with
            the ResultPath.Path attribute. I'm not sure what to do about this
            and want to see what others think. Once finalized, a global
            replace should take care of making the switch.
    """
    def __new__(cls,path):
        try:
            return str.__new__(cls, path.strip('"'))
        except AttributeError:
            return str.__new__(cls,'')

    def __str__(self):
        """ wrap self in quotes, or return the empty string if self == '' """
        if self == '': return ''
        return ''.join(['"',self,'"'])

    def __add__(self,other):
        return FilePath(''.join([self,other]))

def get_tmp_filename(tmp_dir=gettempdir(), prefix="tmp", suffix=".txt",
                     result_constructor=FilePath):
    """ Generate a temporary filename and return as a FilePath object
    
        tmp_dir: the directory to house the tmp_filename (default: '/tmp')
        prefix: string to append to beginning of filename (default: 'tmp')
            Note: It is very useful to have prefix be descriptive of the
            process which is creating the temporary file. For example, if 
            your temp file will be used to build a temporary blast database, 
            you might pass prefix=TempBlastDB
        suffix: the suffix to be appended to the temp filename (default '.txt')
        result_constructor: the constructor used to build the result filename
            (default: cogent.app.parameters.FilePath). Note that joining 
            FilePath objects with one another or with strings, you must use
            the + operator. If this causes trouble, you can pass str as the 
            the result_constructor.
    """
    # check not none
    if not tmp_dir:
        tmp_dir = ""
    # if not current directory, append "/" if not already on path
    elif not tmp_dir.endswith("/"):
        tmp_dir += "/"

    chars = "abcdefghigklmnopqrstuvwxyz"
    picks = chars + chars.upper() + "0123456790"
    return result_constructor(tmp_dir) + result_constructor(prefix) +\
        result_constructor("%s%s" % \
        (''.join([choice(picks) for i in range(20)]),suffix))

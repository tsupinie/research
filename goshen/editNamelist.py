
import re

def _namelistToString(value):
    if type(value) in [ int, float ]:
        str_value = str(value)
    elif type(value) in [ bool ]:
        str_value = ".%s." % str(value).lower()
    else:
        str_value = "'%s'" % str(value)
    return str_value

def _namelistFromString(str_value):
    if str_value[0] == "'" and str_value[-1] == "'": 
        # Remove the quotes and leave the value as a string
        value = str_value[1:-1]
    elif str_value[0] == "." and str_value[-1] == ".":
        # Logical true or false.
        value = (str_value[1:-1] == "true") or (str_value[1:-1] == "TRUE")
    else:
        # Try to convert the value to an integer, if it won't work, convert it to a float
        try:
            value = int(str_value)
        except ValueError:
            value = float(str_value)

    return value

def _searchLine(line, **kwargs):
    for parameter, replacement in kwargs.iteritems():
        # Only pull out the regular expressions if we have reason to believe the parameter might be on this line.
        if line.find(parameter) > 0:
            # Escape any special characters in the parameter
            parameter_escaped = re.escape(parameter)

            # Looks for the parameter name preceded by a space, a comma, or the beginning of the line.
            match = re.search(r"(?<=[\s,])%s(?=[\s]*=[\s]*[\S]+[\s]*,)" % parameter_escaped, line)
            if match is None: match = re.search(r"^%s(?=[\s]*=[\s]*[\S]+[\s]*,)" % parameter_escaped, line)

            if line[0] != '!' and match is not None:
                if type(replacement) in [ list, tuple ]:
                    str_value = ", ".join([ _namelistToString(r) for r in replacement ])
                else:
                    str_value = _namelistToString(replacement)

#               print "Replacing %s with %s" % (parameter, replacement)

                # Looks for a string that matches (parameter = a comma-separated list of either strings, booleans, or numbers) followed by a bareword (another variable name) or a newline.
                line = re.sub(r"%s[\s]*=[\s]*(?:(?:'[^']+'[\s]*,?[\s]*)+|(?:\.(?:true|TRUE|false|FALSE)\.[\s]*,?[\s]*)+|(?:[\dEe.-]+[\s]*,?[\s]*)+?)(?=[A-Za-z][\w(,)]*|" % parameter_escaped + "\n)", 
                    "%s = %s," % (parameter, str_value), line)

    return line

def joinSplicedStrings(spliced_list, delimiter=", "):
    idx_start = 0
    idx = idx_start
    while idx < len(spliced_list):
        if not ("".join(spliced_list[idx_start:(idx + 1)]).count("'") % 2) or idx == len(spliced_list) - 1:
            spliced_list[idx_start:(idx + 1)] = [ delimiter.join(spliced_list[idx_start:(idx + 1)]) ]
            idx_start += 1
            idx = idx_start
        else:
            idx += 1
    return spliced_list

#print joinSplicedStrings(['      iusesng(1,1)=0', 'iusesng(1,2)=1', 'iusesng(1,3)=0', 'iusesng(1,4)=0', '\n'])
#breakme

def editNamelistFile(file_name_src, file_name_dest, **kwargs):
    """
    editNamelistFile
    Purpose:    Takes a FORTRAN-90 namelist file and edits it either using the keys/values of kwargs.  The special key __file_name_values__ 
                    can be used to find all the values from another namelist file and use them.
    Parameters: file_name_src [type=str]
                    The name of the namelist file to read from
                file_name_dest [type=str]
                    The name of the namelist file to write to
                kwargs [type=dict]
                    Expanded dictionary (sequence of key/value pairs, e.g. var1="string", var2=1) of parameters to change in the namelist file.
    Example Use:
                editNamelistFile("namelist.input", "namelist.input", timestep=10, dx=(3000, 1000), dy=(3000, 1000))
                    Take namelist.input and change the timestep variable to "10", change the dx variable to be "3000, 1000", and dy to be "3000, 1000"
                editNamelistFile("namelist.input", "namelist.test.input" __file_name_values__="namelist.old.input")
                    Take namelist.input and replace all the values in it with the values in namelist.old.input.  Put the result in namelist.test.input .
    """
    if '__file_name_values__' in kwargs:
        # We have a file to look in for the values
        file_name_values = kwargs['__file_name_values__']
        del kwargs['__file_name_values__']

        file_values = open(file_name_values, 'r')
        for line in file_values:
            if line[0] != '!' and line.find(",") > -1:
                bits = re.split(r",[\s]*(?=[A-Za-z][\w(,)]*|" + "\n)", line)
                bits = joinSplicedStrings(bits)

                for kv_pair in bits:
                    if kv_pair.find("=") > -1:
                        key, value = [ bit.strip() for bit in kv_pair.split("=") ]

                        if value.find(",") > 0:
                            value = value.split(",")
                            value = joinSplicedStrings(value, ",")
                            value = tuple([ _namelistFromString(v.strip()) for v in value ])
                        else:
                            value = _namelistFromString(value.strip())

                        kwargs[key] = value

        file_values.close()

    if file_name_src != file_name_dest:
        # If the files do not have the same name, open them both at the same time (saves on memory)
        file_src = open(file_name_src, 'r')
        file_dest = open(file_name_dest, 'w')

        for line in file_src:
            line = _searchLine(line, **kwargs)
            file_dest.write(line)

        file_dest.close()
        file_src.close()
    else:
        # If the files do have the same name, open one after the other and load the entire source file into memory.
        file_src = open(file_name_src, 'r')
        buffer = []
        for line in file_src:
            line = _searchLine(line, **kwargs)
            buffer.append(line)
        file_src.close()

        file_dest = open(file_name_dest, 'w')
        for line in buffer:
            file.dest.write(line)
        file_dest.close()

    return

if __name__ == "__main__":
#   editNamelistFile("arps.input", "arps.test.input")
    editNamelistFile("/home/tsupinie/arps5.3/input/arps.input", "ext2arps_new.3km.input", __file_name_values__="ext2arps.3km.input")
#   editNamelistFile("arpsenkf.template.input", "arpsenkf.input", __file_name_values__="arpsenkf.old.input")
#   editNamelistFile("arps.input", "arps.test.input", inifile='/brashear/supinie/cov_infl/30May2004/enf001.hdf000300')

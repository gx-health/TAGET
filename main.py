import sys

if sys.version_info[0] < 3:
    sys.stderr.write("ERROR: python version > 3.0.0 required ! Now is {}.{}.{}\n".format(sys.version_info[0], sys.version_info[1],sys.version_info[2]))
    exit()
    
import TransAnnot


TransAnnot.main()
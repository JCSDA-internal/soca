#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Replace xxxx and XXXX in templace ModelX
Example:
     $  template2model.py mom5cice5 MOM5CICE55'
     will replace xxxx with mom5cice5 and XXXX with MOM5CICE5 
Todo:
     *
"""

import glob
import sys
import numpy as np

if __name__ == '__main__':

    print sys.argv[0],': Replaces xxxx and XXXX with user defined strings in *.*'
    try:
        str2replace=sys.argv[1]      # String to replace
        replaceby=sys.argv[2]      # replacement string
    except:
        print 'Usage example: changestr.py mom5cice5_ soca_'
        sys.exit("Error: Wrong arguments")

    flist=np.append(glob.glob("./**"), glob.glob("./**/**"))
    flist=np.append(flist,glob.glob("./**/**/**"))
    flist.sort()
    print np.shape(flist)
    for fname in flist:
        print fname
        if (fname!=sys.argv[0]):
            try:
                with open(fname, 'r') as file :
                    filedata = file.read()
                    filedata = filedata.replace(str2replace, replaceby)
                with open(fname, 'w') as file:
                    file.write(filedata)
            except:
                print 'Skipping ', fname,'\n'

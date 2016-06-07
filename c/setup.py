from distutils.core import setup, Extension
import os
import numpy as np

module = 'precip_tools'
setup(name=module, version='1.0',  \
      ext_modules=[Extension(module, ['precip_tools.c'],
        include_dirs=[np.get_include()]
        )])
# setup(name=module, version='1.0',  \
#       ext_modules=[Extension(module, ['precip_tools.c']
#         )])
# Move it up to current directory
os.system('cp build/lib.mac*/%s.so ../%s.so'%(module, module))

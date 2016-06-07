from distutils.core import setup, Extension
import numpy
import os
# define the extension module
cos_module_np = Extension('cos_module_np', sources=['cos_module_np.c'],
                          include_dirs=[numpy.get_include()])

# run the setup
setup(ext_modules=[cos_module_np])
module = 'cos_module_np'

os.system('cp build/lib.mac*/%s.so ../%s.so'%(module, module))

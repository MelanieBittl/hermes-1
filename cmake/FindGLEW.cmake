#
# GLEW
#

if(WIN64)
  FIND_LIBRARY(GLEW_LIBRARY NAMES glew32 GLEW PATHS ${GLEW_ROOT}/lib/x64 ${GLEW_ROOT}/lib)
else(WIN64)  
  FIND_LIBRARY(GLEW_LIBRARY NAMES glew32 GLEW PATHS ${GLEW_ROOT}/lib /usr/lib /usr/local/lib /usr/lib64 /usr/local/lib64 /usr/lib/x86_64-linux-gnu )

endif(WIN64)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GLEW DEFAULT_MSG GLEW_LIBRARY)

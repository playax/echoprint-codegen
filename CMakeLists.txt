cmake_minimum_required(VERSION 3.7)
project(echoprint_codegen)

set(CMAKE_CXX_STANDARD 11)

set(OPTFLAGS "-O3 -DBOOST_UBLAS_NDEBUG -DNDEBUG -DUNHASHED_CODES")
set(BOOST_CFLAGS "-I/usr/local/include/boost -I/usr/local/Cellar/boost/1.64.0_1/include/")
set(JNI_CFLAGS "-I/Library/Java/JavaVirtualMachines/jdk1.8.0_60.jdk/Contents/Home/include/ -I/Library/Java/JavaVirtualMachines/jdk1.8.0_60.jdk/Contents/Home/include/darwin/ -I/Library/Java/JavaVirtualMachines/jdk1.8.0_60.jdk/Contents/Home/include/linux/")

set(TAGLIB_CFLAGS "-I/usr/local/Cellar/taglib/1.9.1/include/taglib")
set(TAGLIB_LIBS "-L/usr/local/Cellar/taglib/1.9.1/lib -ltag")


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BOOST_CFLAGS} ${TAGLIB_CFLAGS} ${JNI_CFLAGS} -fPIC ${OPTFLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -fPIC ${OPTFLAGS}")

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${TAGLIB_LIBS} -lz -lpthread ${OPTFLAGS}")
set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${TAGLIB_LIBS} -lz -lpthread ${OPTFLAGS}")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${TAGLIB_LIBS} -lz -lpthread ${OPTFLAGS}")
set(CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} ${TAGLIB_LIBS} -lz -lpthread ${OPTFLAGS}")

set(SOURCE_FILES
        src/AudioBufferInput.cxx
        src/AudioBufferInput.h
        src/AudioStreamInput.cxx
        src/AudioStreamInput.h
        src/Base64.cxx
        src/Base64.h
        src/Codegen.cxx
        src/Codegen.h
        src/Common.h
        src/echoprint.cxx
        src/echoprint.h
        src/File.h
        src/Fingerprint.cxx
        src/Fingerprint.h
        src/functions.cxx
        src/functions.h
        src/main.cxx
        src/MatrixUtility.cxx
        src/MatrixUtility.h
        src/Metadata.cxx
        src/Metadata.h
        src/Params.h
        src/SubbandAnalysis.cxx
        src/SubbandAnalysis.h
        src/Whitening.cxx
        src/Whitening.h
        src/win_funcs.h
        src/win_unistd.h)

add_executable(echoprint_codegen ${SOURCE_FILES})
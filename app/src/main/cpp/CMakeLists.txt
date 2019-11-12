# For more information about using CMake with Android Studio, read the
# documentation: https://d.android.com/studio/projects/add-native-code.html

cmake_minimum_required(VERSION 3.4.1)

add_library(native-lib
            SHARED
            native-lib.cpp )

find_library(log-lib log)

set(TOOLCHAIN $ENV{HOME}/src/gnuradio-android/toolchain)

target_include_directories(native-lib PUBLIC
        ${TOOLCHAIN}/${ANDROID_ABI}/include
        ${TOOLCHAIN}/${ANDROID_ABI}/include/boost-1_69)

if(${ANDROID_ABI} STREQUAL "armeabi-v7a")
    set(BOOST_ARCH "a32")
else()
    set(BOOST_ARCH "a64")
endif()

target_link_libraries(native-lib

        ${TOOLCHAIN}/${ANDROID_ABI}/lib/libboost_filesystem-clang-mt-${BOOST_ARCH}-1_69.a
        ${TOOLCHAIN}/${ANDROID_ABI}/lib/libvolk.a
        ${log-lib})
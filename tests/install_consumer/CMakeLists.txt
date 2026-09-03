cmake_minimum_required(VERSION 3.20)
project(LinearRaccessInstallConsumer LANGUAGES CXX)

find_package(LinearRaccess 0.1 REQUIRED)

add_executable(consumer main.cpp)
target_link_libraries(consumer PRIVATE LinearRaccess::linearraccess)

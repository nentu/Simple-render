cmake_minimum_required(VERSION 3.10)
project(path_tracer CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Оптимизации для релиза
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -ffast-math -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")

# Пытаемся найти Embree 4, затем Embree 3
find_package(embree 4.0 QUIET)
if(NOT embree_FOUND)
    find_package(embree 3.0 QUIET)
    if(NOT embree_FOUND)
        message(FATAL_ERROR "Embree not found. Install with: sudo apt install libembree3-dev")
    else()
        message(STATUS "Found Embree 3")
        # Определяем макрос для кода
        add_definitions(-DEMBREE3)
    endif()
else()
    message(STATUS "Found Embree 4")
endif()

# Исходные файлы
set(SOURCES main.cpp)

# Создаем исполняемый файл
add_executable(path_tracer ${SOURCES})

# Подключаем Embree
target_link_libraries(path_tracer PRIVATE embree)

# Включаем OpenMP для параллельного рендеринга
find_package(OpenMP REQUIRED)
target_link_libraries(path_tracer PRIVATE OpenMP::OpenMP_CXX)

# Предупреждения компилятора
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(path_tracer PRIVATE -Wall -Wextra -Wpedantic)
endif()
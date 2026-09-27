cmake_minimum_required(VERSION 3.16)
project(bToMuHnlSim LANGUAGES CXX)

# ---------------------------------------------------------------------------
# Estandar de C++
#
# Geant4 y ROOT deben compilarse con el MISMO estandar que esta aplicacion, o
# aparecen errores de linking raros (sobre todo por std::string en las
# branches del TTree). Si tu ROOT fue compilado con C++20, cambia esto a 20.
# Para saberlo:  root-config --cflags | grep -o 'std=c++[0-9]*'
# ---------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ standard (debe coincidir con ROOT/Geant4)")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Para clangd / IDEs
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ---------------------------------------------------------------------------
# Geant4
#
# Requiere haber hecho:  source <geant4-install>/bin/geant4.sh
# (o pasar -DGeant4_DIR=<install>/lib64/cmake/Geant4)
#
# ui_all/vis_all no son estrictamente necesarios para esta app (corre en batch,
# sin macros ni visualizacion), pero no estorban y dejan la puerta abierta a
# agregar un visor despues. Si tu instalacion no los tiene, borra esos dos
# argumentos.
# ---------------------------------------------------------------------------
find_package(Geant4 REQUIRED ui_all vis_all)

# ---------------------------------------------------------------------------
# ROOT
#
# Requiere haber hecho:  source <root-install>/bin/thisroot.sh
# (o pasar -DROOT_DIR=<install>/cmake)
#
# Componentes usados: TFile/TTree (MuonShieldFieldMap, TrackingAction,
# bToMuHnlSim) y los typedefs de Rtypes.h.
# ---------------------------------------------------------------------------
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree)

# ---------------------------------------------------------------------------
# Pythia8
#
# Pythia no instala un archivo de configuracion de CMake propio, asi que hay
# que buscarlo a mano. Orden de busqueda:
#   1. pkg-config (si instalaste Pythia con ./configure --prefix=... y tienes
#      pythia8.pc en PKG_CONFIG_PATH)
#   2. las variables de entorno/cache PYTHIA8 o PYTHIA8_DIR
#   3. rutas estandar del sistema
#
# Si no lo encuentra, configura con:
#   cmake -S . -B build -DPYTHIA8_DIR=/ruta/a/pythia8
# ---------------------------------------------------------------------------
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(PYTHIA8_PC QUIET pythia8)
endif()

set(PYTHIA8_HINT_DIRS
    ${PYTHIA8_DIR}
    $ENV{PYTHIA8_DIR}
    $ENV{PYTHIA8}
    ${PYTHIA8_PC_PREFIX}
)

find_path(PYTHIA8_INCLUDE_DIR
  NAMES Pythia8/Pythia.h
  HINTS ${PYTHIA8_HINT_DIRS} ${PYTHIA8_PC_INCLUDE_DIRS}
  PATH_SUFFIXES include
)

find_library(PYTHIA8_LIBRARY
  NAMES pythia8 Pythia8
  HINTS ${PYTHIA8_HINT_DIRS} ${PYTHIA8_PC_LIBRARY_DIRS}
  PATH_SUFFIXES lib lib64
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Pythia8
  REQUIRED_VARS PYTHIA8_LIBRARY PYTHIA8_INCLUDE_DIR
  FAIL_MESSAGE "No se encontro Pythia8. Configura con -DPYTHIA8_DIR=/ruta/a/pythia8 o exporta PYTHIA8=/ruta/a/pythia8"
)

if(NOT TARGET Pythia8::Pythia8)
  add_library(Pythia8::Pythia8 UNKNOWN IMPORTED)
  set_target_properties(Pythia8::Pythia8 PROPERTIES
    IMPORTED_LOCATION "${PYTHIA8_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${PYTHIA8_INCLUDE_DIR}"
  )
endif()

# ---------------------------------------------------------------------------
# Fuentes
#
# Layout plano: todos los .cc y .hh viven en la raiz del repo. Se listan
# explicitamente (en vez de usar GLOB) para que agregar o quitar un archivo
# obligue a reconfigurar, que es lo que uno quiere.
# ---------------------------------------------------------------------------
set(BTOMUHNLSIM_SOURCES
  src/BellMagneticField.cxx
  src/DetectorConstruction.cxx
  src/MuonShieldFieldMap.cxx
  src/MuonShieldGeometry.cxx
  src/PhysicsList.cxx
  src/PrimaryGeneratorAction.cxx
  src/SteppingAction.cxx
  src/TrackingAction.cxx
  bToMuHnlSim.cxx
)

add_executable(bToMuHnlSim ${BTOMUHNLSIM_SOURCES})

target_include_directories(bToMuHnlSim PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(bToMuHnlSim PRIVATE
  ${Geant4_LIBRARIES}
  ROOT::Core
  ROOT::RIO
  ROOT::Tree
  Pythia8::Pythia8
)

# Geant4_INCLUDE_DIRS/Geant4_DEFINITIONS vienen de find_package(Geant4); se
# aplican a mano en vez de con include(${Geant4_USE_FILE}) para no contaminar
# variables globales del proyecto.
target_include_directories(bToMuHnlSim PRIVATE ${Geant4_INCLUDE_DIRS})
target_compile_definitions(bToMuHnlSim PRIVATE ${Geant4_DEFINITIONS})

# ---------------------------------------------------------------------------
# Datos: el mapa de campo TRY_2026.root
#
# El ejecutable lo abre por ruta RELATIVA al directorio desde donde se corre
# (default "files/TRY_2026.root"). Copiar files/ junto al binario permite
# correrlo desde build/ sin pasar la ruta a mano.
# ---------------------------------------------------------------------------
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/files)
  add_custom_command(TARGET bToMuHnlSim POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${CMAKE_CURRENT_SOURCE_DIR}/files
            $<TARGET_FILE_DIR:bToMuHnlSim>/files
    COMMENT "Copiando files/ (mapa de campo TRY_2026.root) junto al ejecutable"
  )
endif()

install(TARGETS bToMuHnlSim RUNTIME DESTINATION bin)

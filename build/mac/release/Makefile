BASE=$(shell brew --prefix)
BOOST_VERSION=
BOOST_INCLUDE = $(BASE)/include
C_PLATFORM=-pthread
GPP=/usr/bin/clang++
C_OPTIONS= -DNDEBUG -std=c++14 -fvisibility=hidden -O3 -flto=thin
LDFLAGS= -flto=thin
BOOST_LIB_VERSION=
BOOST_STATIC=y

include ../../makefile_common

# Set path to Alembic and IlmBase
ABC_INC := -I/home/syoyo/local/alembic-dist/include
ABC_LDFLAGS := -L/home/syoyo/local/alembic-dist/lib -lAlembic

ILMBASE_INC := -I/home/syoyo/local/ilmbase-2.2.0-dist/include/OpenEXR

# Suppress some C++ warnings(in clang).
EXTRA_CXXFLAGS := -Wno-deprecated-register -Weverything -Wno-c++98-compat

all:
	clang++ -std=c++11 -g -o curve-combiner $(ABC_INC) $(ILMBASE_INC) $(EXTRA_CXXFLAGS) curve-combiner.cc $(ABC_LDFLAGS)


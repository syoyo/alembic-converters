# Set path to Alembic and IlmBase
ABC_INC := -I/usr/local/include
ABC_LDFLAGS := -L/usr/local/lib -lAlembic

ILMBASE_INC := -I/usr/local/include/OpenEXR

# Suppress some C++ warnings(in clang).
EXTRA_CXXFLAGS := -Wno-deprecated-register -Weverything

all:
	clang++ -std=c++11 -g -o abc2cyhair $(ABC_INC) $(ILMBASE_INC) $(EXTRA_CXXFLAGS) abc2cyhair.cc cyhair-writer.cc $(ABC_LDFLAGS)


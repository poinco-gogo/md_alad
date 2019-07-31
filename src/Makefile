TARGET=md_alad_ewald.x

CXX=icpc

CXXFLAGS= -O3 -std=c++11 -I../include -I../include/uiuc

LIBS= -mkl

SRCS = Random.cpp dcdplugin.cpp Atom.cpp NAMDBin.cpp PDBUtil.cpp PDB.cpp PSF.cpp LoadParm.cpp Option.cpp System.cpp ComputeBond.cpp ComputeAngle.cpp ComputeLJ.cpp ComputeES.cpp Energy.cpp Integrator.cpp Output.cpp md_alad_ewald.cpp

OBJS = $(SRCS:.cpp=.o)

HEADS= common.hpp Random.hpp Lattice.hpp Atom.hpp NAMDBin.hpp Bond.hpp Angle.hpp PDBUtil.hpp PDB.hpp PSF.hpp LoadParm.hpp Option.hpp System.hpp ComputeBond.hpp ComputeAngle.hpp ComputeLJ.hpp ComputeES.hpp Energy.hpp Integrator.hpp Output.hpp

$(TARGET): $(OBJS) $(HEADS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

depend:
	$(CXX) $(CXXFLAGS) -MM -MG $(SRCS) > Makefile.depend

clean:
	rm -f $(OBJS) $(TARGET)

-include Makefile.depend

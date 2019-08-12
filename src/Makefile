TARGET=md_alad_ewald.x

CXX=icpc

CXXFLAGS= -O3 -qopenmp -std=c++11 -I../include -I../include/uiuc -DHAVE_FFTWD=1

LIBS= -mkl

SRCS = Random.cpp dcdplugin.cpp Atom.cpp NAMDBin.cpp PDBUtil.cpp PDB.cpp PSF.cpp LoadParm.cpp Option.cpp System.cpp ComputeSpline.cpp ComputeBond.cpp ComputeAngle.cpp ComputeDihedral.cpp ComputeImproper.cpp ComputeCmap.cpp ComputeLJ.cpp ComputeES.cpp Energy.cpp Integrator.cpp Output.cpp md_alad_ewald.cpp

OBJS = $(SRCS:.cpp=.o)

HEADS=

$(TARGET): $(OBJS) $(HEADS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

depend:
	$(CXX) $(CXXFLAGS) -MM -MG $(SRCS) > Makefile.depend

clean:
	rm -f $(OBJS) $(TARGET)

-include Makefile.depend

CXXFLAGS =	-O2 -g -Wall -fmessage-length=0 -I /usr/local/include/eigen3 -I ./include

OBJS =		polynomefit.o variable_selection.o polynomefit-error-propagation.o

LIBS = 

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

polynomefit.o : polynomefit.cpp include/monopoly.h
	$(CXX) -o polynomefit.o $(CXXFLAGS) -c polynomefit.cpp

polynomefit: polynomefit.o monopoly.o
	$(CXX) -o polynomefit polynomefit.o monopoly.o
	
variable_selection.o : variable_selection.cpp include/tensor_serie.hh include/monopoly.h
	$(CXX) -o variable_selection.o $(CXXFLAGS) -c variable_selection.cpp

variable_selection: variable_selection.o correlation_tensor_serie.o tensor_serie.o monopoly.o
	$(CXX) -o variable_selection variable_selection.o correlation_tensor_serie.o tensor_serie.o monopoly.o
	
tensor_serie.o: src/tensor_serie.cc include/tensor_serie.hh include/monopoly.h
	$(CXX) -o tensor_serie.o $(CXXFLAGS) -c src/tensor_serie.cc
	
monopoly.o: src/monopoly.cc include/monopoly.h
	$(CXX) -o monopoly.o $(CXXFLAGS) -c src/monopoly.cc

correlation_tensor_serie.o: src/correlation_tensor_serie.cc include/correlation_tensor_serie.hh include/monopoly.h
	$(CXX) -o correlation_tensor_serie.o $(CXXFLAGS) -c src/correlation_tensor_serie.cc
	
full_correlation_tensor_serie.o: src/full_correlation_tensor_serie.cc include/full_correlation_tensor_serie.hh include/tensor_serie.hh include/monopoly.h
	$(CXX) -o full_correlation_tensor_serie.o $(CXXFLAGS) -c src/full_correlation_tensor_serie.cc

polynomefit-error-propagation.o: polynomefit-error-propagation.cpp include/monopoly.h include/full_correlation_tensor_serie.hh
	$(CXX) -o polynomefit-error-propagation.o $(CXXFLAGS) -c polynomefit-error-propagation.cpp
	
polynomefit-error-propagation: polynomefit-error-propagation.o monopoly.o full_correlation_tensor_serie.o tensor_serie.o
	$(CXX) -o polynomefit-error-propagation full_correlation_tensor_serie.o tensor_serie.o polynomefit-error-propagation.o monopoly.o 

treepoly_smooth.o: include/treepoly_smooth.hh src/treepoly_smooth.cc include/tensor_serie.hh
	$(CXX) -o treepoly_smooth.o $(CXXFLAGS) -c src/treepoly_smooth.cc
	
smoothed-tree-polynomial.o: smoothed-tree-polynomial.cpp include/treepoly_smooth.hh
	$(CXX) -o smoothed-tree-polynomial.o $(CXXFLAGS) -c smoothed-tree-polynomial.cpp

smoothed-tree-polynomial: smoothed-tree-polynomial.o tensor_serie.o treepoly_smooth.o monopoly.o
	$(CXX) -o smoothed-tree-polynomial smoothed-tree-polynomial.o tensor_serie.o treepoly_smooth.o monopoly.o
	
tree_regression.o: include/tree_regression.hh src/tree_regression.cc include/datautil.hh include/full_correlation_tensor_serie.hh include/tensor_serie.hh
	$(CXX) -o tree_regression.o $(CXXFLAGS) -c src/tree_regression.cc

tree-polynomial.o: tree-polynomial.cpp include/tree_regression.hh
	$(CXX) -o tree-polynomial.o $(CXXFLAGS) -c tree-polynomial.cpp

tree-polynomial: tree-polynomial.o tree_regression.o tensor_serie.o monopoly.o full_correlation_tensor_serie.o
	$(CXX) -o tree-polynomial tree-polynomial.o tree_regression.o tensor_serie.o monopoly.o full_correlation_tensor_serie.o
	
#fractional_fit.o: src/fractional_fit.cc include/fractional_fit.hh
#	$(CXX) -o fractional_fit.o $(CXXFLAGS) -c include/fracional_fit.hh

all: polynomefit variable_selection polynomefit-error-propagation smoothed-tree-polynomial tree-polynomial

clean:
	rm -f $(OBJS) $(TARGET)

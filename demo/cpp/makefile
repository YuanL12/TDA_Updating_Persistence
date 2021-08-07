CXXFLAGS = -Wall -std=c++17 -Wextra
OUTFLAGS = -O3 -fopenmp
PROFFLAGS = -pg -O3
DBFLAGS = -g -O1
OPFLAGS = -O3 -DBATS_OPCOUNT
CXX = g++
IFLAGS = -I../BATS/include -I../include

default: hello.out

%.out: %.cpp
	$(CXX) -o $@ $< $(IFLAGS) $(CXXFLAGS) $(OUTFLAGS)

%.op: %.cpp
	$(CXX) -o $@ $< $(IFLAGS) $(CXXFLAGS) $(OPFLAGS)

%.prof: %.cpp
	$(CXX) -o $@ $< $(IFLAGS) $(CXXFLAGS) $(PROFFLAGS)

%.db: %.cpp
	$(CXX) -o $@ $< $(IFLAGS) $(CXXFLAGS) $(DBFLAGS)

%: %.out %.prof %.db

.PHONY: clean
clean:
	rm -f *.out *.prof *.db
	rm -rf *.dgm

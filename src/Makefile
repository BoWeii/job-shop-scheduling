CXX=clang++
CXXFLAGS=-O3 -Wall -m64 -shared -std=c++11 -fPIC
PYBIND=`python3 -m pybind11 --includes`
PYCONFIG=`python3-config --extension-suffix`


all:
	${CXX} $(CXXFLAGS) $(PYBIND) _jss.cpp -o _jss$(PYCONFIG) -I/usr/include/python3.8
test: all
	python3 -m pytest -v -s

clean:
	rm -rf *.so __pycache__ .pytest_cache performance.txt

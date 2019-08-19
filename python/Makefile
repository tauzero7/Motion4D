
CXX   = g++
#SWIG  = /usr/local/swig/3.0.6/bin/swig
#SWIG = /usr/bin/swig
SWIG = /usr/local/misc/swig/3.0.12/bin/swig

#INC_PATH = /usr/include/python2.7
INC_PATH = /usr/local/python/3.7.3/include/python3.7m
LIB_PATH = /usr/local/python/3.7.3/lib
CXXFLAGS = -std=c++11 -O2 -fPIC

BASEFILE = m4d

M4D_LIB_DIR = ../lib
M4D_SRC_DIR  = ../src
M4D_EXTRA_DIR = $(M4D_SRC_DIR)/extra
M4D_MATH_DIR = $(M4D_SRC_DIR)/math
M4D_METRIC_DIR = $(M4D_SRC_DIR)/metric
M4D_MOTION_DIR = $(M4D_SRC_DIR)/motion
M4D_LIB_SO = libm4dd.so


BUILD_DIR = compiled/


OBJS = $(BUILD_DIR)$(BASEFILE)_wrap.o 

INCLUDES = -I$(M4D_SRC_DIR) -I$(M4D_EXTRA_DIR) \
		   -I$(M4D_MATH_DIR) -I$(M4D_METRIC_DIR) -I$(M4D_MOTION_DIR)


.PHONY:  show clean

all: _$(BASEFILE).so


_$(BASEFILE).so: $(OBJS) 
	$(CXX) -shared -o $@ $(OBJS) -Wl,-rpath $(M4D_LIB_DIR) $(M4D_LIB_DIR)/$(M4D_LIB_SO) -Wl,-rpath ${LIB_PATH}
	cp $(BUILD_DIR)$(BASEFILE).py .


$(BUILD_DIR)$(BASEFILE)_wrap.o: $(BUILD_DIR)$(BASEFILE)_wrap.cxx $(BASEFILE).h
	$(CXX) $(CXXFLAGS) -c $< -o $@ -I$(INC_PATH) $(INCLUDES)

$(BUILD_DIR)$(BASEFILE)_wrap.cxx: $(BASEFILE).i 
	@echo do swig
	$(SWIG) -o $@  -python -c++ -shadow $< 

show:
	@echo $(BASEFILE)
	@echo $(OBJ_METRIC_FILES)

clean:
	rm -rf $(BASEFILE).o $(BUILD_DIR)$(BASEFILE)_wrap.o $(BUILD_DIR)$(BASEFILE)_wrap.cxx $(BASEFILE).py $(BUILD_DIR)/$(BASEFILE).py _$(BASEFILE).so $(BUILD_DIR)/*.o



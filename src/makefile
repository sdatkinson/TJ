#=======================================================================
#
#  Makefile for TJ SLP Program
#
#=======================================================================



DIM=3

EXECUTABLE=../tj_$(DIM)d



#=======================================================================
#  Flags



DEBUG_FLAG:=-g
#OPT_LEVEL=-O0
OPT_LEVEL=-O3

DIM_FLAG:=-D DIM=$(DIM)


FLAGS=$(DEBUG_FLAG) $(OPT_LEVEL) $(DIM_FLAG)



#=======================================================================
#  Files to build



OBJS = \
./rng_class_$(DIM)d.o \
./main_$(DIM)d.o \
./input_class_$(DIM)d.o \
./parameters_class_$(DIM)d.o \
./packing_class_$(DIM)d.o \
./sphere_class_$(DIM)d.o \
./slp_class_$(DIM)d.o \
./lp_class_$(DIM)d.o \
./glpk_solver_class_$(DIM)d.o \



#=======================================================================
#  Libraries



GSL_LIB    = -L/usr/local/lib -lgsl -lgslcblas 
GLPK_LIB   = -L/usr/local/lib -lglpk
#GSL_LIB    = -L/usr/local/lib -lgsl -lgslcblas -l:libglpk.so.0.5.0
OTHER_LIBS = -lm

LIBS := $(GSL_LIB) $(GLPK_LIB) $(OTHER_LIBS)



#=======================================================================
#  Targets



all: TJ

# Tool invocations
TJ: $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ -o $(EXECUTABLE) $(OBJS) -I$(INC) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Object files	
%_$(DIM)d.o: ./%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ $(FLAGS) -std=c++11 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

# Other Targets
RM := rm -rf
clean:
	-$(RM) $(OBJS) $(EXECUTABLE)
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

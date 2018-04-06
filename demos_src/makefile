
# get options which hold for all files of the project
D=..
include $D/make.global_options


# create lists of file names
c-files_demo1 = demo1.cpp
c-files_demo2 = demo2.cpp
c-files_demo3 = demo3.cpp
c-files_demo4 = demo4.cpp
c-files_demo5 = demo5.cpp

c-files_read = read_config.cpp

o-files-main = $(addprefix ../$(OBJDIR), main.o)
o-files-read = $(addprefix ../$(OBJDIR), read_config.o)

o-files = $(addprefix ../$(OBJDIR), main.o read_config.o)
h-files      = $(sort $(shell echo *.h))

d-files      = $(o-files:.o=.d)



# top level targets
demo1_: 
	@echo "compiling: demo1 --> $(c-files_demo1)"
	@$(CXX) $(CXXFLAGS) -fopenmp $(c-files_demo1) $(OBJ) $(o-files-main)
	@echo "compiling: demo1 --> $(c-files_read)"
	@$(CXX) $(CXXFLAGS) $(c-files_read) $(OBJ) $(o-files-read)


demo2_: 
	@echo "compiling: demo2 --> $(c-files_demo2)"
	@$(CXX) $(CXXFLAGS) $(c-files_demo2) $(OBJ) $(o-files-main)
	@echo "compiling: demo2 --> $(c-files_read)"
	@$(CXX) $(CXXFLAGS) $(c-files_read) $(OBJ) $(o-files-read)


demo3_: 
	@echo "compiling: demo3 --> $(c-files_demo3)"
	@$(CXX) $(CXXFLAGS) $(c-files_demo3) $(OBJ) $(o-files-main)
	@echo "compiling: demo3 --> $(c-files_read)"
	@$(CXX) $(CXXFLAGS) $(c-files_read) $(OBJ) $(o-files-read)


demo4_: 
	@echo "compiling: demo4 --> $(c-files_demo4)"
	@$(CXX) $(CXXFLAGS) $(c-files_demo4) $(OBJ) $(o-files-main)
	@echo "compiling: demo4 --> $(c-files_read)"
	@$(CXX) $(CXXFLAGS) $(c-files_read) $(OBJ) $(o-files-read)


demo5_: 
	@echo "compiling: demo5 --> $(c-files_demo5)"
	@$(CXX) $(CXXFLAGS) $(c-files_demo5) $(OBJ) $(o-files-main)
	@echo "compiling: demo5 --> $(c-files_read)"
	@$(CXX) $(CXXFLAGS) $(c-files_read) $(OBJ) $(o-files-read)

clean: 
	rm -f $(o-files) $(d-files)
	
-include $(d-files)

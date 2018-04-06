
.DEFAULT_GOAL := demo1

include make.global_options


model_demo1: 
	cd $/demos_src && $(MAKE) $(MAKEOPTIONS) demo1_

model_demo2: 
	cd $/demos_src && $(MAKE) $(MAKEOPTIONS) demo2_
	
model_demo3: 
	cd $/demos_src && $(MAKE) $(MAKEOPTIONS) demo3_
	
model_demo4: 
	cd $/demos_src && $(MAKE) $(MAKEOPTIONS) demo4_
	
model_demo5: 
	cd $/demos_src && $(MAKE) $(MAKEOPTIONS) demo5_
	
fastchem_: 
	cd $/fastchem_src && $(MAKE) $(MAKEOPTIONS) fastchem_
	
	
all-fastchem-c: fastchem_

	@echo "fastchem compiling done."

all-l:
	cd $/obj && $(MAKE) $(MAKEOPTIONS) all-l
	
all-l-omp:
	cd $/obj && $(MAKE) $(MAKEOPTIONS) all-l-omp
	
	@echo "linking done."
	
demo1: model_demo1 all-fastchem-c all-l-omp

demo2: model_demo2 all-fastchem-c all-l

demo3: model_demo3 all-fastchem-c all-l

demo4: model_demo4 all-fastchem-c all-l

demo5: model_demo5 all-fastchem-c all-l

	@echo "everything is done and fine. enjoy your day!"

	
#and here we clean the mess
clean: clean-binary clean-model_main clean-chemistry

clean-binary: 
	rm -f $(EXECUTABLE_NAME)
	
clean-model_main: 
	cd $/demos_src && $(MAKE) clean
	
clean-chemistry: 
	cd $/fastchem_src && $(MAKE) clean
	
	@echo "all clean. have a nice day!"
	

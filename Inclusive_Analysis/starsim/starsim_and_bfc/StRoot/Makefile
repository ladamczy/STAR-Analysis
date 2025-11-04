StarGeneratorFilt.so: StrangeHadronFilter_dict.cxx StrangeHadronFilter.cxx $(STAR)/StRoot/StarGenerator/FILT/StarFilterMaker.cxx $(STAR)/StRoot/StChain/StMaker.cxx
	g++ -shared -o $@ `root-config --cflags` -I$(ROOTSYS)/include -I$(STAR)/StRoot -I$(STAR)/StRoot/StChain -I$(STAR)/StRoot/St_base -I$(STAR)/StRoot/Star2Root -I$(STAR_LIB)/../include $^ `root-config --libs`

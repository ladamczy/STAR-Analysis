#!/bin/bash

starver pro

rootcint -f StrangeHadronFilter_dict.cxx -c -I$STAR/StRoot -I$STAR/StRoot/StChain -I$STAR/StRoot/St_base -I$STAR/StRoot/Star2Root *Filter.h StarGenerator/FILT/StarFilterMaker.h StMaker.h StrangeHadronFilter_linkdef.h

gmake -f Makefile

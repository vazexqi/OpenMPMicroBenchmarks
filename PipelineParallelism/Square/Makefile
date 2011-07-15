include ../../make.def

EXEs=PipelineSerial.out PipelineTaskUnordered.out PipelineTaskOrdered.out

all: $(EXEs)

%.$(OBJ): %.cpp
	$(CPP) $(CFLAGS) $(LDFLAGS) -c $<

%.out: %.$(OBJ)
	$(CPP) $(CFLAGS) $(LDFLAGS) -o $@ $<

run: all
	time ./PipelineSerial.out 
	time ./PipelineTaskUnordered.out 
	time ./PipelineTaskOrdered.out

clean:
	$(RM) -f $(EXEs) *.$(OBJ)

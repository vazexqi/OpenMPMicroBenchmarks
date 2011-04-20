include ../../make.def

EXEs=blackscholes_for.out blackscholes_task.out

all: $(EXEs)

%.$(OBJ): %.c
	$(CC) $(CFLAGS) $(LDFLAGS) -c $<

%.out: %.$(OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

run: all
	./blackscholes_for.out 4 in_260K.txt prices.txt
	./blackscholes_task.out 4 in_260K.txt prices.txt
clean:
	$(RM) -f $(EXEs) *.$(OBJ)

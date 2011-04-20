include ../../make.def

EXEs=mmmult_for.out mmmult_task.out

all: $(EXEs)

%.$(OBJ): %.c
	$(CC) $(CFLAGS) $(LDFLAGS) -c $<

%.out: %.$(OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

run: all
	./mmmult_task.out
	./mmmult_for.out

clean:
	$(RM) -f $(EXEs) *.$(OBJ)

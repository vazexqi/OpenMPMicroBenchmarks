include ../../make.def

EXEs=mmmult_for.out mmmult_task.out

all: $(EXEs)

%.$(OBJ): %.c
	$(CC) $(CFLAGS) $(LDFLAGS) -c $<

%.out: %.$(OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

run: all
	./mmmult_for.out
	./mmmult_task.out

clean:
	$(RM) -f $(EXEs) *.$(OBJ)

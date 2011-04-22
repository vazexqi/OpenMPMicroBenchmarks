include ../../make.def

EXEs=scanner_for.out scanner_task.out

all: $(EXEs)

%.$(OBJ): %.c
	$(CC) $(CFLAGS) $(LDFLAGS) -c $<

%.out: %.$(OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

run: all
	./scanner_for.out -scanfile c756hel.in -trainfile1 a10.img -trainfile2 hc.img -stride 1 -startx 110 -starty 220 -endx 172 -endy 260 -objects 10
	./scanner_task.out -scanfile c756hel.in -trainfile1 a10.img -trainfile2 hc.img -stride 1 -startx 110 -starty 220 -endx 172 -endy 260 -objects 10

clean:
	$(RM) -f $(EXEs) *.$(OBJ)

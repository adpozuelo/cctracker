TEST = centers
MODS = common.cuf distances.cuf trajectories.cuf io.cuf
FC = nvfortran
EXE = exe
FCFLAGS = -fast

all: build run verify

build: $(TEST).cuf
	$(FC) $(FCFLAGS) $(MODS) -o $(TEST).$(EXE) $<

run: $(TEST).$(EXE)
	$(RUN) ./$(TEST).$(EXE)

verify:

clean:
	@echo 'Cleaning up...'
	@rm -rf *.$(EXE) *.dwf *.pdb *.mod *.o prof

# Helper read together with SNANA's configured src/Makefile, so this folder
# compiles/links with exactly SNANA's compiler and libraries:
#   make -s -C $SNANA_DIR/src -f Makefile -f <this file> snana-libs
snana-cc:
	@echo $(CC)
snana-libs:
	@echo $(SNLDFLAGS) $(LGSL) $(LCFITSIO) $(LROOT) $(CPPLIB) -lm
.PHONY: snana-cc snana-libs

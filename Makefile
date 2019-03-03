# Makefile - prepared for Pyrat-Bay package
#
# `make` - Build and compile the pyratbay package.
# `make clean` - Remove all compiled (non-source) files that are created.
#
# If you are interested in the commands being run by this makefile, you may add
# "VERBOSE=1" to the end of any `make` command, i.e.:
#
# 		make VERBOSE=1
#
# This will display the exact commands being used for building, etc.

# Set verbosity
#
Q = @
O = > /dev/null
ifdef VERBOSE
	ifeq ("$(origin VERBOSE)", "command line")
		Q =
		O =
	endif
else
	MAKEFLAGS += --no-print-directory
endif

LIBDIR = pyratbay/lib/

# Get the location of this Makefile.
mkfile_dir := $(dir $(lastword $(MAKEFILE_LIST)))

# `make [clean]` should run `make [clean]` on all of the modules.
all: make_pb make_mc3 make_pytips
clean: clean_pb clean_mc3 clean_pytips


make_pb:
	@echo "Building Pyrat-Bay package."
	$(Q) python setup.py build $(O)
	@mv -f build/lib.*/pyratbay/lib/*.so $(LIBDIR)
	@rm -rf build/
	@echo "Successful compilation."
	@echo ""

make_mc3:
	@cd $(mkfile_dir)/modules/MCcubed/ && make

make_pytips:
	@cd $(mkfile_dir)/modules/pytips/ && make


clean_pb:
	@rm -rf $(LIBDIR)*.so
	@echo "Cleaned Pyrat Bay."

clean_mc3:
	@cd $(mkfile_dir)/modules/MCcubed && make clean
	@echo "Cleaned MC3."

clean_pytips:
	@cd $(mkfile_dir)/modules/pytips/ && make clean
	@echo "Cleaned pytips."


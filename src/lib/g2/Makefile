include version.inc
build:
	mkdir -p ../exec ../intel/include/g2_$(VERSION_STRING)_4 ../intel/include/g2_$(VERSION_STRING)_d && \
	( set -xue ; make -f makefile_4_iobuf_intel ) && \
	( set -xue ; make -f makefile_d_iobuf_intel ) && \
	( set -xue ; cd ../unit_test ; make all )

clean:
	rm -rf exec ../intel/include/g2_$(VERSION_STRING)_4 ../intel/include/g2_$(VERSION_STRING)_d *.a && \
	( set -xue ; cd ../unit_test ; make bare )

test: build
	../exec/g2_unit_test_4 && ../exec/g2_unit_test_d

# Common synonyms:
bare: clean
distclean: clean
all: build test

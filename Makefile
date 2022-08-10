CMAKE_COMPILER_OPTS=
CMAKE_TESTS_OPTS=-DENABLE_TESTS=on -DENABLE_FPGA=on
CMAKE_DTESTS_OPTS=${CMAKE_COMPILER_OPTS} -DCMAKE_BUILD_TYPE=debug ${CMAKE_TESTS_OPTS}
CMAKE_OTESTS_OPTS=${CMAKE_COMPILER_OPTS} -DCMAKE_BUILD_TYPE=optim ${CMAKE_TESTS_OPTS}

all: build
	make -C build

install: build
	make -C build install

clean: build
	make -C build clean

distclean:
	rm -rf build builddtests buildotests *.log _x .Xil *.compile_summary *.xo* *.info *.link_summary *.csv *.run_summary .run/ fft.xclbin.*; true

# test: builddtests src/test/googletest/CMakeLists.txt
test: builddtests buildotests src/test/googletest/CMakeLists.txt
	make -C builddtests VERBOSE=1
	make -C buildotests VERBOSE=1
# make -j $(nproc) -C builddtests VERBOSE=1
# make -j $(nproc) -C buildotests VERBOSE=1
#	make -j $(nproc) -C builddtests test VERBOSE=1
#	make -j $(nproc) -C buildotests test VERBOSE=1

build: src/test/googletest/CMakeLists.txt
	mkdir build; cd build; cmake ../src -DCMAKE_INSTALL_PREFIX=/home/ejaco020/tfhe/; cd ..

builddtests:
	rm -rf $@; true; mkdir $@;
	cd $@; cmake ../src ${CMAKE_DTESTS_OPTS};
	cd $@; cmake ../src ${CMAKE_DTESTS_OPTS};
	cd ..

buildotests:
	rm -rf $@; true; mkdir $@;
	cd $@; cmake ../src ${CMAKE_OTESTS_OPTS};
	cd $@; cmake ../src ${CMAKE_OTESTS_OPTS};
	cd ..

src/test/googletest/CMakeLists.txt:
	git submodule init
	git submodule update

alltests:
	make distclean && make test CMAKE_COMPILER_OPTS="-DCMAKE_CXX_COMPILER=clang++-libc++ -DCMAKE_C_COMPILER=clang"
	make distclean && make test CMAKE_COMPILER_OPTS="-DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_COMPILER=clang"
	make distclean && make test CMAKE_COMPILER_OPTS="-DCMAKE_CXX_COMPILER=g++-9 -DCMAKE_C_COMPILER=gcc-9"
	make distclean && make test CMAKE_COMPILER_OPTS="-DCMAKE_CXX_COMPILER=g++-8 -DCMAKE_C_COMPILER=gcc-8"
	make distclean && make test CMAKE_COMPILER_OPTS="-DCMAKE_CXX_COMPILER=g++-7 -DCMAKE_C_COMPILER=gcc-7"

VPP := v++
PLATFORM := xilinx_u280_xdma_201920_3
TARGET := sw_emu
CONFIG_NAME := config.cfg
KERNEL_XO := tGswFFTExternMulToTLwe.xo
KERNEL_SOURCES := fft_transform_reverse.cpp fft_transform.cpp TorusPolynomial_fft.cpp tGswTorus32PolynomialDecompH.cpp IntPolynomial_ifft.cpp tLweFFTClear.cpp tLweFFTAddMulRTo.cpp tLweFromFFTConvert.cpp TorusPolynomial_ifft.cpp
KERNEL_FOLDER := ./src/kernels
PROJECT_NAME := fft

# ifeq ($(TARGET), sw_emu) # sw_emu needs the sources, easiest to generate their xo files to include the sources. hw(_emu) have sources manually included
# 	KERNEL_XO += fft_transform_reverse.xo fft_transform.xo
# endif

VPP_XCLBIN_FLAGS := -l -j $(shell nproc) --profile_kernel data:all:all:all -O1 --platform $(PLATFORM) -t $(TARGET) --input_files $(KERNEL_XO) -o $(PROJECT_NAME).xclbin
VPP_XO_FLAGS := -c -O1 --platform $(PLATFORM) -t $(TARGET) -I$(KERNEL_FOLDER)/include/

xclbin: $(KERNEL_XO)
	$(VPP) $(VPP_XCLBIN_FLAGS)
	emconfigutil --platform $(PLATFORM) --nd 1

%.xo: src/kernels/%.cpp
# ifeq ($(TARGET), sw_emu)
# 	$(VPP) $(VPP_XO_FLAGS) -k $(basename $(notdir $<)) $< -o $@
# else
	$(VPP) $(VPP_XO_FLAGS) -k $(basename $(notdir $<)) $< $(addprefix $(KERNEL_FOLDER)/, $(KERNEL_SOURCES)) -o $@
# endif

runtest: test xclbin
	cp $(PROJECT_NAME).xclbin emconfig.json ./builddtests/test
	XCL_EMULATION_MODE=$(TARGET) ./builddtests/test/test-bootstrapping-fft-fpga
	XCL_EMULATION_MODE=$(TARGET) ./builddtests/test/test-gate-bootstrapping-fpga
# XCL_EMULATION_MODE=$(TARGET) ./builddtests/test/test-decomp-tgsw-fpga
	XCL_EMULATION_MODE=$(TARGET) ./builddtests/test/test-lwe-fpga
	XCL_EMULATION_MODE=$(TARGET) ./builddtests/test/test-multiplication-fpga
	XCL_EMULATION_MODE=$(TARGET) ./builddtests/test/test-addition-boot-fpga
	XCL_EMULATION_MODE=$(TARGET) ./builddtests/test/test-long-run-fpga
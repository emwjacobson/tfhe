CMAKE_COMPILER_OPTS=-DCMAKE_INSTALL_PREFIX=/home/ejaco020/tfhe/ -DCMAKE_C_COMPILER='/tools/Xilinx/Vitis/2020.1/bin/clang_wrapper' -DCMAKE_CXX_COMPILER='/tools/Xilinx/Vitis/2020.1/bin/clang_wrapper'
CMAKE_TESTS_OPTS=-DENABLE_TESTS=on -DENABLE_FFTW=off \
		 -DENABLE_NAYUKI_PORTABLE=on -DENABLE_NAYUKI_AVX=on \
		 -DENABLE_SPQLIOS_AVX=on -DENABLE_SPQLIOS_FMA=on -DENABLE_FPGA_FFT=on
CMAKE_DTESTS_OPTS=${CMAKE_COMPILER_OPTS} -DCMAKE_BUILD_TYPE=debug ${CMAKE_TESTS_OPTS}
CMAKE_OTESTS_OPTS=${CMAKE_COMPILER_OPTS} -DCMAKE_BUILD_TYPE=optim ${CMAKE_TESTS_OPTS}

all: build
	make -C build

install: build
	make -C build install

clean: build
	make -C build clean

distclean:
	rm -rf build builddtests buildotests *.log _x .Xil *.compile_summary *.xo *.info *.link_summary; true

test: builddtests buildotests src/test/googletest/CMakeLists.txt
	make -j $(nproc) -C builddtests VERBOSE=1
	make -j $(nproc) -C buildotests VERBOSE=1
	make -j $(nproc) -C builddtests test VERBOSE=1
	make -j $(nproc) -C buildotests test VERBOSE=1

build: src/test/googletest/CMakeLists.txt
	mkdir build; cd build; cmake ../src ${CMAKE_COMPILER_OPTS} -DENABLE_TESTS=on -DENABLE_FFTW=off -DENABLE_NAYUKI_PORTABLE=on -DENABLE_NAYUKI_AVX=on -DENABLE_SPQLIOS_AVX=on -DENABLE_SPQLIOS_FMA=on -DENABLE_FPGA_FFT=on; cd ..

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
# CONFIG_NAME := config.cfg
KERNEL_XO := fft_transform_reverse.xo
PROJECT_NAME := fft

# VPP_XCLBIN_FLAGS := -l --profile_kernel data:all:all:all -O1 --platform $(PLATFORM) -t $(TARGET) --config $(CONFIG_NAME) $(KERNEL_XO) -o $(PROJECT_NAME).xclbin
VPP_XCLBIN_FLAGS := -l --profile_kernel data:all:all:all -O1 --platform $(PLATFORM) -t $(TARGET) $(KERNEL_XO) -o $(PROJECT_NAME).xclbin
VPP_XO_FLAGS := -c --platform $(PLATFORM) -I src/libtfhe/fft_processors/fpga/

xclbin: $(KERNEL_XO)
	$(VPP) $(VPP_XCLBIN_FLAGS)

%.xo: src/libtfhe/fft_processors/fpga/kernels/%.cpp
	$(VPP) $(VPP_XO_FLAGS) -k $(basename $(notdir $<)) $< -o $@

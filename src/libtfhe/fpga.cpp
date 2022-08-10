#include <iostream>
#include <unistd.h>
#include <fstream>
#include "fpga.h"
#include "tfhe_core.h"
#include "polynomials.h"

FPGA_Processor fpga;

FPGA_Processor::FPGA_Processor() {
    omegaxminus1 = (cplx*) malloc(sizeof(cplx) * Value_2N);
    for (int32_t x=0; x<Value_2N; x++) {
        omegaxminus1[x]=cplx(cos(x*M_PI/Value_N)-1., sin(x*M_PI/Value_N)); // instead of cos(x*M_PI/N)-1. + sin(x*M_PI/N) * 1i
        //exp(i.x.pi/N)-1
    }

    // Initialize OpenCL Environment
    cl_int err;
    std::vector<cl::Device> devices = get_xilinx_devices();
    devices.resize(1);
    cl::Device device = devices[0];
    context = cl::Context(device, NULL, NULL, NULL, &err);
    auto filebuf = read_binary_file("fft.xclbin");
    cl::Program::Binaries bins{{filebuf.data(), filebuf.size()}};

    program = cl::Program(context, devices, bins, NULL, &err);
    q = cl::CommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &err);
    k_tGswFFTExternMulToTLwe = cl::Kernel(program, "tGswFFTExternMulToTLwe", &err);
    printf("Finished loading FPGA kernels\n");
}

FPGA_Processor::~FPGA_Processor() {
    free(omegaxminus1);
}

std::vector<cl::Device> FPGA_Processor::get_xilinx_devices() {
    size_t i;
    cl_int err;
    std::vector<cl::Platform> platforms;
    err = cl::Platform::get(&platforms);
    cl::Platform platform;
    for (i  = 0 ; i < platforms.size(); i++){
        platform = platforms[i];
        std::string platformName = platform.getInfo<CL_PLATFORM_NAME>(&err);
        if (platformName == "Xilinx"){
            std::cout << "INFO: Found Xilinx Platform" << std::endl;
            break;
        }
    }
    if (i == platforms.size()) {
        std::cout << "ERROR: Failed to find Xilinx platform" << std::endl;
        exit(EXIT_FAILURE);
    }

    //Getting ACCELERATOR Devices and selecting 1st such device
    std::vector<cl::Device> devices;
    err = platform.getDevices(CL_DEVICE_TYPE_ACCELERATOR, &devices);
    return devices;
}

std::vector<unsigned char> FPGA_Processor::read_binary_file(const std::string &xclbin_file_name) {
    std::cout << "INFO: Reading " << xclbin_file_name << std::endl;

    if (access(xclbin_file_name.c_str(), R_OK) != 0) {
        printf("ERROR: %s xclbin not available please build\n",
               xclbin_file_name.c_str());
        exit(EXIT_FAILURE);
    }
    //Loading XCL Bin into char buffer
    std::cout << "Loading: '" << xclbin_file_name.c_str() << "'\n";
    std::ifstream bin_file(xclbin_file_name.c_str(), std::ifstream::binary);
    bin_file.seekg(0, bin_file.end);
    auto nb = bin_file.tellg();
    bin_file.seekg(0, bin_file.beg);
    std::vector<unsigned char> buf;
    buf.resize(nb);
    bin_file.read(reinterpret_cast<char*>(buf.data()), nb);
    return buf;
}
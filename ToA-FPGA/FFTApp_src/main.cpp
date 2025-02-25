#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <stdint.h>
#include "xaxidma.h"
#include "xparameters.h"
#include "xil_printf.h"
#include "xaxidma_hw.h"
#include "xil_io.h"
#include "sleep.h"
#include "xtime_l.h"

#include "fft_data.h"



typedef std::complex<float> Complex;



#define FFT_SIZE 1024
#define EPSILON 0.001


static Complex output_data[FFT_SIZE] __attribute__((aligned(64)));


// output compare to header
//void compare_results() {
//	char* buffer = (char*)malloc(10000);
//    int mismatch_count = 0;
//
//    xil_printf("Comparing FFT Output with Reference Data...\n");
//
//    for (int i = 0; i < FFT_SIZE; i++) {
//        float real_diff = fabs(output_data[i].real() - fft_output_real[i]);
//        float imag_diff = fabs(output_data[i].imag() - fft_output_imag[i]);
//
//        if (real_diff > EPSILON || imag_diff > EPSILON) {
//        	xil_printf(buffer,"real_diff: %.4f,imag_diff: %.4f", real_diff, imag_diff);
//            xil_printf("%s\n", buffer);
//        	mismatch_count++;
//        }
//    }
//
//    if (mismatch_count == 0) {
//        xil_printf("All FFT output values match within tolerance (EPSILON = %f).\n", EPSILON);
//    } else {
//        xil_printf("Total mismatches: %d / %d\n", mismatch_count, FFT_SIZE);
//    }
//}

void compare_results() {
    int mismatch_count = 0;

    xil_printf("Comparing FFT Output with Reference Data...\n");

    for (int i = 0; i < FFT_SIZE; i++) {
        float real_pl = output_data[i].real();
        float imag_pl = output_data[i].imag();
        float real_ref = fft_output_real[i];
        float imag_ref = fft_output_imag[i];

        float real_diff = fabs(real_pl - real_ref);
        float imag_diff = fabs(imag_pl - imag_ref);

        if (real_diff > EPSILON || imag_diff > EPSILON) {
            xil_printf("Mismatch at index %d!\n", i);
            mismatch_count++;
        }
    }

    if (mismatch_count == 0) {
        xil_printf("All FFT output values match within tolerance.\n");
    } else {
        xil_printf("Total mismatches: %d / %d\n", mismatch_count, FFT_SIZE);
    }
}



int init_dma(XAxiDma *dma_instance, u16 device_id) {
    XAxiDma_Config *dma_config;
    dma_config = XAxiDma_LookupConfig(device_id);
    if (!dma_config) {
        xil_printf("No DMA config found for Device ID %d.\n", device_id);
        return XST_FAILURE;
    }
    if (XAxiDma_CfgInitialize(dma_instance, dma_config) != XST_SUCCESS) {
        xil_printf("DMA initialization failed for Device ID %d.\n", device_id);
        return XST_FAILURE;
    }
    if (XAxiDma_HasSg(dma_instance)) {
        xil_printf("Device ID %d configured as SG mode. SG mode is not supported.\n", device_id);
        return XST_FAILURE;
    }
    xil_printf("DMA initialized successfully for Device ID %d.\n", device_id);
    return XST_SUCCESS;
}

//void print_complex_data(const char* stream_name, const Complex* data, int size) {
//    char buffer[10000];
//    xil_printf("Complex data from %s:\n", stream_name);
//    for (int i = 0; i < size; ++i) {
//        sprintf(buffer, "[%d]: %.4f + %.4fi", i, data[i].real(), data[i].imag());
//        xil_printf("%s\n", buffer);
//    }
//}

void print_complex_data(const char* stream_name, const Complex* data, int size) {
    char* buffer = (char*)malloc(10000);
    if (!buffer) {
        xil_printf("Memory allocation failed for buffer!\n");
        return;
    }

    xil_printf("Complex data from %s:\n", stream_name);
    for (int i = 0; i < size; ++i) {
        sprintf(buffer, "[%d]: %.4f + %.4fi", i, data[i].real(), data[i].imag());
        xil_printf("%s\n", buffer);
    }

    free(buffer);
}


int main() {
    XTime start_time, end_time1;
    XAxiDma AxiDma_FFT_in,AxiDma_FFT_out;
    int status, status1, status2;

    xil_printf("Initializing DMA for FFT...\n");
    status = init_dma(&AxiDma_FFT_in, XPAR_AXI_DMA_3_DEVICE_ID);
    if (status != XST_SUCCESS) {
        xil_printf("Failed to initialize DMA for FFT in.\n");
        return XST_FAILURE;
    }
    status = init_dma(&AxiDma_FFT_out, XPAR_AXI_DMA_2_DEVICE_ID);
    if (status != XST_SUCCESS) {
        xil_printf("Failed to initialize DMA for FFT out.\n");
        return XST_FAILURE;
    }

    XAxiDma_IntrDisable(&AxiDma_FFT_in, XAXIDMA_IRQ_ALL_MASK, XAXIDMA_DMA_TO_DEVICE);
    XAxiDma_IntrDisable(&AxiDma_FFT_out, XAXIDMA_IRQ_ALL_MASK, XAXIDMA_DEVICE_TO_DMA);

//    Complex *input_data = (Complex *)malloc(FFT_SIZE * sizeof(Complex));
//    Complex *output_data = (Complex *)malloc(FFT_SIZE * sizeof(Complex));
    static Complex input_data[FFT_SIZE] __attribute__((aligned(64)));
//    static Complex output_data[FFT_SIZE] __attribute__((aligned(64)));

    if (!input_data || !output_data) {
        xil_printf("Memory allocation failed!\n");
        return XST_FAILURE;
    }

// impulse input
//    for (int i = 0; i < FFT_SIZE; i++) {
//        if (i == 0)
//            input_data[i] = Complex(1.0f, 0.0f);
//        else
//            input_data[i] = Complex(0.0f, 0.0f);
//    }

// rect input
	#define RECT_WIDTH 16

	for (int i = 0; i < FFT_SIZE; i++) {
		if (i < RECT_WIDTH)
			input_data[i] = Complex(1.0f, 0.0f);
		else
			input_data[i] = Complex(0.0f, 0.0f);
	}

// input from header
	for (int i = 0; i < FFT_SIZE; i++) {
	    input_data[i] = Complex(fft_input_real[i], fft_input_imag[i]);
	}

    xil_printf("Sending impulse response to FFT...\n");

    XTime_GetTime(&start_time);

    Xil_DCacheFlushRange((UINTPTR)&output_data[0], sizeof(Complex) * FFT_SIZE);

    status1 = XAxiDma_SimpleTransfer(&AxiDma_FFT_out, (UINTPTR)&output_data[0], sizeof(Complex) * FFT_SIZE, XAXIDMA_DEVICE_TO_DMA);
    if (status1 != XST_SUCCESS) {
        xil_printf("DMA transfer for FFT output failed.\n");
        return XST_FAILURE;
    }


    Xil_DCacheFlushRange((UINTPTR)&input_data[0], sizeof(Complex) * FFT_SIZE);

    status2 = XAxiDma_SimpleTransfer(&AxiDma_FFT_in, (UINTPTR)&input_data[0], sizeof(Complex) * FFT_SIZE, XAXIDMA_DMA_TO_DEVICE);
    if (status2 != XST_SUCCESS) {
        xil_printf("DMA transfer for FFT input failed.\n");
        return XST_FAILURE;
    }

    while (XAxiDma_Busy(&AxiDma_FFT_in, XAXIDMA_DMA_TO_DEVICE)) {
        xil_printf("Waiting for FFT input DMA transfer to complete...\n");
    }

    u32 dma_status = XAxiDma_ReadReg(AxiDma_FFT_in.RegBase + XAXIDMA_TX_OFFSET, XAXIDMA_SR_OFFSET);
    xil_printf("DMA Input Status Register: 0x%08x\n", dma_status);


    xil_printf("Receiving FFT output...\n");
    while (XAxiDma_Busy(&AxiDma_FFT_out, XAXIDMA_DEVICE_TO_DMA)) {
        xil_printf("Waiting for FFT output DMA transfer to complete...\n");
    }

    Xil_DCacheInvalidateRange((UINTPTR)output_data, sizeof(Complex) * FFT_SIZE);
    XTime_GetTime(&end_time1);


    print_complex_data("FFT Output", output_data, FFT_SIZE);

    double elapsed_time = 1.0 * (end_time1 - start_time) / COUNTS_PER_SECOND;
    printf("Elapsed time: %f seconds\n", elapsed_time);

//    free(input_data);
//    free(output_data);

    xil_printf("FFT processing completed.\n");

//    compare_results();
    return 0;
}


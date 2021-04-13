#include <stdio.h>
#include <complex.h>
#include <math.h>

// Helper function for the function fft()
void fft_(complex* inp, int N, int isfft)
{
    // Base case
    if (N <= 1) return;

    // Storing odd and even parts in different arrays
    complex* odd = malloc(N/2 * sizeof(complex));
    complex* even = malloc(N/2 * sizeof(complex));

    for (int i = 0; i < N / 2; i++){
        even[i] = inp[2 * i];
        odd[i] = inp[2 * i + 1];
    }

    // Recursively get fft/ifft of the smaller arrays
    fft_(even, N/2, isfft);
    fft_(odd, N/2, isfft);

    // Updating original array with the transform
    for (int j = 0; j < N/2; j++)
    {
        complex w = cexp(-I*2*M_PI*j*isfft / N);
        inp[j] = even[j] + w*odd[j];
        inp[j + N/2] = even[j] - w*odd[j];
    }

    // Freeing memory no longer needed
    free(even);
    free(odd);
}

// The parameter isfft determines whether to take the fourier transform or the inverse
// fourier transform of the given array. isfft=1 for fourier transform and isfft= -1 for inverse
void fft(complex* inp, int N, int isfft){
    fft_(inp, N, isfft);

    // For inverse fft, we scale the array
    if(isfft==-1){
        for(int i=0;i<N;i++)  inp[i] = inp[i]/N;
    }
    return;
}


int main()
{
    // N should be a power of 2 for this algorithm to work
    int N = (1<<20);
    complex x[N];

    FILE *finp, *fout1, *fout2;
    // We take data from a file for the signal
    finp = fopen("x.dat", "r");
    int len = 0;
    double a;
    while(!feof(finp) && len<n){
            fscanf(finp, "%lf", &a);
            X[len] = CMPLX(a,0);
            len++;
    }

    // Storing the fourier transform data in a dat file
    fft(x, N, 1);
    fout1 = fopen("fourier_fft.dat","w");
    for(int i=0;i<N;i++)
    {
        fprintf(fout1, "%1f+%1fi\n", creal(x[i]), cimag(x[i]));
    }

    // Setting the isfft parameter to -1 gives the ifft. We store this too in a .dat file
    fft(x, N, -1);
    fout2 = fopen("Reconst_ifft.dat","w");
    for(int i=0;i<N;i++)
    {
        fprintf(fout2,"%lf \n",creal(x[i]));
    }
    fclose(finp);
    fclose(fout1);
    fclose(fout2);
    return 0;
}

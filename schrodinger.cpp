#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>

#define N 200

using namespace std;

int main()
{
    int j, t;

    double lambda, s, k, pot[N + 1], sigma, ciclos, j1, norma, norma1;
    complex<double> alpha[N + 1], fonda[N + 1], beta[N + 1], x[N + 1], unimag;
    FILE *f1, *f2;

    //Abrimos los archivos
    f1=fopen("fonda.txt","w");
    f2=fopen("norma.txt","w");

    // Damos los parámetros iniciales:
    // Los ciclos van de 1 a N/4
    ciclos=N/4;

    // Cambiaremos la altura de la barrera de potencial:
    lambda = 0.3;

    // Generamos s, k, el potencial y el valor de la función de onda.
    k = 2.0 * M_PI * double(ciclos) / double(N);
    s = 1.0 / (4 * k * k);

    for (j = 0; j <= N; j++) {
        if ((j >= 0.4 * N) && (j <= 0.6 * N))
            pot[j] = lambda * k * k;
        else
            pot[j] = 0.0;
    }

    norma=0.;
    unimag=complex<double>(0.0, 1.0);

    // Imponemos las condiciones de contorno.
    fonda[0] = complex<double>(0.0, 0.0);
    fonda[N] = complex<double>(0.0, 0.0);

    // La funcion de onda en forma general es fonda[j] = exp(complex<double>(0.0, 1.0) * k * j) * exp((-8.0 * pow((4.0 * j - N), 2)) / (N * N));

    for (j = 1; j <= N - 1; j++) {
        j1 = 1.0 * j;
        fonda[j] = exp(unimag * k * j1) * exp(-8.0 * pow((4.0 * j1 - double(N)), 2) / pow(double(N), 2));
        norma=norma + pow(imag(fonda[j]), 2) + pow(real(fonda[j]), 2);
    }

    // Normalizamos la función de onda
    //for (j=0; j<=N; j++)
        //fonda[j]=fonda[j]/pow(norma, 1./2);

    // Escribimos en el fichero la función de onda en el instante t=0
    for (j=0; j<=N; j++)
    {
        fprintf(f1, "%i%c\t%e%c\t%e%c\t%e%c\t%e\n", j, 44, real(fonda[j]), 44, imag(fonda[j]), 44, pot[j], 44, norm(fonda[j]));
    }
    fprintf(f1, "\n");
   
    // Obtenemos alpha y beta en orden decreciente.
    // Comenzamos con alpha. Como no depende del tiempo la calculamos 1 vez.

    alpha[N - 2] = complex<double>(0.0, 0.0);

    for ( j = N-2; j > 0; j--)
    {
        alpha[j-1]=-1./(-2.+2.*(unimag/s)+alpha[j]-pot[j]);
    }

    // Comenzamos la recurrencia.

    for (t = 0; t <= 1000; t++)
    {
        // Calculamos ahora beta:

        beta[N-2]=complex<double>(0.0, 0.0);

        for ( j = N-2; j > 0; j--)
          {
           beta[j-1]=(1./(-2.+2.*(unimag/s)+alpha[j]-pot[j]))*(4.*unimag*(fonda[j]/s)-beta[j]);
          }

        // Imponemos condiciones de contorno.

        x[0] = complex<double>(0.0, 0.0);
        x[N] = complex<double>(0.0, 0.0);

        for (j = 0; j <= N - 2; j++)
            x[j + 1] = alpha[j] * x[j] + beta[j];

        // Calculamos el valor de la función de onda y su norma.

        norma=0.;

        for (j = 0; j <= N; j++)
        {
            fonda[j] = x[j] - fonda[j];
            norma=norma + (pow(real(fonda[j]),2) + pow(imag(fonda[j]),2));
        }

        fonda[0] = complex<double>(0.0, 0.0);
        fonda[N] = complex<double>(0.0, 0.0);


        for (j=0; j<=N; j++)
        {
        fprintf(f1, "%i%c\t%e%c\t%e%c\t%e%c\t%e\n", j, 44, real(fonda[j]), 44, imag(fonda[j]), 44, pot[j], 44, norm(fonda[j]));
        }

        fprintf(f1, "\n");
        fprintf(f2, "%i\t%e\n", t, norma);
    }


    return 0;
}

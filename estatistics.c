#include<stdio.h>
#include "estatistics.h"

double media(double numeros[], int tamanho){
    double soma = 0.0;
    for  (int i = 0; i < tamanho; i ++){
        soma += numeros[i];
    }
    return soma/tamanho;
}

double mediana(double numeros[], int tamanho) {
    double mediana = 0.0;
    for (int i = 0; i < tamanho - 1; i++) {
        for (int j = 0; j < tamanho - i - 1; j++) {
            if (numeros[j] > numeros[j + 1]) {
                double temp = numeros[j];
                numeros[j] = numeros[j + 1];
                numeros[j + 1] = temp;
            }
        }
    }
    if (tamanho % 2 == 0) {
        mediana = (numeros[tamanho/2 - 1] + numeros[tamanho/2]) / 2.0;
    } else {
        mediana = numeros[tamanho/2];
    }
    return mediana;
}

double moda(double numeros[], int tamanho) {
    int maxCont = 0;
    double valorModa = numeros[0];

    for (int i = 0; i < tamanho; i++) {
        int cont = 0;
        for (int j = 0; j < tamanho; j++) {
            if (numeros[j] == numeros[i])
                cont++;
        }
        if (cont > maxCont) {
            maxCont = cont;
            valorModa = numeros[i];
        }
    }

    return valorModa;
}


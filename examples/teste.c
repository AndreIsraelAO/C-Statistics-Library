#include <stdio.h>
#include "../include/estatistics.h"

int main() {
    double numeros[] = {1, 2, 2, 3, 4};
    int tamanho = 5;

    printf("Média: %.2f\n", media(numeros, tamanho));
    printf("Mediana: %.2f\n", mediana(numeros, tamanho));
    printf("Moda: %.2f\n", moda(numeros, tamanho));

    return 0;
}

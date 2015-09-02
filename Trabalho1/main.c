#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "matEntrada.h"
#include "interfaces.h"
#define MAX_DOUBLE 2

int main (int argc, char* argv[]){

    // Variáveis auxiliares para leitura do arquivo
    int i, j;
    int auxInt;
    short auxShort;
    double auxDouble;
    long auxLong;
    char* caminhoMat = telaInicial(-1);
    double* resultado = NULL;
    double* x0 = NULL;

    // Ponteiro para o arquivo que será aberto
    FILE* arqMat = NULL;
    // Cria a estrutura que armazenará as matrizes
    MAT_ENTRADA* m = criarMatEntrada();

    arqMat = fopen(caminhoMat,"r");

    if(arqMat == NULL){
        printf("Erro ao abrir arquivo.");
    }
    else{
        // Le as primeiras informações do arquivo
        fscanf(arqMat, "%d", &auxInt);
        fscanf(arqMat, "%hd", &auxShort);
        fscanf(arqMat, "%lf", &auxDouble);
        fscanf(arqMat, "%ld", &auxLong);
        // Salva os valores na estrutura
        inicializarValsMatEntrada(m,auxInt,auxShort,auxDouble,auxLong,0);
        // Aloca espaço para as matrizes A e B e diagonalAux, agora que sabemos sua dimensão
        alocarMatA(m);
        alocarMatB(m);
        alocarDiagonalAux(m);
        // Insere os elementos na Matriz A
        for(i=0;i<getOrdem(m);i++){
            for(j=0;j<getOrdem(m);j++){
                fscanf(arqMat,"%lf",&auxDouble);
                inserirElemMatA(m,auxDouble,i,j);
            }
        }
        // Insere os elementos na Matriz B
        for(i=0;i<getOrdem(m);i++){
            fscanf(arqMat,"%lf",&auxDouble);
            inserirElemMatB(m,auxDouble,i);
        }
        // Fecha o arquivo agora que não precisamos mais dele
        fclose(arqMat);

        //imprimirInfosMatEntrada(m);
    }

    x0 = (double*)malloc(getOrdem(m)*sizeof(double));

    srand((unsigned int)time(NULL));
    for(i=0;i<getOrdem(m);i++){
        x0[i] = ((float)rand()/(float)(RAND_MAX))*MAX_DOUBLE;
    }

    prepararMatrizes(m);
    resultado = jacobiRichardson(m,x0);
    imprimirResultado(m,resultado);

    // Libera a memória alocada pela estrutura
    destruirMatEntrada(m);

    return 0;
}

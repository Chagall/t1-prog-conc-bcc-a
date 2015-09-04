#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "interfaces.h"

char* getCaminhoArquivo(short opcao){

    // Caminhos até cada arquivo das matrizes
    const char caminhoMat250[] = "matrizes/matriz250.txt";
    const char caminhoMat500[] = "matrizes/matriz500.txt";
    const char caminhoMat1000[] = "matrizes/matriz1000.txt";
    const char caminhoMat1500[] = "matrizes/matriz1500.txt";
    const char caminhoMat2000[] = "matrizes/matriz2000.txt";
    const char caminhoMat3000[] = "matrizes/matriz3000.txt";
    const char caminhoMat4000[] = "matrizes/matriz4000.txt";
    char* caminho = NULL;

    switch(opcao){
        case 1:
            caminho = (char*)malloc((strlen(caminhoMat250))*sizeof(char));
            strcpy(caminho,caminhoMat250);
            break;
        case 2:
            caminho = (char*)malloc((strlen(caminhoMat500))*sizeof(char));
            strcpy(caminho,caminhoMat500);
            break;
        case 3:
            caminho = (char*)malloc((strlen(caminhoMat1000))*sizeof(char));
            strcpy(caminho,caminhoMat1000);
            break;
        case 4:
            caminho = (char*)malloc((strlen(caminhoMat1500))*sizeof(char));
            strcpy(caminho,caminhoMat1500);
            break;
        case 5:
            caminho = (char*)malloc((strlen(caminhoMat2000))*sizeof(char));
            strcpy(caminho,caminhoMat2000);
            break;
        case 6:
            caminho = (char*)malloc((strlen(caminhoMat3000))*sizeof(char));
            strcpy(caminho,caminhoMat3000);
            break;
        case 7:
            caminho = (char*)malloc((strlen(caminhoMat4000))*sizeof(char));
            strcpy(caminho,caminhoMat4000);
            break;
    }

    return(caminho);
}

char* telaInicial(short opInicial, short* opEscolhida){

    short opcao = opInicial;
    char* caminho = NULL;

    if(opcao == -1){
        printf("\t\tSSC0143 - Programacao Concorrente - Turma A\n");
        printf("\t\t Trabalho 01 - Metodo de Jacobi-Richardson\n\n");
        printf("\tAlunos: Fabio Alves Martins Pereira (Num.USP: 7987435)\n\n");
        printf("\t\t\t\tMenu Principal\n\n");
        printf("\t1 - Carregar Matriz de 250 elementos\n");
        printf("\t2 - Carregar Matriz de 500 elementos\n");
        printf("\t3 - Carregar Matriz de 1000 elementos\n");
        printf("\t4 - Carregar Matriz de 1500 elementos\n");
        printf("\t5 - Carregar Matriz de 2000 elementos\n");
        printf("\t6 - Carregar Matriz de 3000 elementos\n");
        printf("\t7 - Carregar Matriz de 4000 elementos\n\n");
        printf("\tEscolha uma opcao: ");
        scanf("%hd", &opcao);
        printf("\n");
    }

    *opEscolhida = opcao;
    caminho = getCaminhoArquivo(opcao);

    return(caminho);
}

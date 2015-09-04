#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "matEntrada.h"
#include "interfaces.h"
#define MAX_DOUBLE 2    // Auxiliar para gerar os valores do vetor chute x0
#define NUM_EXEC 10     // N�mero de vezes que o algoritmo ser� executado

int main (int argc, char* argv[]){

    // Vari�veis auxiliares para leitura do arquivo
    int i, j, k;        // Vars. controladoras de itera��es
    int auxInt;         // Armazena temporariamente a Ordem da Matriz
    short auxShort;     // Armazena temporariamente a fila da matriz a ser avaliada
    double auxDouble;   // Armazena temporariamente o erro m�x. permitido e vals. das matrizes
    long auxLong;       // Armazena temporariamente o n�m. max e num real de itera��es
    short opEscolhida;  // Op��o escolhida na tela inicial
    char* caminhoMat = telaInicial(-1,&opEscolhida); // Armazena o caminho para o arquivo da matriz
    double* resultado = NULL;   // Armazena o resultado do m�todo de jacobi richardson
    double* resIt = NULL;    // Armazena o resultado das NUM_EXEC itera��es do m�todo jacobi richardson
    double* x0 = NULL;       // Vetor x0 de chute inicial
    FILE* arqMat = NULL;     // Ponteiro para o arquivo de matrizes que ser� aberto

    clock_t begin, end;      // Clock de inicio e fim da execucao das iteracoes do jacobi richardson

    double* time_spent = (double*)calloc(NUM_EXEC,sizeof(double)); // Tempo gasto para rodar as iteracoes (em segundos)
    double mediaTimeSpent = 0;    // Media dos tempos gastos
    double desvioPadrao;

    // Cria a estrutura que armazenar� as matrizes
    MAT_ENTRADA* m = criarMatEntrada();

    arqMat = fopen(caminhoMat,"r");

    if(arqMat == NULL){
        printf("Erro ao abrir arquivo.");
    }
    else{
        // Le as primeiras informa��es do arquivo
        fscanf(arqMat, "%d", &auxInt);      // L� a ordem da matriz
        fscanf(arqMat, "%hd", &auxShort);   // L� a fila que ser� avaliada
        fscanf(arqMat, "%lf", &auxDouble);  // L� o erro permitido
        fscanf(arqMat, "%ld", &auxLong);    // L� o n�m. max de itera��es

        // Salva os valores na estrutura
        inicializarValsMatEntrada(m,auxInt,auxShort,auxDouble,auxLong,0);

        // Aloca espa�o para as matrizes A e B e diagonalAux, agora que sabemos a dimens�o
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
        // Fecha o arquivo agora que n�o precisamos mais dele
        fclose(arqMat);
    }

    resIt = (double*)calloc(getOrdem(m),sizeof(double));
    auxLong = 0;

    // Prepara as matrizes A e B dividindo os valores pela diagonal da A
    prepararMatrizes(m);

    // 10 itera��es do m�todo de jacobiRichardson
    for(i=0;i<NUM_EXEC;i++){
        //printf("Iniciando execucao %d do Jacobi Richardson.\n", i+1);
        // Aqui geramos o vetor de chute inicial x0
        // Cada pos x0[i] recebe um valores double entre 0 e MAX_DOUBLE
        x0 = (double*)malloc(getOrdem(m)*sizeof(double));
        srand((unsigned int)time(NULL));
        for(k=0;k<getOrdem(m);k++){
            x0[k] = ((float)rand()/(float)(RAND_MAX))*MAX_DOUBLE;
        }

        begin = clock();    // Come�a a contar o tempo de execucao
        resultado = jacobiRichardson(m,x0); // Execucao do metodo
        end = clock();      // Termina de contar contar o tempo de execucao
        time_spent[i] = (double)(end - begin) / CLOCKS_PER_SEC;    // Calcula o valor em segundos
        mediaTimeSpent += time_spent[i];

        salvarSaidasItermediarias(m,resultado,opEscolhida);

        auxLong += getRealIt(m);    // Salva o n� de iteracoes da execu��o "i"
        for(j=0;j<getOrdem(m);j++){
            resIt[j] += resultado[j]; // Salva o vetor resultado da execu��o "i"
        }

        free(x0);
        x0 = NULL;
        free(resultado);
        resultado = NULL;

        //printf("Fim da execucao %d.\n\n", i+1);
    }

    // Pega a media das 10 itera��es
    mediaTimeSpent = mediaTimeSpent/NUM_EXEC;
    desvioPadrao = calcDesvioPadrao(time_spent,mediaTimeSpent,NUM_EXEC);
    auxLong = auxLong/NUM_EXEC;
    for(i=0;i<getOrdem(m);i++){
        resIt[i] = resIt[i]/NUM_EXEC;
    }

    // Salva e imprime na tela
    salvarSaidaFinal(m,resIt,auxLong,opEscolhida);
    imprimirResultado(m,resIt,auxLong);

    for(i=0;i<NUM_EXEC;i++){
        printf("Tempo de execucao[%d]: %lfs.\n", i, time_spent[i]);
    }
    printf("Media do tempo de execucao: %lfs.\n", mediaTimeSpent);
    printf("Desvio padrao: %lfs.\n",desvioPadrao);

    // Libera a mem�ria alocada pela estrutura
    destruirMatEntrada(m);
    free(time_spent);

    return 0;
}


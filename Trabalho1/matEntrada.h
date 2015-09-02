#ifndef MATENTRADA_H_INCLUDED
#define MATENTRADA_H_INCLUDED

struct mat_entrada{
    int ordemMat;       // Ordem da Matriz de entrada
    short filaAval;     // Fila a ser avaliada
    double erroPerm;    // Erro máximo permitido
    long maxIt;         // Número máximo de iterações
    long realIt;        // Número de iterações realizadas
    double** matA;         // Matriz A
    double* matB;          // Matriz B
    double* diagonalAux;     // Salva o valor da diagonal da Matriz A
};

typedef struct mat_entrada MAT_ENTRADA;

// -- Funções de Manipulação da Estrutura

// Funções de alocação de memória
MAT_ENTRADA* criarMatEntrada(); // Aloca o espaço de memória para a estrutura
void alocarMatA(MAT_ENTRADA* m);    // Aloca o espaço de memória para a matriz A
void alocarMatB(MAT_ENTRADA* m);    // Aloca o espaço de memória para a matriz B
void alocarDiagonalAux(MAT_ENTRADA* m); // Aloca o espaço de memória para a diagonalAux

// Funções de inserção de valores (lidos dos arquivos)
void inserirOrdem(MAT_ENTRADA* m, int ordem);       // Insere a ordem das matrizes
void inserirFilaAval(MAT_ENTRADA* m, short fila);   // Insere qual fila será avaliada
void inserirErroPerm(MAT_ENTRADA* m, double erro);  // Insere o erro permitido
void inserirMaxIt(MAT_ENTRADA* m, long maxIt);      // Insere o num. maximo de iterações
void inserirRealIt(MAT_ENTRADA* m, long realIt);    // Insere o num. real de iterações realizadas
void inicializarValsMatEntrada(MAT_ENTRADA* m, int ordem, short fila, double erro, long maxIt, long realIt);    // Inicializa os valores da estrutura
void inserirElemMatA(MAT_ENTRADA* m, double elem, int posX, int posY); // Insere um elemento da matriz A
void inserirElemMatB(MAT_ENTRADA* m, double elem, int pos);    // Insere um elemento da matriz B
void inserirElemDiagonalAux(MAT_ENTRADA* m, double elem, int pos);    // Insere um elemento da diagonalAux


// Funções para recuperar valores
int getOrdem(MAT_ENTRADA* m);       // Retorna a ordem da matriz
short getFilaAval(MAT_ENTRADA* m);  // Retorna a fila que será avaliada
double getErroPerm(MAT_ENTRADA* m); // Retorna o erro permitido
long getMaxIt(MAT_ENTRADA* m);      // Retorna o máximo de iterações
long getRealIt(MAT_ENTRADA* m);     // Rretorna o número de iterações realizadas
double getElemMatA(MAT_ENTRADA* m, int posX, int posY);    // Retorna um elemento da matriz A
double getElemMatB(MAT_ENTRADA* m, int pos);   // Retorna um elemento da matriz B
double getElemDiagonalAux(MAT_ENTRADA* m, int pos); // Retorna um elemento da diagonalAux

// Funções de liberação de memória
void desalocarMatA(MAT_ENTRADA* m); // Libera o espaço de memória alocado pela matriz A
void desalocarMatB(MAT_ENTRADA* m); // Libera o espaço de memória alocado pela matriz B
void desalocarDiagonalAux(MAT_ENTRADA* m);  // Libera o espaço de memória alocado pela diagonalAux
void destruirMatEntrada(MAT_ENTRADA* m);   // Libera o espaço de memória alocado por toda a estrutura

// Funções de impressão
void imprimirMatA(MAT_ENTRADA* m);  // Imprime os elementos da Matriz A
void imprimirMatB(MAT_ENTRADA* m);  // Imprime os elementos da Matriz B
void imprimirInfosMatEntrada(MAT_ENTRADA* m);   // Imprime algumas infos. da Estrutura
void imprimirResultado(MAT_ENTRADA* m, double* res);

// Funções especificas do Método de Jacobi-Richardson
void prepararMatrizes(MAT_ENTRADA* m);
double* jacobiRichardson(MAT_ENTRADA* m, double *x0);
double calcErro(double *xk, double *xkMaisUm, int tam);
double calcNormaVetor(double *vet, int tam);
double* subtracaoVetores(double *vet1, double *vet2, int tam);


#endif // MATENTRADA_H_INCLUDED

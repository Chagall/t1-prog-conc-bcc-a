#ifndef MATENTRADA_H_INCLUDED
#define MATENTRADA_H_INCLUDED

struct mat_entrada{
    int ordemMat;       // Ordem da Matriz de entrada
    short filaAval;     // Fila a ser avaliada
    double erroPerm;    // Erro m�ximo permitido
    long maxIt;         // N�mero m�ximo de itera��es
    int** matA;         // Matriz A
    int* matB;          // Matriz B
};

typedef struct mat_entrada MAT_ENTRADA;

// -- Fun��es de Manipula��o da Estrutura

// Fun��es de aloca��o de mem�ria
MAT_ENTRADA* criarMatEntrada(); // Aloca o espa�o de mem�ria para a estrutura
void alocarMatA(MAT_ENTRADA* m);    // Aloca o espa�o de mem�ria para a matriz A
void alocarMatB(MAT_ENTRADA* m);    // Aloca o espa�o de mem�ria para a matriz B

// Fun��es de inser��o de valores (lidos dos arquivos)
void inserirOrdem(MAT_ENTRADA* m, int ordem);       // Insere a ordem das matrizes
void inserirFilaAval(MAT_ENTRADA* m, short fila);   // Insere qual fila ser� avaliada
void inserirErroPerm(MAT_ENTRADA* m, double erro);  // Insere o erro permitido
void inserirMaxIt(MAT_ENTRADA* m, long maxIt);      // Insere o num. maximo de itera��es
void inicializarValsMatEntrada(MAT_ENTRADA* m, int ordem, short fila, double erro, long maxIt);    // Inicializa os valores da estrutura
void inserirElemMatA(MAT_ENTRADA* m, int elem, int posX, int posY); // Insere um elemento da matriz A
void inserirElemMatB(MAT_ENTRADA* m, int elem, int pos);    // Insere um elemento da matriz B

// Fun��es para recuperar valores
int getOrdem(MAT_ENTRADA* m);       // Retorna a ordem da matriz
short getFilaAval(MAT_ENTRADA* m);  // Retorna a fila que ser� avaliada
double getErroPerm(MAT_ENTRADA* m); // Retorna o erro permitido
long getMaxIt(MAT_ENTRADA* m);      // Retorna o m�ximo de itera��es
int getElemMatA(MAT_ENTRADA* m, int posX, int posY);    // Retorna um elemento da matriz A
int getElemMatB(MAT_ENTRADA* m, int pos);   // Retorna um elemento da matriz B

// Fun��es de libera��o de mem�ria
void desalocarMatA(MAT_ENTRADA* m); // Libera o espa�o de mem�ria alocado pela matriz A
void desalocarMatB(MAT_ENTRADA* m); // Libera o espa�o de mem�ria alocado pela matriz B
void destruirMatEntrada(MAT_ENTRADA* m);   // Libera o espa�o de mem�ria alocado por toda a estrutura

// Fun��es de impress�o
void imprimirMatA(MAT_ENTRADA* m);  // Imprime os elementos da Matriz A
void imprimirMatB(MAT_ENTRADA* m);  // Imprime os elementos da Matriz B
void imprimirInfosMatEntrada(MAT_ENTRADA* m);   // Imprime algumas infos. da Estrutura

#endif // MATENTRADA_H_INCLUDED

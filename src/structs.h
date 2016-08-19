#ifndef STRUCTS_H
#define STRUCTS_H

#include <iostream>
#include <vector>
#include <time.h>
#include <string.h>

using namespace std;

//estrutura utilizada para armazenar, distancia, demanda e penalidades
struct celula{
	int id; //identifica o numero do n�
	double valor; //valor da celula
};

struct no{
	int id;//id do respectivo no
	int hub;//id do hub no qual est� alocado
};

//estrutura de uma instancia a ser carregada
struct instancia{
	string arquivo; //caminho do arquivo de entrada
	double FO_otimo; //valor da solu��o OTIMA para a determinada instancia
};

struct solucao{
	vector <int> hubs; //contem n�s que s�o hubs
	vector <int> hubs_bin;
	vector <no> alocacao; //contem lista ordenada de n�s n�o concentradores e suas aloca��es
};

//variaveis problema (Leitura de arquivo e parametros de regulagem)
struct DATA{

    //marcadores de tempo
    clock_t inicio;
    clock_t final;
    double tempo;

    double alvo;

	int nos; //quantidade de n�s da instancia trabalhada

	double alpha;
	double ex;

	double peso_cf; //peso da penalidade do custo fixo de instala�ao
	double peso_od; //peso da penalidade do O_i e D_i
	double peso_dist; //peso da penalidade da distancia

	int numFixoHubs;//N�mero fixo de hubs a ser instalado

	vector < vector<double> > cordenadas; //cordenadas dos n�s

	vector < vector<celula> > distancia; //matriz com distancia entre os n�s
	vector < vector<celula> > distanciaOrdenada; //vetores mais proximos do hub identificado pela linha

	vector < vector<double> > demanda; //matriz de demanda entre os n�s
	vector < vector<celula> > demandaOrigemOrdenada;  //matriz demanda com origem em i (W_ij)
	vector < vector<celula> > demandaDestinoOrdenada; //matriz demanda com com destino em i (W_ji)
	vector < celula > demandaTotalOrdenada;

	vector < int > pib;
	vector < int > populacao;
	vector < string > municipios;

	vector < int > ordemPib;

	vector <double> O;
	vector <double> D;
	vector <celula> O_ordenado;
	vector <celula> D_ordenado;

	vector < vector<celula> > demandaOD; //soma de demandas OD

	vector < double > custoIntalacao; //custo fixo de instala��o de cada hub
	vector < celula > custoInstalacaoOrdenado;

	vector < double > capacidade;
	vector < double > outrosValores;

	vector < celula > hubsPromissores; //armazena penalidades e os nos mais promissores a se tornar hub
    vector < celula > hubsPromissoresGrasp;

	//matriz binarias
	vector < vector<int> > ligacao; //z liga��es matriz inteira

	float T;
};

#endif

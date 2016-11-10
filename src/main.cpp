/*
* File:   main.cpp
* Author: Rodrigo Brito
*
* Aplica��o de Heur�stica Construtiva
* Hub inicializado com o mais proximo de todos atraves de ordenamento de matriz
*/
//#define CYCLE_HUB
//#define GRASP_ILS
//#define STOP_OPTIMAL

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <algorithm>//lower bound
#include <cfloat>
#include "structs.h"
#include "tsp.h"
#define INF DBL_MAX
using namespace std;

ofstream arquivoSaida;

//Trata operacoes relacionadas ao TSP
TSP * tsp;

//func�es e procedimentos
void resize(DATA * dados);
//double FO(solucao * s, DATA * dados);
void imprimeDATA(DATA * dados);
void imprimeVetor(vector < vector<double> >, DATA * dados);
void imprimeVetor(vector < vector<int> >, DATA * dados);
void imprimeSolucao(solucao * s);
void imprimeFluxoArestas(DATA * dados, solucao * s);
void imprimirCidades(DATA * dados, solucao * s);
void ordenaDistancia(DATA * dados);
void ordenaDemanda(DATA * dados);
void ordenaCustoInstalacao(DATA * dados);
void penalizaNos(DATA * dados);
void ordenaPenalidade(DATA * dados);
void ordenaPromissoresGRASP(DATA * dados);
int sorteioPosicao(int inicio, int fim);

//pre processamento
void inicializaSolucao_(DATA * dados, solucao * s);
void zerarSolucao(solucao * s);
void calcula_fluxo(DATA * dados);

//func�es ordenamento distancia
int partition(vector<celula>&, int, int);
int quickSort(vector<celula>&, int, int);
double Calcula_FO(DATA * dados, solucao *s);
double Calcula_FO2(DATA * dados, solucao *s);

//fun�oes de tempo
void iniciarCronometro(DATA * dados);
void finalizarCronometro(DATA * dados);

//fun�oes para arquivo de entrada e saida
DATA * reduzirDados(DATA * dados, int numeroCidades);
void lerArquivo(DATA * dados, char * arquivo);
void lerArquivoDistancia(DATA * dados, char * arquivo);
void lerArquivoPibPopulacao(DATA * dados, char * arquivo);
void lerArquivoOrdemPib(DATA * dados, char * arquivo);
void lerArquivoMunicipios(DATA * dados, char * arquivo);
void gerarDemanda(DATA * dados);
void gerarCustoFixo(DATA * dados);
void somaDemanda(DATA * dados);
void salvarResultado(DATA * dados, double FO, solucao * s);
void salvarMensagem(int i);
void lerInstacias(vector < instancia > * instancias, char * arquivo);

void instanciaIndividual( char * instancia, double alpha, double ex, int maxInter );
void processarInstancias(vector < instancia > * instancias);
void processarGoogleMapsMinas(int numCidades, float a, float b);

//manipula��o de solu��es
int buscaBin(int x, int e, int d, vector<int> v);//busca binaria de nos
bool isHub(int no, solucao * s); //verifica se dado no � hub
int hubMaisProximo(DATA * dados, solucao * s, int no);
void alocarNos(DATA * dados, solucao * s);

//controle de hubs
void addHub(int no, solucao * s);
void removeHub(int no, solucao * s);
void trocarAlocacao(int hub, int no, solucao * s);
int posInsercaoHub(int x, int e, int d, vector<int> v);

//Buscas locais
void buscaLocal_Shift_(DATA * dados, solucao * s, solucao * s_star, double * FO_star);
void buscaLocal_Shift(DATA * dados, solucao * s, solucao * s_star, double * FO_star);
void buscaLocal_Shift_P(DATA * dados, solucao * s, solucao * s_star, double * FO_star);

void buscaLocal_DeslocamentoAlocacao(DATA * dados, solucao * s, solucao * s_star, double * FO_star);
void buscaLocal_TrocaFuncao(DATA * dados, solucao * s, solucao * s_star, double * FO_star);
void buscaLocal_TrocaFuncao_(DATA * dados, solucao * s, solucao * s_star, double * FO_star);
void buscaLocal_HubPromissor_(DATA * dados, solucao * s, solucao * s_star, double * FO_star);
void buscaLocal_AdicionaHub_(DATA * dados, solucao * s, solucao * s_star, double * FO_star);
void buscaLocal_RemoveHub(DATA * dados, solucao * s, solucao * s_star, double * FO_star);
void buscaLocal_RemoveHub_(DATA * dados, solucao * s, solucao * s_star, double * FO_star);

void cicloBuscas(DATA * dados, solucao * s);

void VND(DATA * dados, solucao * s, solucao * s_star, double * FO_star);
void VNS(DATA * dados, solucao * s, solucao * s_star, double * FO_star);
void VNS_GRASP(DATA * dados, solucao * s, solucao * s_star, double * FO_star);
void VNS_TROCA(DATA * dados, solucao * s, solucao * s_star, double * FO_star);
void ILS(DATA * dados, solucao * s1, solucao * s_star, double * FO_star, int maxInter);

//pertuba��es aleat�rias
void pertubar_JucaoSolucoes(DATA * dados, solucao * s);
void pertubar_AdicionaHub(DATA * dados, solucao * s);

int main(int argc, char* argv[]){
    char * instancia;
    double alpha;
    double ex = -1;
    int maxIter = 0;
    if(argc == 5){
        instancia = argv[1];
        alpha = atof(argv[2]);
        ex = atof(argv[3]);
        maxIter = atoi(argv[4]);
    } else if (argc == 4){
        instancia = argv[1];
        alpha = atof(argv[2]);
        ex = atof(argv[3]);
    } else if (argc == 3){
        instancia = argv[1];
        alpha = atof(argv[2]);
    } else {
        cout<<"Parâmetros incorretos!"<<endl;
        cout<<"Uso: ./executavel instancia alfa ex maxIter (ex e maxIter são opcionais)"<<endl;
        cout<<"maxIter = Maximo de iteracoes do ILS, caso seja utilizado"<<endl;
        cout<<"Ex: ./main instancia.txt 0.2 0.7 50"<<endl;
        exit(EXIT_FAILURE);
    }

    #ifdef CYCLE_HUB
        cout<<"Cycle Hub Location Problem - Single Allocation"<<endl;
    #else
        cout<<"Hub Location Problem - Single Allocation"<<endl;
    #endif // CYCLE_HUB


    srand((unsigned)time(0));
    //processa isntancias
    instanciaIndividual(instancia, alpha, ex, maxIter);
    return 0;
}

void zerarSolucao(solucao * s, int nos){
    s->alocacao.resize(0);
    s->alocacao.resize(nos);
    s->hubs.resize(0);
    s->hubs_bin.resize(0);
    s->hubs_bin.resize(nos,0);
}

void instanciaIndividual(char * instancia, double alpha, double ex, int maxInter){
    DATA * dados; //estrutura temporia para armazenar os dados de cada instancia
    solucao * s1 = new solucao; //solu��o tempor�ria
    double * FO_star = new double; //melhor FO corrente da instancia
    solucao * s_star = new solucao; //melhor solu��o corrente
    dados = new DATA; //aloca espa�o para pacote com informa��o da instancia

    cout<<instancia<<endl;
    dados->alvo = 0;
    dados->ex = ex;

    //dados->numFixoHubs = 3;

    //leitura de arquivo da instancia i
    lerArquivo(dados, instancia);//pode passar ex
    tsp = new TSP( &dados->distancia );

    //definindo pesos para pontua��o no ordenamento dos hubs promissores
    dados->peso_cf = dados->nos;
    dados->peso_dist = 1; //estava 1.5
    dados->peso_od = alpha;//estava 0
    dados->alpha = alpha;

    calcula_fluxo(dados); //separa O_i D_i

    ordenaDistancia(dados);
    ordenaDemanda(dados);
    ordenaCustoInstalacao(dados);

    penalizaNos(dados);//faz a contagem do ordenamento para pontuar penalidades

    ordenaPenalidade(dados);//ordena para ver os mais promissores a se tornar hub

    zerarSolucao(s1 , dados->nos);
    zerarSolucao(s_star, dados->nos);
    *FO_star = INF;

    //inicia cronometro para contagem de tempo
    iniciarCronometro(dados);

    //GRASP Contrutivo
    inicializaSolucao_(dados, s1);

    *FO_star = Calcula_FO(dados, s1);
    *s_star = *s1;

    #ifdef GRASP_ILS
        ILS(dados, s1, s_star, FO_star, maxInter);
    #else
        //VNS com pertuba��o de jun��o de solu��es
        VNS_TROCA(dados, s1, s_star, FO_star);
    #endif // GRASP_ILS

    //finaliza contagem de tempo
    finalizarCronometro(dados);


    double FO = Calcula_FO(dados, s1);

    imprimeFluxoArestas(dados, s1);

    //salvando resultado em arquivo de saida
    salvarResultado(dados,FO,s1);
    delete tsp;
    delete FO_star;
    delete s_star;
    delete dados;
    delete s1;

}


void processarInstancias(vector < instancia > * instancias){
    DATA * dados; //estrutura temporia para armazenar os dados de cada instancia

    char * arquivo; //arquivo de entrada de cada instancia selecionada

    solucao * s1 = new solucao; //solu��o tempor�ria

    double * FO_star = new double; //melhor FO corrente da instancia

    solucao * s_star = new solucao; //melhor solu��o corrente


    for(int i = 0; i< instancias->size(); i++){ //percorre toda as listas de instancias
		cout << "--Intancia "<<i<<"--" << endl;
        dados = new DATA; //aloca espa�o para pacote com informa��o da instancia
        arquivo = &instancias->at(i).arquivo[0];//recebe endere�o da instancia i
        cout<<arquivo<<endl;
        dados->alvo = instancias->at(i).FO_otimo;//recebe FO �tima da instancia i

        //dados->numFixoHubs = 3;

        //leitura de arquivo da instancia i
        lerArquivo(dados, arquivo);
        tsp = new TSP( &dados->distancia );

        //definindo pesos para pontua��o no ordenamento dos hubs promissores
        dados->peso_cf = dados->nos;
        dados->peso_dist = 1; //estava 1.5
        dados->peso_od = dados->alpha;//estava 0

        calcula_fluxo(dados); //separa O_i D_i

        ordenaDistancia(dados);
        ordenaDemanda(dados);
        ordenaCustoInstalacao(dados);

        penalizaNos(dados);//faz a contagem do ordenamento para pontuar penalidades

        ordenaPenalidade(dados);//ordena para ver os mais promissores a se tornar hub

        zerarSolucao(s1 , dados->nos);
        zerarSolucao(s_star, dados->nos);
        *FO_star = INF;

        //inicia cronometro para contagem de tempo
        iniciarCronometro(dados);

        //GRASP Contrutivo
        #ifdef DEBUG
        cout<<"1 - Contrucao GRASP"<<endl;
        #endif // DEBUG
        inicializaSolucao_(dados, s1);
        //pertubar_JucaoSolucoes(dados, s1); //pertuba
        //inicializaSolucao(dados, s1);

        *FO_star = Calcula_FO(dados, s1);
        *s_star = *s1;
    	//Busca Shift
//    	buscaLocal_Shift_(dados, s1, s_star, FO_star);

		//Busca Troca Fun��o
//		buscaLocal_TrocaFuncao_(dados, s1, s_star, FO_star);

		//Busca Insere Hub
//		buscaLocal_AdicionaHub_(dados, s1, s_star, FO_star);
		//Remove Hub
//		buscaLocal_RemoveHub_(dados, s1, s_star, FO_star);

        //Ciclo de buscas - Shift, Troca de Fun��o, inser��o de Hub e remo��o de Hub
        //VND(dados, s1, s_star, FO_star);

        //Estrutura do VNS
//        VNS(dados, s1, s_star, FO_star);

        //VNS com pertuba��o de jun��o de solu��es
        #ifdef DEBUG
        cout<<"2 - VNS - TROCA"<<endl;
        #endif // DEBUG
        //VND(dados, s1, s_star, FO_star);

        //finaliza contagem de tempo
        finalizarCronometro(dados);


        double FO = Calcula_FO(dados, s1);
//        printf("Custo (FO) = %18.4f \n", FO);
//        imprimeSolucao(s1);

//        cout<<"-----"<<endl;
//        printf("Tempo: %4.4f ",dados->tempo);
//        cout<<"-----"<<endl;

        //salvando resultado em arquivo de saida
        #ifdef DEBUG
        cout<<"3 - Salvando resultado\n"<<endl;
        #endif // DEBUG
        salvarResultado(dados,FO,s1);
        delete tsp;
    }
    delete FO_star;
    delete s_star;
    delete dados;
    delete s1;
}


void processarGoogleMapsMinas(int numCidades, float a, float b){
    DATA * dados = new DATA; //aloca espa�o para pacote com informa��o da instancia
    char * arquivo; //arquivo de entrada de cada instancia selecionada
    solucao * s1 = new solucao; //solu��o tempor�ria
    double * FO_star = new double; //melhor FO corrente da instancia
    solucao * s_star = new solucao; //melhor solu��o corrente

    dados->nos = 853;
    resize(dados);

    //leitura de arquivo da instancia i
    char * arqDistancia = "instancia/DISTANCIA.txt"; //Matriz de dist�ncia das cidades de Minas Gerais
    lerArquivoDistancia(dados, arqDistancia);
    //cout<<"Distancia: "<<dados->distancia[852][851].valor<<endl;;

    //Leitura de PIB e Popula��o
    char * arqPibPopulacao = "instancia/PIB_POPULACAO.txt";
    lerArquivoPibPopulacao(dados, arqPibPopulacao);
    //cout<<"Ultimo pib e pop: "<<dados->pib[852]<<" / "<<dados->populacao[852]<<endl;

    //Leitura de arquivo com o nome das cidades
    char * arqMunicipios = "instancia/MUNICIPIOS.txt";
    lerArquivoMunicipios(dados, arqMunicipios);

    //Leitura de arquivo que contem a ordem de municipios por PIB
    char * arqOrdemPib = "instancia/ORDEM_PIB.txt";
    lerArquivoOrdemPib(dados, arqOrdemPib);

    //Gerar demanda pelo calculo do Ricardo/Bruno
    gerarDemanda(dados);

    //Armazena em vetor separado a somada das demandas de origem e destino para ordena��o
    somaDemanda(dados);

    //Gerar demanda pelo calculo do Ricardo/Bruno
    gerarCustoFixo(dados);

    //definindo pesos para pontua��o no ordenamento dos hubs promissores
    //dados->peso_cf = dados->nos;
    //dados->peso_dist = 1; //estava 1.5
    //dados->peso_od = dados->alpha;//estava 0

    calcula_fluxo(dados); //separa O_i D_i

    //ordenaDistancia(dados);
    //ordenaDemanda(dados);
    //ordenaCustoInstalacao(dados);

    //penalizaNos(dados);//faz a contagem do ordenamento para pontuar penalidades

    //ordenaPenalidade(dados);//ordena para ver os mais promissores a se tornar hub

    zerarSolucao(s1 , dados->nos);
    zerarSolucao(s_star, dados->nos);
    *FO_star = INF;

    //Reduz a instancia geral em cem dados
    DATA * dados_temp = reduzirDados(dados, numCidades);
    dados_temp->alvo = -DBL_MAX;
    //inicia cronometro para contagem de tempo
    iniciarCronometro(dados_temp);
    //GRASP Contrutivo
    dados_temp->alpha = a;
    dados_temp->T = b;

	for (int i = 0; i < dados_temp->nos; i++) {
		dados_temp->custoIntalacao.at(i) = 99999999 * float(dados_temp->pib[i]/dados_temp->populacao[i]);
		long double daux = 0.0;
		for (int j = 0; j < dados_temp->nos; j++) {
			if(i!=j)
			daux += dados_temp->demanda[i][j];
			else
			dados_temp->demanda[i][j] = 0;
		}
		dados_temp->O[i] = daux;
		daux = 0.0;
		for (int j = 0; j < dados_temp->nos; j++) {
			if(i!=j)
			daux += dados_temp->demanda[j][i];
		}
		dados_temp->D[i] = daux;
	}

    #ifdef DEBUG
    cout<<"1 - Contrucao GRASP"<<endl;
    #endif // DEBUG
    inicializaSolucao_(dados_temp, s1);


    //pertubar_JucaoSolucoes(dados, s1); //pertuba
    //inicializaSolucao(dados, s1);

    *FO_star = Calcula_FO(dados_temp, s1);
    *s_star = *s1;
    //Busca Shift
//    	buscaLocal_Shift_(dados, s1, s_star, FO_star);

    //Busca Troca Fun��o
//		buscaLocal_TrocaFuncao_(dados, s1, s_star, FO_star);

    //Busca Insere Hub
//		buscaLocal_AdicionaHub_(dados, s1, s_star, FO_star);
    //Remove Hub
//		buscaLocal_RemoveHub_(dados, s1, s_star, FO_star);

    //Ciclo de buscas - Shift, Troca de Fun��o, inser��o de Hub e remo��o de Hub
    //VND(dados, s1, s_star, FO_star);

    //Estrutura do VNS
//        VNS(dados, s1, s_star, FO_star);

    //VNS com pertuba��o de jun��o de solu��es
    #ifdef DEBUG
    cout<<"2 - VNS - TROCA"<<endl;
    #endif // DEBUG
    VNS_TROCA(dados_temp, s1, s_star, FO_star);

    //finaliza contagem de tempo
    finalizarCronometro(dados_temp);

    double FO = Calcula_FO(dados_temp, s1);
    #ifdef DEBUG
    printf("Custo (FO) = %18.4f \n", FO);
    //imprimeSolucao(s1);

    cout<<"-----"<<endl;
    printf("Tempo: %4.4f ",dados_temp->tempo);
    cout<<"-----"<<endl;

    //salvando resultado em arquivo de saida
    cout<<"3 - Salvando resultado\n"<<endl;
    #endif // DEBUG
    salvarResultado(dados_temp,FO,s1);

    delete arquivo;
    delete dados_temp;
    delete FO_star;
    delete s_star;
    delete dados;
    delete s1;

}

void salvarPlot(DATA * dados, solucao * s){

//    FILE *arquivoSaida;
//    arquivoSaida = fopen("PLOT.txt","a");
//    fprintf(arquivoSaida, "//Cidades = %d\n", dados->nos);
//    string hubs = "", ligacoes = "";
//    for(int i = 0; i < s->hubs.size(); i++){
//        hubs +=  "\t'"+dados->municipios[s->hubs[i]]+"'";
//        if(i < s->hubs.size()-1){
//            hubs += ",\n";
//        }
//    }
//    for(int i = 0; i < s->alocacao.size(); i++){
//        if(!isHub(s->alocacao[i].id, s)){
//            ligacoes += "\t{'no':'"+dados->municipios[s->alocacao[i].id]+"', 'concentrador':'"+dados->municipios[s->alocacao[i].hub]+"'}";
//        }
//        if(i < s->alocacao.size()-1 && !isHub(s->alocacao[i].id, s)){
//            ligacoes += ",\n";
//        }
//    }
//    fprintf(arquivoSaida,"solucao = { 'hubs':[\n%s\n],\n'ligacoes':[\n%s\n] };", hubs.c_str(), ligacoes.c_str());
//    fprintf(arquivoSaida, "\n\n");
//    fclose(arquivoSaida);
	    FILE *arquivoSaida;
	    arquivoSaida = fopen("alocacao.txt","a");
	    fprintf(arquivoSaida, "//Cidades = %d Alpha = %f  T=%f", dados->nos, dados->alpha, dados->T);
	    for(int i = 0; i < int(s->alocacao.size()); i++){
	        fprintf(arquivoSaida, "\n%d ", i);
	        fprintf(arquivoSaida, "\t%d ",s->alocacao[i].hub);
	    }
	    fprintf(arquivoSaida, "\n\n");
	    fclose(arquivoSaida);
}

void salvarResultado(DATA * dados, double FO, solucao * s){

    FILE *arquivoSaida;
    arquivoSaida = fopen("output/results.txt","a");
    double gap = (FO - dados->alvo)/dados->alvo;
    #ifdef CYCLE_HUB
    fprintf(arquivoSaida, "NOS=%d %.2f HUBS=%d TEMPO=%.4f FO=%18.4f\n", dados->nos, dados->alpha, s->hubs.size(), dados->tempo, FO, dados->alvo, gap);
    #else
    fprintf(arquivoSaida, "NOS=%d %.2f HUBS=%d TEMPO=%.4f FO=%18.4f OTIMO=%18.4f GAP=%3.4f\n", dados->nos, dados->alpha, s->hubs.size(), dados->tempo, FO, dados->alvo, gap);
    #endif // CYCLE_HUB
    //fprintf(arquivoSaida, " -- CIDADES --\n");
    /*for(int i = 0; i < s->hubs.size(); i++){
        fprintf(arquivoSaida, dados->municipios[s->hubs[i]].c_str());
        fprintf(arquivoSaida, "\t[%d]\n",s->hubs[i]);
    }*/
    //fprintf(arquivoSaida, " ---- \n\n");
    fclose(arquivoSaida);
    //salvarPlot(dados, s);
    #ifdef DEBUG
        imprimeSolucao(s);
        printf("TEMPO = %.4f\n", dados->tempo);
        printf("FO = %18.4f\n", Calcula_FO(dados, s));
    #endif // DEBUG
}

void salvarMensagem(int i){
    FILE *arquivoSaida;
    arquivoSaida = fopen("OUTPUT.txt","a");
    fprintf(arquivoSaida, "-- VALENDO -- %d \n", i);
    fclose(arquivoSaida);
}
void salvarResultado_antigo(DATA * dados, double FO){
    double gap = ((FO-dados->alvo)/dados->alvo*100);

    if(gap<0.0001)
        gap=0;

    arquivoSaida.open ("OUTPUT.txt",ios::app);
    arquivoSaida << "Inst_"<<dados->nos<<"_"<<dados->alpha<< setprecision(14)
    <<"\tFO: "<<FO
    <<"\tOtimo: "<<dados->alvo<<setprecision(4)
    <<"\tGAP: "<<gap<<"%"
    <<"\tTempo: "<<dados->tempo<<"\n";
    arquivoSaida.close();

}
void buscaLocal_DeslocamentoAlocacao(DATA * dados, solucao * s, solucao * s_star, double * FO_star){
	solucao * s_nova = new solucao;
	*s_nova = *s;
	double FO = *FO_star;
	//percorrer todos todos os n�s
	for(int i = 0; i < dados->nos; i++){
        if(!isHub(i, s_nova)){//se n�o for hub
            *s_nova = *s;
            for(int j = 0; j<s_nova->hubs.size(); j++){//troca aloca��o em todos os hubs
                s_nova->alocacao.at(i).hub = s_nova->hubs.at(j); //aloca no atual no hub;
                FO = Calcula_FO(dados, s_nova);
                if(FO < *FO_star){
                    *s = *s_nova;
                    *s_star = *s_nova;
                    *FO_star = FO;
                }
            }
        }
	}
}

void buscaLocal_Shift_(DATA * dados, solucao * s, solucao * s_star, double * FO_star){
//	cout<<"Busca local deslocamente de aloca��o"<<endl;
	solucao * s_nova = new solucao;
	*s_nova = *s;
	double FO_antes, FO_depois;
	int hub_solucao;
	int melhorAloc;
	//percorrer todos todos os n�s
	for(int i = 0; i < dados->nos; i++){
        if(!isHub(i, s_nova)){//se n�o for hub
            *s_nova = *s;
            FO_antes = Calcula_FO(dados,s_nova);
            hub_solucao = s_nova->alocacao.at(i).hub; //pega o hub que i est� alocado
            melhorAloc = hub_solucao; //melhor aloca��o � a atual aloca��o, at� o momento
            for(int j = 0; j<s_nova->hubs.size(); j++){//percorre todos os hubs
                if(s_nova->hubs.at(j) != hub_solucao){ //verifica se j� n�o � o hub atual
                    s_nova->alocacao.at(i).hub = s_nova->hubs.at(j); //aloca no atual no hub;
                    FO_depois = Calcula_FO(dados, s_nova);
                    if(FO_depois < FO_antes){
                        FO_antes = FO_depois;
                        melhorAloc = s_nova->hubs.at(j);
                        if(FO_depois < *FO_star){
                            *FO_star = FO_depois;
                            *s_star = *s_nova;
                        }
                    }
                    *s_nova = *s;
                }
            }
            s->alocacao.at(i).hub = melhorAloc;
        }
	}
}

void buscaLocal_AdicionaHub_(DATA * dados, solucao * s, solucao * s_star, double * FO_star){
//    cout<<"Busca Local Adiciona Hub"<<endl;
//    cout<<"-----"<<endl;
    solucao * s_nova = new solucao; //variavel temporaria para testes
    *s_nova = *s; //faz copia da variavel original
    double FO = Calcula_FO(dados, s_nova);
    double novaFO = FO;
//    printf("FO inicial = %18.4f \n", FO);

    for (int n = 0; n < dados->nos; n++) { //percorre todos os nos
		if (!isHub(n, s_nova)) { //se o no n�o for hub
            addHub(n, s_nova); //adiciona o no N a solucao temporaria
            alocarNos(dados, s_nova); //aloca no ao hub mais proximo
            novaFO = Calcula_FO(dados, s_nova); //calcula FO total
            if(novaFO < FO){ //se houver melhora, a solu��o temporaria passa ser a atual e sai do loo
                *s = *s_nova;
                FO = novaFO;
                if(novaFO < *FO_star){
                    *s_star = *s_nova;
                    *FO_star = novaFO;
                }
            }else{//caso nao houver melhora, reseta a solu��o temporaria
                *s_nova = *s;
            }
		}
    }
//    printf("FO final = %18.4f \n", FO);
    delete s_nova;
}

void buscaLocal_AdicionaHub(DATA * dados, solucao * s, solucao * s_star, double * FO_star){
    solucao * s_nova = new solucao; //variavel temporaria para testes
    *s_nova = *s; //faz copia da variavel original
    double FO = Calcula_FO(dados, s_nova);//Calcula FO da atual solucao
    double melhorFO = FO; //inicializa Melhor FO com FO atual

    bool melhora = false; //verifica se houve melhora no processa

    for (int n = 0; n < dados->nos; n++) { //percorre todos os nos
		if (!isHub(n, s_nova)) { //se o no n�o for hub
            addHub(n, s_nova); //adiciona o no N a solucao temporaria
            alocarNos(dados, s_nova); //aloca no ao hub mais proximo
            FO = Calcula_FO(dados, s_nova); //calcula FO total
            if(FO < melhorFO){ //se houver melhora, a solu��o temporaria passa ser a atual e sai do loop
                melhora = true;
                melhorFO = FO;
                *s = *s_nova;
                break;
            }else{//caso nao houver melhora, reseta a solu��o temporaria
                *s_nova = *s;
            }
		}
    }
    if(melhora){//se houve rmelhora, atualiza star
        *s_star = *s_nova;
        *FO_star = melhorFO;
    }
    delete s_nova;
}


void buscaLocal_RemoveHub_(DATA * dados, solucao * s, solucao * s_star, double * FO_star){
//    cout<<"Busca Local Remover Hub"<<endl;
//    cout<<"-----"<<endl;
    solucao * s_nova = new solucao;
    *s_nova = *s;
    double FO = Calcula_FO(dados, s_nova);
//    printf("FO inicial = %18.4f \n", FO);

    double novaFO = FO;

    for (int n = 0; n < dados->nos; n++) { //percorre todos os nos
        #ifdef CYCLE_HUB
		if (isHub(n, s_nova) && s_nova->hubs.size() > 3) { //se o n� for hub e tiver mais de um hub
		#else
		if (isHub(n, s_nova) && s_nova->hubs.size() > 1) { //se o n� for hub e tiver mais de um hub
		#endif // CYCLE_HUB
            removeHub(n, s_nova);//remove da lista de hubs e desaloca
            alocarNos(dados, s_nova);//aloca nos ao mais proximo
            novaFO = Calcula_FO(dados, s_nova);//calcula FO
            if(novaFO < FO){//veirfica se houve melhora e salva solu��o caso sim
                *s = *s_nova;
                FO = novaFO;
                if(novaFO < *FO_star){
                    *s_star = *s_nova;
                    *FO_star = novaFO;
                }
            }else{//caso contrario, limpa solu��o
                *s_nova = *s;
            }
		}
    }
//    printf("FO final = %18.4f \n", FO);

    delete s_nova;
}
void buscaLocal_RemoveHub(DATA * dados, solucao * s, solucao * s_star, double * FO_star){
    solucao * s_nova = new solucao;
    *s_nova = *s;
    double FO = Calcula_FO(dados, s_nova);

    double melhorFO = FO;
    bool melhora = false;

    for (int n = 0; n < dados->nos; n++) { //percorre todos os nos
		if (isHub(n, s_nova)) { //se o n� for hub
            removeHub(n, s_nova);//remove da lista de hubs e desaloca
            alocarNos(dados, s_nova);//aloca nos ao mais proximo
            FO = Calcula_FO(dados, s_nova);//calcula FO
            if(FO < melhorFO){//veirfica se houve melhora e salva solu��o caso sim
                melhora = true;
                melhorFO = FO;
                *s = *s_nova;
                break;
            }else{//caso contrario, limpa solu��o
                *s_nova = *s;
            }
		}
    }
    if(melhora){
        *s_star = *s_nova;
        *FO_star = melhorFO;
    }

    delete s_nova;
}

void buscaLocal_TrocaFuncao_(DATA * dados, solucao * s, solucao * s_star, double * FO_star){
//	cout<<"Busca Local troca de funcao"<<endl;
//    cout<<"-----"<<endl;
    solucao * s_nova = new solucao; //solucao temporaria para testar combinacoes
    *s_nova = *s;
    double FO_antes, FO_depois; //computa fo inicial
//    printf("FO inicial = %18.4f \n", FO);
    for (int i = 0; i < s_nova->hubs.size(); i++) {//percorre todos os hubs

        *s_nova = *s; //reseta solu��o para novos testes
        FO_antes = Calcula_FO(dados, s_nova); //a cada teste de hub computa nova FO
        double melhorNo = s_nova->hubs[i]; //melhor no � o atual de in�cio

        for (int j = 0; j < dados->nos; j++) {// percorre os nos
            //if ( (s_nova->alocacao.at(n).hub == j) && (!isHub(n, s_nova)) ) {//se n � alocado a j e n n�o for hub
            if(s_nova->alocacao.at(j).hub == s_nova->hubs.at(i) && !isHub(j, s_nova)){//todos os nos que n�o forem concentradores e estiverem alocados ao no selecionado
                int ehHub = (s_nova->hubs_bin[j]);
                int noRemovido = s_nova->hubs.at(i);
                int jAlocado = s_nova->alocacao.at(j).hub;
                int noAdicionado = j;

                trocarAlocacao(s_nova->hubs.at(i),j,s_nova);//Faz troca e aloca��o

                FO_depois = Calcula_FO(dados, s_nova);

                if(FO_depois < FO_antes){//verifica se houve melhora para atualizar a solu��o

                    FO_antes = FO_depois; // passa a ser melgor FP
                    melhorNo = j;

                    if(FO_depois < *FO_star){
                        *FO_star = FO_depois;
                        *s_star = *s_nova;
                    }
                }
                *s_nova = *s; //reseta solu��o para novos testes
            }
        }
        if(melhorNo != s_nova->hubs[i]){//se encontrou melhor candidato faz troca definitiva
            trocarAlocacao(s_nova->hubs.at(i),melhorNo,s);//Faz troca e aloca��o definitiva
        }
    }

    delete s_nova;
}

void buscaLocal_HubPromissor_(DATA * dados, solucao * s, solucao * s_star, double * FO_star){
    solucao * s_nova = new solucao; //solucao temporaria para testar combinacoes
    *s_nova = *s;
    double FO = *FO_star; //computa fo inicial

    vector<int> hubsPromissores;
    for(int i = 0; i<dados->nos; i++){
        if(dados->hubsPromissoresGrasp[i].valor > 0){
            hubsPromissores.push_back(dados->hubsPromissoresGrasp[i].id);
        }else{
            break;
        }
    }

    for (int i = 0; i < s_nova->hubs.size(); i++) {//percorre todos os hubs
        // hub � p s_nova->hubs[i]
        for (int j = 0; j < hubsPromissores.size(); j++) {// percorre os hubs promissores
            //if ( (s_nova->alocacao.at(n).hub == j) && (!isHub(n, s_nova)) ) {//se n � alocado a j e n n�o for hub
            if(!isHub(hubsPromissores[j], s_nova)){//todos os nos que n�o forem concentradores e estiverem alocados ao no selecionado

                //s_nova->hubs[i] passa a ser n� e no j passa a ser nhub
                addHub(hubsPromissores[j],s_nova);//adiciona o hub promissor
                removeHub(s_nova->hubs.at(i), s_nova);//remove o hub atual testado para troca

                alocarNos(dados, s_nova);

                FO = Calcula_FO(dados, s_nova);

                if(FO < *FO_star){//verifica se houve melhora para atualizar a solu��o
                    *s = *s_nova;
                    *s_star = *s_nova;
                    *FO_star = FO;
                }else{
                    *s_nova = *s;
                }
            }
        }
    }
    delete s_nova;
}



void buscaLocal_TrocaFuncao(DATA * dados, solucao * s, solucao * s_star, double * FO_star){
    solucao * s_nova = new solucao; //solucao temporaria para testar combinacoes
    solucao * s_melhor = new solucao; //guarda a melhor solucao das combinacoes
    *s_nova = *s;
    double FO = Calcula_FO(dados, s_nova); //computa fo inicial
    double melhorFO = FO;

    bool melhora = false;

    for (int j = 0; j < dados->nos; j++) {//percorre todos os nos
		if (isHub(j, s_nova)) { //se ele for hub processa
			for (int n = 0; n < dados->nos; n++) {
                if ( (s_nova->alocacao.at(n).hub == j) && (!isHub(n, s_nova)) ) {//se n � alocado a j e n n�o for hub
                    //n passa a ser hub e j passa a ser n�
                    addHub(n,s_nova);
                    removeHub(j, s_nova);
                    alocarNos(dados, s_nova);

                    FO = Calcula_FO(dados, s_nova);

                    if(FO < melhorFO){//verifica se houve melhora para atualizar a solu��o
                        melhora = true;
                        *s = *s_nova;
                        melhorFO = FO;
                        break;
                    }else{//caso contrario limpa solu��o temporaria
                        *s_nova = *s;
                    }
				}
			}
		}
    }
    if(melhora){
        *s_star = *s_nova;
        *FO_star = melhorFO;
    }
    delete s_nova;
}
//Busca local com melhor aprimorante
void buscaLocal_Shift(DATA * dados, solucao * s, solucao * s_star, double * FO_star){
    solucao * s_nova = new solucao;
    *s_nova = *s;
    double FO = Calcula_FO(dados, s_nova);

    //melhores valores da rodada
    no melhorAlocacao;
    double melhorAlocacaoFO = FO;
    bool melhora=true;
    bool trocar=false; //quando houver uma melhora no processo geral
    //for(int h = 0; h<dados->nos; h++){
    while(melhora){
    	melhora = false;
	    for(int i = 0; i<dados->nos; i++){
	        int melhorHub = -1;
	        double melhorFO = FO;

	        //percorre lista de hubs para testes
	        double novaFO;
	        for(int j = 0; j<s_nova->hubs.size(); j++){

	            s_nova->alocacao[i].hub = s_nova->hubs[j];

                double novaFO = Calcula_FO(dados, s_nova);
	            if(novaFO < melhorFO){
	                melhorHub = s_nova->hubs[j];
	                melhorFO = novaFO;
	            }
	        }
	        //se houver uma solu��o melhor
	        if(melhorHub != -1 && melhorFO < melhorAlocacaoFO){
	        	//guarda as informa��es da melhor altera��o
	        	melhorAlocacao.id = i;
	        	melhorAlocacao.hub = melhorHub;
	            melhorAlocacaoFO = melhorFO;

	            //sinaliza que houve uma melhora
	            melhora = true;
	            trocar = true;
	            break;
	        }
	        *s_nova = *s;//zera para tentar em novo no
	    }

	    //usa a melhor re-aloca��o dos testes e repete at� nao melhorar mais
	    if(trocar){
	    	s_nova->alocacao[ melhorAlocacao.id ] = melhorAlocacao;
			FO = melhorAlocacaoFO;
			*s = *s_nova;
			break;
	    }
    }

    delete s_nova;
}
//Busca local com primeira de melhora
void buscaLocal_Shift_P(DATA * dados, solucao * s, solucao * s_star, double * FO_star){
    solucao * s_nova = new solucao;
    *s_nova = *s;
    double melhorFO = Calcula_FO(dados, s_nova);

    bool melhora = false; //registra se houve melhora na busca

    for(int i = 0; i<dados->nos; i++){//percorre todos os nos
        if(!isHub(i, s)){ //n�o pode ser hub
            for(int j = 0; j<s_nova->hubs.size(); j++){//percorre lista de hubs
                *s_nova = *s;
                s_nova->alocacao[i].hub = s_nova->hubs[j];//troca o hub do n� 'i'
                double novaFO = Calcula_FO(dados, s_nova);//calcula FO

                if(novaFO < melhorFO){//se melhorou ent�o troca
                    melhora = true;
                    *s = *s_nova;
                    melhorFO = novaFO;
                    break; //achou melhoria ent�o PARAR
                }
            }
        }
    }

    //se houver melhora ele altera a STAR
    if(melhora){
        *s_star = *s_nova;
        *FO_star = melhorFO;
    }
    //printf("FO final = %18.4f \n", FO);
    delete s_nova;
}

int sorteioPosicao(int inicio, int fim){
    int random_integer;
    int lowest=inicio, highest=fim-1;
    int range=(highest-lowest)+1;

    random_integer = lowest+int(range*(rand()/(RAND_MAX + 1.0)));

	return random_integer;
}

int hubMaisProximo(DATA* dados, solucao * s, int no){
	int prox = no;
	for (int i = 0; i < dados->nos; i++){
		if (isHub(dados->distanciaOrdenada[no][i].id, s)){
			prox = dados->distanciaOrdenada[no][i].id;
			break;
		}
	}
	return prox;
}

int hubMaisProximo2(DATA* dados, solucao * s, int no){
	int prox = s->hubs[0];
	for (int i = 1; i < s->hubs.size(); i++){
		if (dados->distancia[s->hubs[i]][no].valor < dados->distancia[no][prox].valor){
			prox = s->hubs[i];
		}
	}
	return prox;
}

//aloca todos os n�s ao hub mais proximo da lista de hubs da solu��o
void alocarNos(DATA* dados, solucao* s){
	s->alocacao.resize(dados->nos);
	for (int i = 0; i < dados->nos; i++){
		s->alocacao[i].id = i;
		s->alocacao[i].hub = hubMaisProximo2(dados, s, s->alocacao[i].id);
	}
}

//busca binaria de determinado valor em um vetor de inteiros
int buscaBin(int x, int e, int d, vector<int> v) {
	if (e > d) return -1;
	else {
		int m = (e + d) / 2;
		if (v[m] == x) return m;
		if (v[m] < x)
			return buscaBin(x, m + 1, d, v);
		else
			return buscaBin(x, e, m - 1, v);
	}
}

bool isHub(int no, solucao * s){
	if(s->hubs_bin[no] == 1){
		return true;
	}else{
		return false;
	}
}

bool isHub_ANTIGA(int no, solucao * s){
	bool encontrou = false;
	for(int i = 0; i< s->hubs.size(); i++){
		if(s->hubs[i] == no){
			encontrou = true;
			break;
		}
	}
	return encontrou;
}

void iniciarCronometro(DATA * dados){
    dados->inicio = clock();
}

void finalizarCronometro(DATA * dados){
    dados->final = clock();
    dados->tempo = double(dados->final - dados->inicio)/CLOCKS_PER_SEC;
}

bool isHubBin(int no, solucao * s){
	if (s->hubs.size() == 0){
		return false;
	}
	else if (buscaBin(no, 0, s->hubs.size() - 1, s->hubs) != -1){
		return true;
	}
	else{
		return false;
	}
}

void VND(DATA * dados, solucao * s, solucao * s_star, double * FO_star) {
	int vizinhanca = 1;
	solucao * s_nova = new solucao;
	double FO_depois;
	double FO_antes = Calcula_FO(dados, s);
	while (vizinhanca <= 4) { //explorando todas vizinhancas
		switch (vizinhanca) {
            case 1:
                *s_nova = *s; //reseta solucao
                FO_antes = Calcula_FO(dados, s_nova);
                //buscaLocal_Shift_P(dados, s_nova, s_star, FO_star);//aplica busca local
                buscaLocal_Shift_(dados, s_nova, s_star, FO_star);

                FO_depois = Calcula_FO(dados, s_nova);//calcula nova FO

                if (FO_depois < FO_antes) {//verifica se houve melhoria

                    //guarda nova FO
                    FO_antes = FO_depois;
                    //nova solucao passa ser a atual
                    *s = *s_nova;

                    //se chegar na FO otima para
                    #ifdef STOP_OPTIMAL
                    if (FO_depois <= dados->alvo + 0.1) {
//                        cout<<"ALCAN�O OTIMO"<<endl;
                        break;
                    }
                    #endif // STOP_OPTIMAL
                    vizinhanca = 1;

                }else{
                    vizinhanca++;
                }
            case 2: //Busca local - vizinhanca 2

                *s_nova = *s; //reseta solucao
                FO_antes = Calcula_FO(dados, s_nova);

                buscaLocal_TrocaFuncao_(dados, s_nova, s_star, FO_star);//aplica busca local

                FO_depois = Calcula_FO(dados, s_nova);//calcula nova FO

                if (FO_depois < FO_antes) {//verifica se houve melhoria

                    //guarda nova FO
                    FO_antes = FO_depois;
                    //nova solucao passa ser a atual
                    *s = *s_nova;

                    //se chegar na FO otima para
                    #ifdef STOP_OPTIMAL
                    if (FO_depois <= dados->alvo + 0.1) {
//                        cout<<"ALCAN�O OTIMO"<<endl;
                        break;
                    }
                    #endif // STOP_OPTIMAL
                    vizinhanca = 1;

                }else{
                    vizinhanca++;
                }
            case 3://Busca local - vizinhanca 3
                *s_nova = *s; //reseta solucao
                FO_antes = Calcula_FO(dados, s_nova);

                buscaLocal_AdicionaHub_(dados, s_nova, s_star, FO_star); //aplica busca local

                FO_depois = Calcula_FO(dados, s_nova);//calcula nova FO

                if (FO_depois < FO_antes) {//verifica se houve melhoria

                    //guarda nova FO
                    FO_antes = FO_depois;
                    //nova solucao passa ser a atual
                    *s = *s_nova;

                    //se chegar na FO otima para
                    #ifdef STOP_OPTIMAL
                    if (FO_depois <= dados->alvo + 0.1) {
//                        cout<<"ALCAN�O OTIMO"<<endl;
                        break;
                    }
                    #endif // STOP_OPTIMAL
                    vizinhanca = 1;

                }else{
                    vizinhanca++;
                }
            case 4://Busca local - vizinhanca 4
                #ifdef CYCLE_HUB
                if(s->hubs.size() > 3){
                #else
                if(s->hubs.size() > 1){
                #endif // CYCLE_HUB
                    *s_nova = *s; //reseta solucao
                    FO_antes = Calcula_FO(dados, s_nova);

                    buscaLocal_RemoveHub_(dados, s_nova, s_star, FO_star);//aplica busca local

                    FO_depois = Calcula_FO(dados, s_nova);//calcula nova FO

                    if (FO_depois < FO_antes) {//verifica se houve melhoria

                        //guarda nova FO
                        FO_antes = FO_depois;
                        //nova solucao passa ser a atual
                        *s = *s_nova;

                        //se chegar na FO otima para
                        #ifdef STOP_OPTIMAL
                        if (FO_depois <= dados->alvo + 0.1) {
//                            cout<<"ALCAN�O OTIMO"<<endl;
                            break;
                        }
                        #endif // STOP_OPTIMAL
                        vizinhanca = 1;

                    }else{
                        vizinhanca++;
                    }
                }else{
                     vizinhanca++;
                }
		}
		if (FO_depois <= dados->alvo + 0.1) {
			break;
		}
	}
	delete s_nova;
}


double Calcula_FO2(DATA * dados, solucao *s) {

	long double FO = 0;
	for (int i = 0; i <s->hubs.size(); i++){
		FO += dados->custoIntalacao[s->hubs[i]]; // total custo fixo de instala��o
	}

	for (int j = 0; j<s->alocacao.size(); j++){
		//custo entre n� e seu hub
		FO = FO + dados->distancia[s->alocacao[j].id][s->alocacao[j].hub].valor
			* (dados->O[s->alocacao[j].id] + dados->D[s->alocacao[j].id]);
				//TODO Perguntar o bruno sobre isso
				//double dif = 2*dados->distancia[s->alocacao[j].id][s->alocacao[j].hub].valor - dados->T;
				//if(dif>0){
					//FO *= 100;
				//}

        #ifdef CYCLE_HUB
        for (int w = j + 1; w<s->alocacao.size(); w++){
            FO = FO + dados->demanda[s->alocacao[j].id][s->alocacao[w].id] * dados->alpha *
                dados->distancia[s->alocacao[j].hub][s->alocacao[w].hub].valor +
                dados->demanda[s->alocacao[w].id][s->alocacao[j].id] * dados->alpha *
                tsp->getSubTourDistance(&s->hubs, s->alocacao[w].hub, s->alocacao[j].hub, false);
        }
		#else//custo entre os hubs (Hub Location Problem)
        for (int w = j + 1; w<s->alocacao.size(); w++){
            FO = FO + dados->demanda[s->alocacao[j].id][s->alocacao[w].id] * dados->alpha *
                dados->distancia[s->alocacao[j].hub][s->alocacao[w].hub].valor +
                dados->demanda[s->alocacao[w].id][s->alocacao[j].id] * dados->alpha *
                dados->distancia[s->alocacao[w].hub][s->alocacao[j].hub].valor; //Distância direta entre os dois HUBS
        }
		#endif // CYCLE_HUB
	}
	return FO;

}

double Calcula_FO(DATA * dados, solucao *s) {

	long double FO = 0;
	for (int i = 0; i <s->hubs.size(); i++){
		FO += dados->custoIntalacao[s->hubs[i]]; // total custo fixo de instala��o
	}

	for (int j = 0; j<s->alocacao.size(); j++){
		//custo entre n� e seu hub
		FO = FO + dados->distancia[s->alocacao[j].id][s->alocacao[j].hub].valor
			* (dados->O[s->alocacao[j].id] + dados->D[s->alocacao[j].id]);
				//TODO Perguntar o bruno sobre isso
				//double dif = 2*dados->distancia[s->alocacao[j].id][s->alocacao[j].hub].valor - dados->T;
				//if(dif>0){
					//FO *= 100;
				//}

        #ifdef CYCLE_HUB
        for (int w = j + 1; w<s->alocacao.size(); w++){
            FO = FO + dados->demanda[s->alocacao[j].id][s->alocacao[w].id] * dados->alpha *
                tsp->getSubTourDistance(&s->hubs, s->alocacao[j].hub, s->alocacao[w].hub, false) +
                dados->demanda[s->alocacao[w].id][s->alocacao[j].id] * dados->alpha *
                tsp->getSubTourDistance(&s->hubs, s->alocacao[w].hub, s->alocacao[j].hub, false);
        }
		#else//custo entre os hubs (Hub Location Problem)
        for (int w = j + 1; w<s->alocacao.size(); w++){
            FO = FO + dados->demanda[s->alocacao[j].id][s->alocacao[w].id] * dados->alpha *
                dados->distancia[s->alocacao[j].hub][s->alocacao[w].hub].valor +
                dados->demanda[s->alocacao[w].id][s->alocacao[j].id] * dados->alpha *
                dados->distancia[s->alocacao[w].hub][s->alocacao[j].hub].valor; //Distância direta entre os dois HUBS
        }
		#endif // CYCLE_HUB
	}
	return FO;

}

void imprimeFluxoArestas(DATA * dados, solucao * s){
    imprimeSolucao(s);
    FILE *arquivoSaida;
    printf("FO = %.4f\n", Calcula_FO(dados, s));
    arquivoSaida = fopen("output/fluxo.txt","a");
    fprintf(arquivoSaida, "----------------\n");
    fprintf(arquivoSaida, "N = %d ALFA = %.2f\n", dados->nos, dados->alpha);
    fprintf(arquivoSaida, "----------------\n");
    fprintf(arquivoSaida, "-- Alocacoes ---\n");
    #ifdef CYCLE_HUB
        vector<double> arestas( s->hubs.size() , 0);
        cout << "-- FLUXO: Arestas Alocacao --"<<endl;
        for (int j = 0; j<s->alocacao.size(); j++){
            //Fluxo entre no nao-concentrador e hub
            cout << s->alocacao[j].id << " " << s->alocacao[j].hub << " " << dados->O[s->alocacao[j].id] + dados->D[s->alocacao[j].id] << endl;
            fprintf(arquivoSaida, "%d %d %.4f\n", s->alocacao[j].id+1, s->alocacao[j].hub+1, dados->O[s->alocacao[j].id] + dados->D[s->alocacao[j].id]);

            //Quantidade de arestas e a mesma que a quantidade de concentradores
            for (int w = j + 1; w<s->alocacao.size(); w++){
                //Fluxo de ida
                vector<int> arestas_percorridas_origem = tsp->getEdges( &s->hubs, s->alocacao[j].hub, s->alocacao[w].hub );
                for(int k = 0; k < arestas_percorridas_origem.size(); k++){
                    arestas.at(arestas_percorridas_origem.at(k)) += dados->demanda[s->alocacao[j].id][s->alocacao[w].id];
                }
                //Fluxo de volta
                vector<int> arestas_percorridas_destino = tsp->getEdges( &s->hubs, s->alocacao[w].hub, s->alocacao[j].hub );
                for(int k = 0; k < arestas_percorridas_destino.size(); k++){
                    arestas.at(arestas_percorridas_destino.at(k)) += dados->demanda[s->alocacao[w].id][s->alocacao[j].id];
                }
            }
        }

        cout << "-- FLUXO: Anel concentradores" << endl;
        fprintf(arquivoSaida, "-- Concentradores ---\n");
        for( int m = 0; m < arestas.size(); m++ ){
            if(m == arestas.size() - 1){
                cout << s->hubs.at(m) <<" "<< s->hubs.at(0)<<" "<<arestas.at(m)<<endl;
                fprintf(arquivoSaida, "%d %d %.4f\n", s->hubs.at(m)+1, s->hubs.at(0)+1, arestas.at(m));
            }else{
                cout << s->hubs.at(m) <<" "<< s->hubs.at(m+1)<<" "<<arestas.at(m)<<endl;
                fprintf(arquivoSaida, "%d %d %.4f\n", s->hubs.at(m)+1, s->hubs.at(m+1)+1, arestas.at(m));
            }
        }
    #else
        vector< vector<double> > arestas( dados->nos ,  vector<double>(dados->nos, 0) );

        cout << "-- FLUXO: Arestas Alocacao --"<<endl;
        for (int j = 0; j<s->alocacao.size(); j++){
            //Fluxo entre no nao-concentrador e hub
            cout << s->alocacao[j].id << " " << s->alocacao[j].hub << " " << dados->O[s->alocacao[j].id] + dados->D[s->alocacao[j].id] << endl;
            fprintf(arquivoSaida, "%d %d %.4f\n", s->alocacao[j].id+1, s->alocacao[j].hub+1, dados->O[s->alocacao[j].id] + dados->D[s->alocacao[j].id]);

            //Quantidade de arestas e a mesma que a quantidade de concentradores
            for (int w = j + 1; w<s->alocacao.size(); w++){
                arestas[ s->alocacao[j].hub ][ s->alocacao[w].hub ] += dados->demanda[s->alocacao[j].id][s->alocacao[w].id];
                arestas[ s->alocacao[j].hub ][ s->alocacao[w].hub ] += dados->demanda[s->alocacao[w].id][s->alocacao[j].id];
            }
        }
        cout << "-- FLUXO: Anel concentradores" << endl;
        fprintf(arquivoSaida, "-- Concentradores ---\n");
        for( int m = 0; m < s->hubs.size(); m++ ){
            for( int n = m+1; n < s->hubs.size(); n++ ){
                cout << s->hubs.at(m)+1<<" "<< s->hubs.at(n)+1<<" "<<arestas[s->hubs.at(m)][s->hubs.at(n)]+arestas[s->hubs.at(n)][s->hubs.at(m)]<<endl;
                fprintf(arquivoSaida, "%d %d %.4f\n", s->hubs.at(m)+1, s->hubs.at(n)+1, arestas[s->hubs.at(m)][s->hubs.at(n)]+arestas[s->hubs.at(n)][s->hubs.at(m)]);
            }
        }
    #endif // CYCLE_HUB
    fclose(arquivoSaida);
}

void imprimeSolucao(solucao * s){
	cout << "Hubs = ";
	for (int i = 0; i<s->hubs.size(); i++){
		cout << s->hubs[i] << " ";
	}
	cout << "\nNos[hub aloc] = ";
	for (int i = 0; i<s->alocacao.size(); i++){
		cout << s->alocacao[i].id << "[" << s->alocacao[i].hub << "] ";
	}
	cout << endl;
}

void imprimirCidades(DATA * dados, solucao * s){
    cout << "Cidades = '";
	for (int i = 0; i < s->hubs.size(); i++){
		cout << dados->municipios[s->hubs[i]] << "','";
	}
	cout<<endl;
}

void resize(DATA * dados){

	dados->cordenadas.resize(dados->nos, vector<double>(2, 0));
	dados->distancia.resize(dados->nos, vector<celula>(dados->nos));
	dados->distanciaOrdenada.resize(dados->nos, vector<celula>(dados->nos));

	dados->demanda.resize(dados->nos, vector<double>(dados->nos, 0));
	dados->demandaDestinoOrdenada.resize(dados->nos, vector<celula>(dados->nos));
	dados->demandaOrigemOrdenada.resize(dados->nos, vector<celula>(dados->nos));
	dados->demandaOD.resize(dados->nos, vector<celula>(dados->nos));

	dados->demandaTotalOrdenada.resize(dados->nos);

	dados->O.resize(dados->nos, 0);
	dados->D.resize(dados->nos, 0);

	dados->O_ordenado.resize(dados->nos);
	dados->D_ordenado.resize(dados->nos);

    dados->populacao.resize(dados->nos, 0);
	dados->pib.resize(dados->nos, 0);
	dados->municipios.resize(dados->nos);

	dados->ligacao.resize(dados->nos, vector<int>(dados->nos, 0));

	dados->custoIntalacao.resize(dados->nos, 0);
	dados->custoInstalacaoOrdenado.resize(dados->nos);

	dados->capacidade.resize(dados->nos, 0);
	dados->outrosValores.resize(4, 0);

	dados->hubsPromissores.resize(dados->nos);
	dados->hubsPromissoresGrasp.resize(dados->nos);
}

void somaDemanda(DATA * dados){
	for (int i = 0; i < dados->nos; i++){
		for (int j = 0; j < dados->nos; j++){
			celula c3;
			c3.id = j;
			c3.valor = dados->demandaDestinoOrdenada[i][j].valor + dados->demandaOrigemOrdenada[i][j].valor;
			dados->demandaOD[i][j] = c3;
		}
	}
}

void lerArquivoDistancia(DATA * dados, char * arquivo){
    ifstream entrada;
    int aux;

	entrada.open(arquivo);

	if (!entrada.is_open()){
		cout << "Erro ao abrir o arquivo" << endl;
		entrada.close();
		system("pause");
		exit(EXIT_FAILURE);
	}

	//gravando parametros da distancia C_ij
	for (int i = 0; i<dados->nos; i++){
		for (int j = 0; j<dados->nos; j++){
			entrada >> aux;
			celula c1;
			//guarda referencia de qual no a celula se refere
			c1.id = j;
			//distancia euclidiana entre dois n�s do problema
			c1.valor = aux / 1000;
			dados->distancia[i][j] = c1;
		}
	}
	entrada.close();
}

void lerArquivoMunicipios(DATA * dados, char * arquivo){

    ifstream entrada;
    int aux;

	entrada.open(arquivo);

	if (!entrada.is_open()){
		cout << "Erro ao abrir o arquivo" << endl;
		entrada.close();
		system("pause");
		exit(EXIT_FAILURE);
	}

    string municipio;

	//gravando parametros da distancia C_ij
	for (int i = 0; i<dados->nos; i++){
        getline(entrada,municipio);
        municipio.erase(std::remove(municipio.begin(), municipio.end(), '\n'), municipio.end());
        municipio.erase(std::remove(municipio.begin(), municipio.end(), '\r'), municipio.end());
        dados->municipios[i] = municipio;
	}
	entrada.close();
}

void lerArquivoPibPopulacao(DATA * dados, char * arquivo){
    ifstream entrada;
	entrada.open(arquivo);

	if (!entrada.is_open()){
		cout << "Erro ao abrir o arquivo" << endl;
		entrada.close();
		system("pause");
		exit(EXIT_FAILURE);
	}
    int pib, populacao;
	//gravando parametros da distancia C_ij
	for (int i = 0; i<dados->nos; i++){
        entrada >> pib;
        entrada >> populacao;
        dados->pib[i] = pib;
        dados->populacao[i] = populacao;
	}
	entrada.close();
}

void lerArquivoOrdemPib(DATA * dados, char * arquivo){
    ifstream entrada;
	entrada.open(arquivo);

	if (!entrada.is_open()){
		cout << "Erro ao abrir o arquivo" << endl;
		entrada.close();
		exit(EXIT_FAILURE);
	}
    int indice;
	//gravando parametros da distancia C_ij
	for (int i = 0; i<dados->nos; i++){
        entrada >> indice;
        dados->ordemPib.push_back(indice);
	}
	entrada.close();
}
void gerarCustoFixo(DATA * dados){
    for(int i = 0; i<dados->nos; i++){
		dados->custoIntalacao.at(i) = 99999999 * float(dados->pib[i]/dados->populacao[i]);
		celula c1;
		c1.id = i;
		c1.valor = dados->custoIntalacao.at(i);
		//dados->custoInstalacaoOrdenado.at(i) = c1;
    }
}
void gerarDemanda(DATA * dados){
    long double demanda = 0;
    for(int i = 0; i<dados->nos; i++){
        for(int j = 0; j<dados->nos; j++){
        	if(i==j){dados->demanda[i][j] = 0;}
            //TODO 300 � referente a metros ou kilometros?
            double distancia = dados->distancia[i][j].valor;
            //if(distancia > 0.3){
                double popi = dados->populacao[i]/1000;
                double popj = dados->populacao[j]/1000;
                double pibi = dados->pib[i]/1000;
                double pibj = dados->pib[j]/1000;
                double expo = exp(-0.001*distancia);

                //demanda = dados->populacao[i] * dados->populacao[j] * (dados->pib[i] * dados->pib[j]);
                //demanda = (popi * popj * (pibi * pibj) * expo)/9000;
                demanda = (popi * popj * (pibi * pibj) * 1)/pow(distancia,2);
                if(demanda < 0){
                    //cout<<"Exp: "<<exp(-0.001*distancia)<<endl;
                    cout<<"ERRO: Demanda = "<<demanda<<endl;
                    exit(EXIT_FAILURE);
                }
                dados->demanda[i][j] = demanda/10000;
                //preenchendo matrizes de ordenamento
                celula c1;
                c1.id = j;
                c1.valor = demanda/10000;

                celula c2;
                c2.id = i;
                c2.valor = demanda/10000;

                dados->demandaOrigemOrdenada[i][j] = c1; //matriz demanda origem em i (W_ij)
                dados->demandaDestinoOrdenada[j][i] = c2; //matriz demanda destiino em i (W_ji)

//            }else{
//                demanda = 0;
//                dados->demanda[i][j] = demanda;
//
//                celula c1;
//                c1.id = j;
//                c1.valor = demanda;
//
//                celula c2;
//                c2.id = i;
//                c2.valor = demanda;
//
//                dados->demandaOrigemOrdenada[i][j] = c1; //matriz demanda origem em i (W_ij)
//                dados->demandaDestinoOrdenada[j][i] = c2; //matriz demanda destiino em i (W_ji)
//            }
        }
    }
}

DATA * reduzirDados(DATA * dados, int numeroCidades){
    //TODO reduzir dados
    DATA * dadosReduzido = new DATA;
    if(numeroCidades<=dados->nos){

        dadosReduzido->nos = numeroCidades;
        resize(dadosReduzido);

        dadosReduzido->alpha = dados->alpha;
        dadosReduzido->numFixoHubs = dados->numFixoHubs;
        dadosReduzido->alvo = dados->alvo;
        dadosReduzido->ordemPib = dados->ordemPib;

        for(int i = 0; i < numeroCidades; i++){
            int cidade_i = dados->ordemPib[i];

            dadosReduzido->pib[i] = dados->pib[cidade_i];
            dadosReduzido->populacao[i] = dados->populacao[cidade_i];
            dadosReduzido->municipios[i] = dados->municipios[cidade_i];
            dadosReduzido->D[i] = dados->D[cidade_i];
            dadosReduzido->O[i] = dados->O[cidade_i];

            double ex;
            if (numeroCidades <= 50) ex = 3;
            else if (numeroCidades <= 100) ex = 3.1;
            else if (numeroCidades <= 200) ex = 2.5;
            else if (numeroCidades <= 300) ex = 1.5;
            else ex = 0.7;

            dadosReduzido->custoIntalacao[i] = dados->custoIntalacao[i] * ex;
            for(int j = 0; j < numeroCidades; j++){
                int cidade_j = dados->ordemPib[j];

                dadosReduzido->distancia[i][j] = dados->distancia[cidade_i][cidade_j];
                dadosReduzido->demanda[i][j] = dados->demanda[cidade_i][cidade_j];
            }
        }
    }else{
        cout<<"ERRO: numero invalido para redimensionar instancia"<<endl;
        system("PAUSE");
    }
    return dadosReduzido;
}

void lerArquivo(DATA * dados, char * arquivo){

	ifstream entrada;

	entrada.open(arquivo);
	//entrada.open(dados->nomeArquivo); //nao ta lendo

	if (!entrada.is_open()){
		cout << "Erro ao abrir o arquivo de INSTANCIA" << endl;
		entrada.close();
		exit(EXIT_FAILURE);
	}

	entrada >> dados->nos; //recebe primeiro parametro (numero de nos)
	entrada >> dados->alpha; //economia de escala


	//define tamanho de acordo com tamanho dos nos (RESIZE)
	resize(dados);

	double aux; //variavel auxiliar para paegar valores de entrada

	double ex;
	if (dados->nos < 70) ex = 1;
	else if (dados->nos < 170) ex = 2;
	else ex = 5;

	//Seta EX apenas se for passado na leitura
    if(dados->ex != -1){
        ex = dados->ex;
    }

	//pega matriz de custo fixo de instala��o F_kk
	for (int i = 0; i < dados->nos; i++){
		entrada >> dados->custoIntalacao.at(i);
		dados->custoIntalacao.at(i) = ex * dados->custoIntalacao.at(i) / 10000;
		celula c1;
		c1.id = i;
		c1.valor = dados->custoIntalacao.at(i);

		dados->custoInstalacaoOrdenado.at(i) = c1;

		//        cout<<dados->custoIntalacao.at(i)<<endl;
	}

	//gravando paramentros de demanda dos dos W_ij
	for (int i = 0; i<dados->nos; i++){
		for (int j = 0; j<dados->nos; j++){
			entrada >> aux;
			aux = aux / 100;
			dados->demanda[i][j] = aux;


			//preenchendo matrizes de ordenamento
			celula c1;
			c1.id = j;
			c1.valor = aux;

			celula c2;
			c2.id = i;
			c2.valor = aux;


			dados->demandaOrigemOrdenada[i][j] = c1; //matriz demanda origem em i (W_ij)
			dados->demandaDestinoOrdenada[j][i] = c2; //matriz demanda destiino em i (W_ji)


		}
	}


	somaDemanda(dados);


	//gravando parametros da distancia C_ij
	for (int i = 0; i<dados->nos; i++){
		for (int j = 0; j<dados->nos; j++){
			entrada >> aux;
			aux = aux / 100;
			celula c1;
			//guarda referencia de qual no a celula se refere
			c1.id = j;
			//distancia euclidiana entre dois n�s do problema
			c1.valor = aux;
			dados->distancia[i][j] = c1;
		}
	}

	entrada.close();
}

void lerInstacias(vector < instancia > * instancias, char * arquivo){

	ifstream entrada;
	string enderecoBase, nomeArquivo;
	int qtdeIntancias;
	double FO;

	entrada.open(arquivo);

	if (!entrada.is_open()){
		cout << "Erro ao abrir o arquivo com valores OTIMO" << endl;
		entrada.close();
		exit(EXIT_FAILURE);
	}

	entrada >> qtdeIntancias;//quantidade de intancias no arquivo
	entrada >> enderecoBase; //Endere�o base da pasta onde se encotra os arquivos de dados

	//leitura de linha com otimo e endere�o do arquivo
	for (int i = 0; i < qtdeIntancias; i++){
		instancia inst;

		entrada >> nomeArquivo;
		entrada >> FO;

		string nomeFinal = enderecoBase+nomeArquivo;

		inst.arquivo = nomeFinal;
		inst.FO_otimo = FO;

		//adiciona arquivo a BASE DE INSTANCIAS
		instancias->push_back(inst);
	}

	entrada.close();

}

void imprimeDATA(DATA* dados){

	cout << "--------\nCusto Fixo" << endl;
	cout << "--------" << endl;
	for (int i = 0; i < dados->nos; i++){
		cout << dados->custoIntalacao[i] << endl;
	}
	cout << endl;
	cout << "--------\nPIB" << endl;
	cout << "--------" << endl;
	for (int i = 0; i < dados->nos; i++){
		cout << dados->pib[i] << endl;
	}
	cout << endl;


	cout << "--------\nPop" << endl;
	cout << "--------" << endl;
	for (int i = 0; i < dados->nos; i++){
		cout << dados->populacao[i] << endl;
	}
	cout << endl;

	cout << "--------\nDistancia C_ij" << endl;
	cout << "--------" << endl;
	//impress�o tempor�ria de dados
	for (int i = 0; i < dados->nos; i++){
		for (int j = 0; j < dados->nos; j++){
			cout << dados->distancia[i][j].valor << " ";
		}
		cout << endl;
	}

	cout << "--------\nDemanda Origem O_i" << endl;
	cout << "--------" << endl;
	for (int i = 0; i < dados->nos; i++){
		cout << dados->O[i] << " ";

	}
	cout << endl;


	cout << "--------\nDemanda Destino D_i" << endl;
	cout << "--------" << endl;
	for (int i = 0; i < dados->nos; i++){
		cout << dados->D[i] << " ";

	}
	cout << endl;

	cout << "--------\nDemanda W_ij" << endl;
	cout << "--------" << endl;
	for (int i = 0; i < dados->nos; i++){
		for (int j = 0; j < dados->nos; j++){
			cout << dados->demanda[i][j] << " ";
		}
		cout << endl;
	}

}

void calcula_fluxo(DATA * dados){
	//Fluxo de origem
	double somaLinha = 0;
	for (int i = 0; i<dados->nos; i++){
		for (int j = 0; j<dados->nos; j++){
			somaLinha = somaLinha + dados->demanda[i][j];
		}
		dados->O[i] = somaLinha;

		dados->O_ordenado[i].id = i;
		dados->O_ordenado[i].valor = somaLinha;
		somaLinha = 0;
	}
	//Fluxo de destino
	double somaColuna = 0;
	for (int i = 0; i<dados->nos; i++){
		for (int j = 0; j<dados->nos; j++){
			somaColuna = somaColuna + dados->demanda[j][i];
		}
		dados->D[i] = somaColuna;
		dados->D_ordenado[i].id = i;
		dados->D_ordenado[i].valor = somaColuna;
		somaColuna = 0;
	}
}


void addHub(int no, solucao * s){
    if (!isHub(no, s)){
        s->hubs.push_back(no);
        s->hubs_bin[no] = 1;
        s->alocacao[no].hub = no;
        #ifdef CYCLE_HUB
        if( s->hubs.size() > 3 ){
            tsp->solve( &(s->hubs) );
        }
        #endif // CYCLE_HUB
    }else{
        cout<<"Erro: Tentou adicionar concentrador repetido."<<endl;
        exit (EXIT_FAILURE);
	}
}

void removeHub(int no, solucao * s){
    #ifdef CYCLE_HUB
        if( s->hubs.size() <= 3 ){
            cout<<"Erro: Tentou remover concentrador em ciclo de menos de 3 nós."<<endl;
            exit (EXIT_FAILURE);
        }
    #endif // CYCLE_HUB
	if(isHub(no, s)){
		for(int i = 0; i<s->hubs.size(); i++){
			if(s->hubs[i] == no){ //busca sequencial
				s->hubs.erase(s->hubs.begin()+i); //remove da lista de hubs
				s->hubs_bin[no] = 0; //remove da lista binaria
				break;
			}
		}
		#ifdef CYCLE_HUB
        if( s->hubs.size() > 3 ){
            tsp->solve( &(s->hubs) );
        }
        #endif // CYCLE_HUB
	}else{
        cout<<"Erro: Tentou remover no concentrador nao configurado."<<endl;
        exit (EXIT_FAILURE);
        cin.ignore();
	}
}

void trocarAlocacao(int hub, int no, solucao * s){
    if(isHub(hub, s) && !isHub(no, s)){//verificando se a troca � feita entre um concentrador e um n�o concentrador
        addHub(no, s);
        removeHub(hub, s);
        for(int i = 0; i < s->hubs.size(); i++){
            if(s->hubs[i] == hub){
                s->hubs[i] = no;
                s->hubs_bin[i] = 0;
                break;
            }
        }
        s->hubs_bin[no] = 1;
        s->alocacao[no].hub = no;

        for(int i = 0; i < s->alocacao.size(); i++){//percorre todos os nos
            if(s->alocacao.at(i).hub == hub){
                s->alocacao.at(i).hub = no; //passa todas as aloca��es do hub para o no
            }
        }
    }else{
        cout<<"Erro: Troca de alocacao nao sendo concentrador."<<endl;
        exit (EXIT_FAILURE);
    }
}


void inicializaSolucao_(DATA * dados, solucao * s){
//    cout<<"Contruindo solucao GRASP..."<<endl;

	double lambida = 0.05; //Valor de aceita��o, incide na variabilidade no sorteio de L
	double FO_nova = DBL_MAX; //vairavel auxiliar para armazenar as FO

    vector<double> FO_solucoes(dados->nos, 0); //guarda a FO das solu��es finais

    solucao * s_nova = new solucao; //solu��o auxiliar para testes

    solucao * melhorSolucao = new solucao; //armazena a melhor solu��o de todas as rodadas

    double melhorFO = DBL_MAX;// FO da melhor solu��o para compara��o

    vector<solucao*> solucoes; //guarda as N solu��es geradas iniciando com o hub unico dos N nos poss�veis

	vector<double> custoMarginal; //Vetor com custo marginais da inser��o dos n�s da base testada
	double custoMarginalMax;//custo marginal m�ximo
	double custoMarginalMin;//custo marginal minimo

    vector<int> L; //vetor de n�s com melhoria da solu��o

    vector<int> NL_bin(dados->nos, 0); //vetor binario com os nos que pioraram a FO

    //constroi lista de solu��es iniciais para testes
	for (int i = 0; i<dados->nos; i++){

        solucao * s_temp = new solucao;

        s_temp->alocacao.resize(dados->nos);
        s_temp->hubs_bin.resize(dados->nos,0);
        s_temp->hubs.resize(0);

        addHub(i,s_temp);//cada uma inicia com um hub diferente

        for (int j = 0; j < dados->nos; j++){
            s_temp->alocacao[j].id = j; //indice do no
            //s_temp->alocacao[j].hub = dados->hubsPromissores[i].id;//todos os n�s alocados ao primeiro da lista
            s_temp->alocacao[j].hub = i;//todos os n�s alocados ao primeiro da lista
        }

        solucoes.push_back(s_temp);  //aloca espa�o no velor
        FO_solucoes.at(i) = Calcula_FO(dados, s_temp); //pega o valor da solu��o inicial

	}

    for(int n = 0; n<solucoes.size(); n++){ //Controla a possibilidade de cada no come�ar como solu��o

        //controi solucao com unico hub, sequenciamente setado por N
        //inicializa com todos os nos alocados ao primeiro hub da lista de promissores

        //aloca todo mundo no hub inicial

        //zera NL e Inclui os hubs em L (apenas para entrar no while, dentro do processo o L � reconstruido)
        custoMarginal.clear();
        custoMarginal.resize(dados->nos, 0); //zera custo marginal
        for(int i = 0; i < dados->nos; i++){
            NL_bin[i] = 0;
        }
        //enquanto tiver solu��es em L
        do{
            custoMarginalMax = -DBL_MAX;
            #ifdef CYCLE_HUB
            custoMarginalMin = DBL_MAX;
            #else
            custoMarginalMin = 0;
            #endif // CYCLE_HUB
            for(int j = 0; j < dados->nos; j++){ //percorre os N nos a serem testados
                custoMarginal[j] = DBL_MAX;
                *s_nova = *solucoes.at(n); //armazena a solu��o inicial em variavel temporaria
                if(!isHub(j, s_nova) && NL_bin[j] == 0){// se n�o for hub da solu��o e n�o estar em NL

                    addHub(j, s_nova); //adiciona hub na solu��o temporaria
                    alocarNos(dados, s_nova); //aloca nos ao hub mais proximo
                    FO_nova = Calcula_FO(dados, s_nova);
                    custoMarginal[j] = FO_nova - FO_solucoes.at(n); //calculo do custo marginal
                    if(custoMarginal[j] < 0 || solucoes.at(n)->hubs.size() < 3){
                        if(custoMarginal[j] < custoMarginalMin){
                            custoMarginalMin = custoMarginal[j];
                        }
                        if(custoMarginal[j] > custoMarginalMax){
                            custoMarginalMax = custoMarginal[j];
                        }
                        NL_bin[j] = 0;
                    }else{
                        NL_bin[j] = 1;//caso houver piora, ele � adicionado ao NL
                    }
                }
            }
            //Limpa a lista L para cria��o de nova lista a partir dos custos marginais capiturados
            L.clear();
            #ifdef CYCLE_HUB
            if(custoMarginalMin < 0 || solucoes.at(n)->hubs.size() < 3){
            #else
            if(custoMarginalMin < 0){//se houve melhoria em algum caso, procede
            #endif // CYCLE_HUB
                double margem = custoMarginalMin + lambida*(custoMarginalMax - custoMarginalMin); //calcula a margem para compara��o
                for(int j = 0; j < dados->nos; j++){
                    #ifdef CYCLE_HUB
                    if((custoMarginal[j] < 0 || solucoes.at(n)->hubs.size() < 3) && (custoMarginal[j] <= margem)){//se houve melhoria ou tiver menos que 3 hubs
                    #else
                    if((custoMarginal[j] < 0) && (custoMarginal[j] <= margem)){//se houver melhoria dentro da margem
                    #endif // CYCLE_HUB
                        L.push_back(j);//nesse caso � incluso nos hubs para sorteio
                    }
                }
            }

            if(L.size()>0){//se a lista n�o estiver vazia, � feito um sorteio dos itens dentro da margem
                int tamanhoL = L.size();
                int noSorteiado = sorteioPosicao(0,L.size()); //sorteio da posi��o (0 a N-1)
                addHub(L.at(noSorteiado), solucoes.at(n));
                alocarNos(dados, solucoes.at(n));
                FO_solucoes[n] = Calcula_FO(dados, solucoes.at(n)); //atualiza a FO da solu��o

            }
            //caso ela for a melhor at� o momento, ela � armazenada
            if(FO_solucoes.at(n) < melhorFO){
                melhorFO = FO_solucoes[n];
                *melhorSolucao = *solucoes.at(n);
            }
        }while(L.size() > 0);
    }
    *s = *melhorSolucao;
    #ifdef CYCLE_HUB
    while(s->hubs.size() < 3){
        pertubar_AdicionaHub(dados, s);
    }
    #endif // CYCLE_HUB
    for(int i = 0; i<solucoes.size(); i++){
        delete solucoes[i];
    }
    delete s_nova;
    delete melhorSolucao;
}

void imprimeVetor(vector< vector<double> > v, DATA * dados){
	for (int i = 0; i < dados->nos; i++){
		for (int j = 0; j < v[i].size(); j++){
			cout << v[i][j] << " ";
		}
		cout << endl;
	}
}
void imprimeVetor(vector< vector<int> > v, DATA * dados){
	for (int i = 0; i < dados->nos; i++){
		for (int j = 0; j < v[i].size(); j++){
			cout << v[i][j] << " ";
		}
		cout << endl;
	}
}

//metodo padrao auxiliar de ordenamento do QuickSort (Sub processo de ordenamento)
int partition(vector<celula>& vetor, int inicio, int fim){
	celula pivot = vetor[fim];
	int bottom = inicio - 1;
	int top = fim;

	bool notdone = true;
	while (notdone){
		while (notdone){
			bottom += 1;
			if (bottom == top){
				notdone = false;
				break;
			}
			if (vetor[bottom].valor > pivot.valor){
				vetor[top] = vetor[bottom];
				break;
			}
		}
		while (notdone){
			top = top - 1;
			if (top == bottom){
				notdone = false;
				break;
			}
			if (vetor[top].valor < pivot.valor){
				vetor[bottom] = vetor[top];
				break;
			}
		}
	}
	vetor[top] = pivot;
	return top;
}
int sort(vector<celula>& vetor, int inicio, int fim){
	if (inicio < fim){
		int split = partition(vetor, inicio, fim);   //recursion
		sort(vetor, inicio, split - 1);
		sort(vetor, split + 1, fim);
	}
	else{
		return 0;
	}
}

void ordenaDistancia(DATA * dados){ //oredena matriz de distancia euclidiana lida de arquivo
	dados->distanciaOrdenada = dados->distancia;//faz uma copia da matriz de distancia original antes de ordenar
	for (int i = 0; i<dados->nos; i++){//la�o para percorrer todas as linhas da matriz
		sort(dados->distanciaOrdenada[i], 0, dados->nos - 1);//ordena linha por linha da matriz, as linhas s�o compostas de celulas que carregam seu indice inicial
	}
}

void ordenaPenalidade(DATA * dados){
	sort(dados->hubsPromissores, 0, dados->nos - 1);
}

void ordenaPromissoresGRASP(DATA * dados){
    sort(dados->hubsPromissoresGrasp, 0, dados->nos - 1);
}

void ordenaCustoInstalacao(DATA * dados){
	sort(dados->custoInstalacaoOrdenado, 0, dados->nos - 1);
}

void ordenaDemanda(DATA * dados){
	for (int i = 0; i<dados->nos; i++){
		sort(dados->demandaDestinoOrdenada[i], 0, dados->nos - 1);
		sort(dados->demandaOrigemOrdenada[i], 0, dados->nos - 1);
		sort(dados->demandaOD[i], 0, dados->nos - 1);
	}
	sort(dados->O_ordenado, 0, dados->nos - 1);
	sort(dados->D_ordenado, 0, dados->nos - 1);
}

void penalizaNos(DATA * dados){

	for (int i = 0; i<dados->nos; i++){
		dados->hubsPromissores[i].id = i;
		dados->hubsPromissores[i].valor = 0;
	}

	for (int i = 0; i<dados->nos; i++){
		dados->hubsPromissores[dados->custoInstalacaoOrdenado[i].id].valor += i * dados->peso_cf;
		dados->hubsPromissores[dados->O_ordenado[i].id].valor += i * dados->peso_od;
		dados->hubsPromissores[dados->D_ordenado[i].id].valor += i * dados->peso_od;
		for (int j = 0; j<dados->nos; j++){
			//dados->hubsPromissores[dados->demandaOD[i][j].id].valor += j;
			dados->hubsPromissores[dados->distanciaOrdenada[i][j].id].valor += j * dados->peso_dist;
			//dados->hubsPromissores[dados->distanciaOrdenada[i][j].id].penalidade += j;
			//dados->hubsPromissores[dados->demandaOrigemOrdenada[i][j].id].valor += j;
			//dados->hubsPromissores[dados->demandaDestinoOrdenada[i][j].id].valor += j;
		}
	}
}

void pertubar_RemoveHub(DATA * dados, solucao * s){
    #ifdef CYCLE_HUB
    if(s->hubs.size() > 3){
    #else
    if(s->hubs.size() > 1){
    #endif // CYCLE_HUB
        int hub_sorteio = s->hubs.at( sorteioPosicao(0,s->hubs.size()) ); //seleciona um hub aleatorio
        removeHub(hub_sorteio, s);//remove hub sorteado
        alocarNos(dados, s);//aloca ao mais proximo
    }
}
void pertubar_AdicionaHub(DATA * dados, solucao * s){
    int no_sorteio;
    do{
        no_sorteio = sorteioPosicao(0,dados->nos); //seleciona um no aleatorio
    }while(isHub(no_sorteio,s));//se j� for um hub sorteia denovo
    addHub(no_sorteio, s);//adiciona no sorteado
    alocarNos(dados, s);//aloca ao mais proximo
}
void pertubar_trocaFuncao(DATA * dados, solucao * s){
    int hub_sorteio = s->hubs.at( sorteioPosicao(0,s->hubs.size()) ); //seleciona um hub aleatorio
    int no_sorteio;
    do{
        no_sorteio = sorteioPosicao(0,dados->nos); //seleciona um no aleatorio
    }while(isHub(no_sorteio, s));//se j� for um hub sorteia denovo
    trocarAlocacao(hub_sorteio,no_sorteio,s);
}
void pertubar_trocaAlocacao(DATA * dados, solucao * s){
    if(s->hubs.size() > 1){
        int no_sorteio;
        do{
            no_sorteio = sorteioPosicao(0,dados->nos); //seleciona um no aleatorio
        }while(isHub(no_sorteio, s));//se j� for um hub sorteia denovo

        //seleciona novo hub para ser alocado ao n� escolhido
        int hub_sorteio;
        do{
            hub_sorteio = s->hubs.at( sorteioPosicao(0,s->hubs.size()) );
        }while(s->alocacao[no_sorteio].hub == hub_sorteio);//nao pode ser o mesmo hub ja alocado
        //executa a realoca��o
        s->alocacao[no_sorteio].hub = hub_sorteio;
    }
}

void pertubar_JucaoSolucoes(DATA * dados, solucao * s){
//    cout<<"Contruindo solucao GRASP..."<<endl;

	double lambida = 0.05; //Valor de aceita��o, incide na variabilidade no sorteio de L
	double FO_nova = 1000000000; //vairavel auxiliar para armazenar as FO

    vector<double> FO_solucoes(dados->nos, 0); //guarda a FO das solu��es finais

    solucao * s_nova = new solucao; //solu��o auxiliar para testes

    solucao * melhorSolucao = new solucao; //armazena a melhor solu��o de todas as rodadas

    double melhorFO = 1000000000;// FO da melhor solu��o para compara��o

    vector<solucao*> solucoes; //guarda as N solu��es geradas iniciando com o hub unico dos N nos poss�veis

	vector<double> custoMarginal; //Vetor com custo marginais da inser��o dos n�s da base testada
	double custoMarginalMax;//custo marginal m�ximo
	double custoMarginalMin;//custo marginal minimo

    vector<int> L; //vetor de n�s com melhoria da solu��o

    vector<int> NL_bin(dados->nos, 0); //vetor binario com os nos que pioraram a FO

    vector<int> nosSorteio;

    for(int i = 0; (i<dados->nos) && (nosSorteio.size() < dados->nos/2); i++){
        if(!isHub(dados->hubsPromissores[i].id, s)){
            nosSorteio.push_back(dados->hubsPromissores[i].id);
            addHub(dados->hubsPromissores[i].id, s);
        }
    }

    //Faz uma lista de solu��es removendo cada uma dos hubs por vez
	for (int i = 0; i < s->hubs.size(); i++){
        solucao * s_temp = new solucao;
        *s_temp = *s;
        removeHub(s_temp->hubs[i], s_temp);
        alocarNos(dados, s_temp);
        solucoes.push_back(s_temp);  //aloca espa�o no velor
        FO_solucoes.at(i) = Calcula_FO(dados, s_temp); //pega o valor da solu��o inicial
	}

    for(int n = 0; n<solucoes.size(); n++){ //Controla a possibilidade de cada no come�ar como solu��o


        //controi solucao com unico hub, sequenciamente setado por N
        //inicializa com todos os nos alocados ao primeiro hub da lista de promissores

        //aloca todo mundo no hub inicial

        //zera NL e Inclui os hubs em L (apenas para entrar no while, dentro do processo o L � reconstruido)
        custoMarginal.clear();
        custoMarginal.resize(dados->nos, 0); //zera custo marginal
        for(int i = 0; i < dados->nos; i++){
            NL_bin[i] = 0;
        }
        //enquanto tiver solu��es em L
        do{
            custoMarginalMax = -1000000000;
            custoMarginalMin = 0;
            for(int j = 0; j < dados->nos; j++){ //percorre os N nos a serem testados
                custoMarginal[j] = 0;
                *s_nova = *solucoes.at(n); //armazena a solu��o inicial em variavel temporaria
                if(isHub(j, s_nova) && NL_bin[j] == 0){// se for hub da solu��o e n�o estar em NL

                    removeHub(j, s_nova); //remove hub da solu��o temporaria
                    alocarNos(dados, s_nova); //aloca nos ao hub mais proximo
                    FO_nova = Calcula_FO(dados, s_nova);
//                    imprimeSolucao(s_nova);
//                    cout<<"FO Nova = "<<FO_nova<<endl;
                    custoMarginal[j] = FO_nova - FO_solucoes.at(n); //calculo do custo marginal
                    if(custoMarginal[j] < 0){
                        NL_bin[j] = 0;
                        if(custoMarginal[j] < custoMarginalMin){
                            custoMarginalMin = custoMarginal[j];
                        }
                        if(custoMarginal[j] > custoMarginalMax){
                            custoMarginalMax = custoMarginal[j];
                        }
                    }else{
                       NL_bin[j] = 1;//caso houver piora, ele � adicionado ao NL
                    }
                }
            }
            //Limpa a lista L para cria��o de nova lista a partir dos custos marginais capiturados
            L.clear();
            if(custoMarginalMin < 0){//se houve melhoria em algum caso, procede
                double margem = custoMarginalMin + lambida*(custoMarginalMax - custoMarginalMin); //calcula a margem para compara��o
//                cout<<"MARGEM = "<<margem<<endl;
//                cout<<"CM Min = "<<custoMarginalMin<<endl;
//                cout<<"CM Max = "<<custoMarginalMax<<endl;
                for(int j = 0; j < dados->nos; j++){
//                    cout<<"No "<<j<<" = "<<custoMarginal[j]<<endl;
                    if((custoMarginal[j] < 0) && (custoMarginal[j] <= margem)){//se houve melhoria e est� dentro da margem
//                      cout<<"No "<<j<<" aceito"<<endl;
                        L.push_back(j);//nesse caso � incluso nos hubs para sorteio
                    }
                }
            }

            if(L.size()>0){//se a lista n�o estiver vazia, � feito um sorteio dos itens dentro da margem
                int tamanhoL = L.size();
                int noSorteiado = sorteioPosicao(0,L.size()); //sorteio da posi��o (0 a N-1)
                removeHub(L.at(noSorteiado), solucoes.at(n));
                alocarNos(dados, solucoes.at(n));
                FO_solucoes[n] = Calcula_FO(dados, solucoes.at(n)); //atualiza a FO da solu��o
            }
            //caso ela for a melhor at� o momento, ela � armazenada
            if(FO_solucoes.at(n) < melhorFO){
                melhorFO = FO_solucoes[n];
                *melhorSolucao = *solucoes.at(n);
            }
            //if (FO_solucoes[n] <= dados->alvo + 0.1) {
//                cout<<"Alcancou otimo"<<endl;
             //   break;
            //}
//            imprimeSolucao(solucoes[n]);
//            cout<<"FO = "<<FO_solucoes[n]<<endl;
//            cout<<"MELHOR = "<<melhorFO<<endl;
//            cin.ignore();
        }while(L.size() > 0);
        //if (melhorFO <= dados->alvo + 0.1) {
//            cout<<"Alcancou otimo"<<endl;
            //break;
        //}

    }
    *s = *melhorSolucao;
    //printf("Custo (FO) = %18.4f \n", melhorFO);
//    imprimeSolucao(s);
//    ordenaPromissoresGRASP(dados);
//    for (int i = 0; i<6; i++){
//        int indice = dados->nos - 1 - i;
//        cout<<"No "<<dados->hubsPromissoresGrasp[indice].id<<" - "<<dados->hubsPromissoresGrasp[indice].valor<<endl;
//	}
//	cin.ignore();
    for(int i = 0; i<solucoes.size(); i++){
        delete solucoes[i];
    }
    delete s_nova;
    delete melhorSolucao;
}

void VNS_GRASP(DATA * dados, solucao * s, solucao * s_star, double * FO_star){
    solucao * s_nova = new solucao;
	double FO_depois;
	double FO_antes = Calcula_FO(dados, s);
	int numIteracoes = 0;

	while(numIteracoes<5){
        *s_nova = *s;
        pertubar_JucaoSolucoes(dados, s_nova);//faz jun��o de solu��es
        VND(dados, s_nova, s_star, FO_star); //aplica busca local

        FO_depois = Calcula_FO(dados, s_nova);//calcula nova FO

        if (FO_depois < FO_antes) {//verifica se houve melhoria
            //guarda nova FO
            FO_antes = FO_depois;
            //nova solucao passa ser a atual
            *s = *s_nova;
            numIteracoes = 0;
            if(FO_depois < *FO_star){
                *FO_star = FO_depois;
                *s_star = *s_nova;
            }
        }else{
            numIteracoes++;
        }
	}
	*s = *s_star;
	delete s_nova;
}

void VNS_TROCA(DATA * dados, solucao * s, solucao * s_star, double * FO_star){
    solucao * s_nova = new solucao;
	double FO_depois;
	double FO_antes = Calcula_FO(dados, s);
	int numIteracoes = 0;
	while(numIteracoes<20){
        *s_nova = *s;
        pertubar_trocaFuncao(dados, s_nova);//Troca fun��o de n�s
        VND(dados, s_nova, s_star, FO_star); //aplica busca local

        FO_depois = Calcula_FO(dados, s_nova);//calcula nova FO
        if (FO_depois < FO_antes) {//verifica se houve melhoria
            //guarda nova FO
            FO_antes = FO_depois;
            //nova solucao passa ser a atual
            *s = *s_nova;
            if(FO_depois < *FO_star){
                *FO_star = FO_depois;
                *s_star = *s_nova;
            }
            numIteracoes = 0;
        }else{
            VNS(dados, s_nova, s_star, FO_star);
            FO_depois = Calcula_FO(dados, s_nova);
            if (FO_depois < FO_antes) {//verifica se houve melhoria
                //guarda nova FO
                FO_antes = FO_depois;
                //nova solucao passa ser a atual
                *s = *s_nova;
                if(FO_depois < *FO_star){
                    *FO_star = FO_depois;
                    *s_star = *s_nova;
                }
                numIteracoes = 0;
            }else{
                numIteracoes++;
            }
        }
        if (Calcula_FO(dados, s_nova) <= dados->alvo + 0.1) {
            break;
        }
	}
	*s = *s_star;
	delete s_nova;
}

void ILS(DATA * dados, solucao * s, solucao * s_star, double * FO_star, int maxInter){

    //Variáveis de apoio
    solucao * s_nova = new solucao;
    double FO_depois;
	double FO_antes = Calcula_FO(dados, s);

    if(maxInter <= 0){
        cout << "ILS precisa de no minimo 1 iteracao para executar" << endl;
        exit (EXIT_FAILURE);
    }
    for(int i = 0; i < maxInter; i++){
        //Realiza sorteio para saber qual pertubação utilizar
        int pertubacao = rand() % 3;

        //Efetua copia temporária de solução
        *s_nova = *s;
        FO_antes = Calcula_FO(dados, s_nova);

        //aplica pertubação sorteada
        switch(pertubacao){
            case 0:
                pertubar_trocaFuncao(dados, s_nova);
                pertubar_trocaFuncao(dados, s_nova);
                pertubar_RemoveHub(dados, s_nova);
            break;
            case 1:
                pertubar_trocaFuncao(dados, s_nova);
                pertubar_trocaFuncao(dados, s_nova);
                pertubar_trocaFuncao(dados, s_nova);
                pertubar_trocaFuncao(dados, s_nova);
                pertubar_trocaFuncao(dados, s_nova);
                pertubar_AdicionaHub(dados, s_nova);
            break;
            case 2:
                pertubar_trocaFuncao(dados, s_nova);
                pertubar_trocaFuncao(dados, s_nova);
                pertubar_trocaFuncao(dados, s_nova);
                pertubar_trocaFuncao(dados, s_nova);
                pertubar_trocaFuncao(dados, s_nova);
                pertubar_trocaFuncao(dados, s_nova);
                pertubar_trocaFuncao(dados, s_nova);
                pertubar_AdicionaHub(dados, s_nova);
                pertubar_AdicionaHub(dados, s_nova);
                pertubar_AdicionaHub(dados, s_nova);
                pertubar_RemoveHub(dados, s_nova);
                pertubar_RemoveHub(dados, s_nova);
            break;
        }

        //Aplica busca local
        VND(dados, s_nova, s_star, FO_star);

        FO_depois = Calcula_FO(dados, s_nova);//calcula nova FO
        if (FO_depois < FO_antes) {//verifica se houve melhoria
            //guarda nova FO
            FO_antes = FO_depois;
            //nova solucao passa ser a atual
            *s = *s_nova;
        }
    }
    if( Calcula_FO(dados, s) < Calcula_FO(dados, s_star) ){
        *FO_star = Calcula_FO(dados, s_star);
        *s_star = *s;
    }
}

void VNS(DATA * dados, solucao * s, solucao * s_star, double * FO_star){
	int vizinhanca = 1;
	solucao * s_nova = new solucao;
	double FO_depois;
	double FO_antes = Calcula_FO(dados, s);
	while (vizinhanca <= 4) { //explorando todas vizinhancas
		switch (vizinhanca) {
            case 1:
                //cout<<"VNS - Troca Alocacao"<<endl;
                *s_nova = *s; //reseta solucao
                FO_antes = Calcula_FO(dados, s_nova);
                pertubar_trocaAlocacao(dados, s_nova);//Troca uma aloca��o aleat�ria
                VND(dados, s_nova, s_star, FO_star); //aplica busca local
                FO_depois = Calcula_FO(dados, s_nova);//calcula nova FO

                if (FO_depois < FO_antes) {//verifica se houve melhoria

                    //guarda nova FO
                    FO_antes = FO_depois;
                    //nova solucao passa ser a atual
                    *s = *s_nova;

                    //se chegar na FO otima para
                    #ifdef STOP_OPTIMAL
                    if (FO_depois <= dados->alvo + 0.1) {
                        cout<<"ALCANCOU OTIMO"<<endl;
                        vizinhanca = 5;
                        break;
                    }
                    #endif // STOP_OPTIMAL
                    vizinhanca = 1;
                }else{
                    vizinhanca = 2;
                }
            case 2:
                //cout<<"VNS - Troca Funcao"<<endl;
                *s_nova = *s; //reseta solucao
                FO_antes = Calcula_FO(dados, s_nova);
                pertubar_trocaFuncao(dados, s_nova);//Troca fun��o de n�s
                VND(dados, s_nova, s_star, FO_star); //aplica busca local
                FO_depois = Calcula_FO(dados, s_nova);//calcula nova FO

                if (FO_depois < FO_antes) {//verifica se houve melhoria

                    //guarda nova FO
                    FO_antes = FO_depois;
                    //nova solucao passa ser a atual
                    *s = *s_nova;



                    //se chegar na FO otima para
                    #ifdef STOP_OPTIMAL
                    if (FO_depois <= dados->alvo + 0.1) {
                        cout<<"ALCANCOU OTIMO"<<endl;
                        vizinhanca = 5;
                        break;
                    }
                    #endif // STOP_OPTIMAL
                    vizinhanca = 1;
                }else{
                    vizinhanca = 3;
                }
            case 3:
                //cout<<"VNS - Adiciona concentrador"<<endl;
                *s_nova = *s; //reseta solucao
                FO_antes = Calcula_FO(dados, s_nova);
                pertubar_AdicionaHub(dados, s_nova); //pertuba
                VND(dados, s_nova, s_star, FO_star); //aplica busca local

                FO_depois = Calcula_FO(dados, s_nova);//calcula nova FO

                if (FO_depois < FO_antes) {//verifica se houve melhoria

                    //guarda nova FO
                    FO_antes = FO_depois;
                    //nova solucao passa ser a atual
                    *s = *s_nova;



                    //se chegar na FO otima para
                    #ifdef STOP_OPTIMAL
                    if (FO_depois <= dados->alvo + 0.1) {
                        cout<<"ALCANCOU OTIMO"<<endl;
                        vizinhanca = 5;
                        break;
                    }
                    #endif // STOP_OPTIMAL
                    vizinhanca = 1;
                }else{
                    vizinhanca = 4;
                }
            case 4:
                //cout<<"VNS - Remove Hub"<<endl;
                *s_nova = *s; //reseta solucao
                FO_antes = Calcula_FO(dados, s_nova);

                pertubar_RemoveHub(dados, s_nova); //pertuba

                VND(dados, s_nova, s_star, FO_star); //aplica busca local

                FO_depois = Calcula_FO(dados, s_nova);//calcula nova FO

                if (FO_depois < FO_antes) {//verifica se houve melhoria

                    //guarda nova FO
                    FO_antes = FO_depois;
                    //nova solucao passa ser a atual
                    *s = *s_nova;
                    //se chegar na FO otima para
                    #ifdef STOP_OPTIMAL
                    if (FO_depois <= dados->alvo + 0.1) {
                        cout<<"ALCANCOU OTIMO"<<endl;
                        vizinhanca = 5;
                        break;
                    }
                    #endif // STOP_OPTIMAL
                    vizinhanca = 1;
                }else{
                    vizinhanca = 5;
                }
		}
		#ifdef STOP_OPTIMAL
		if (FO_depois <= dados->alvo + 0.1) {
            cout<<"ALCANCOU OTIMO"<<endl;
            break;
        }
        #endif // STOP_OPTIMAL
	}
	delete s_nova;
}

/**
*   TSP - Travelling salesman problem
*   Authors: Rodrigo Brito (contato@rodrigobrito.net) e Bruno Gomes (bruno.nonato@ifmg.edu.br)
*   Project URI: https://github.com/rodrigo-brito/C-Vector-TSP
*/

#include "tsp.h"

TSP::TSP(vector< vector<celula> > * d) {
    distance = d;
}

void TSP::solve( vector<int> * tour ) {
    if(tour->size() > 4 ){//TSP only for more than 4 elements
        int rval = 0;
        int semente = rand();
        double szeit, optval, *mybnd, *mytimebound;
        int ncount, success, foundtour, hit_timebound = 0;
        int *in_tour = (int *) NULL;
        int *out_tour = (int *) NULL;
        CCrandstate rstate;
        char *probname = (char *) NULL;
        static int run_silently = 1;
        CCutil_sprand(semente, &rstate);
        mybnd = (double *) NULL;
        mytimebound = (double *) NULL;
        ncount = tour->size();
        int ecount = (ncount * (ncount - 1)) / 2; //Total de arestas distintas
        int *elist = new int[ecount * 2];
        int *elen = new int[ecount];
        int edge = 0;
        int edge_peso = 0;
        for (int kk = 0; kk < tour->size(); kk++) { //Percorre todas as combinações de linhas
            for (int m = kk + 1; m < tour->size(); m++) { //Por colunas
                if (kk != m) { //Desconsidera a diagonal
                    elist[edge] = tour->at(kk);
                    elist[edge + 1] = tour->at(m);
                    elen[edge_peso]	= distance->at( elist[edge] )[ elist[edge + 1] ].valor;
                    edge_peso++;
                    edge = edge + 2;
                }
            }
        }
        out_tour = CC_SAFE_MALLOC (ncount, int);
        probname = CCtsp_problabel(" ");
        rval = CCtsp_solve_sparse(ncount, ecount, elist, elen, in_tour,
                out_tour, mybnd, &optval, &success, &foundtour, probname,
                mytimebound, &hit_timebound, run_silently, &rstate);
        for (int kk = 0; kk < ncount; kk++) {
            tour->at(kk) = out_tour[kk];
        }
        szeit = CCutil_zeit();
        CC_IFFREE (elist, int);
        CC_IFFREE (elen, int);
        CC_IFFREE (out_tour, int);
        CC_IFFREE (probname, char);
    }else if( tour->size() == 4 ){
        //Tour for 4 elements
        minTour(tour);
    }else{
        cout << "WARNING: TSP Class Don't work for less than 4 elements" << endl;
    }
}

double TSP::getCost( vector<int> * tour ){
    double cost = distance->at(tour->at( tour->size()-1 ))[ tour->at(0) ].valor; //Inicia com o custo da última posição voltando para o início (Fechamento do Tour)
    for(int i = 0; i < tour->size()-1; i++){ //Percorre as demais posições
        cost += distance->at( tour->at(i) )[ tour->at(i+1) ].valor;
    }
    return cost;
}

void TSP::minTour( vector<int> * tour ){
    if(tour->size() == 4){
        double min_cost = getCost(tour);
        vector<int> min_tour = *tour;

        //Possibilidades de rotas
        for(int i = 1; i < 3; i++){
            vector<int> temp_tour = *tour;
            int temp = temp_tour.at(i);
            temp_tour.at(i) = temp_tour.at(i+1);
            temp_tour.at(i+1) = temp;

            double cost = getCost( &temp_tour );
            if(cost < min_cost){
                min_cost = cost;
                min_tour = temp_tour;
            }
        }
        *tour = min_tour;
    }else{
        cout << "WARNING: Min Tour used only for 4 elements" << endl;
    }
}

double TSP::getSubTourDistance( vector<int> * tour, int origin, int destiny ){
    double distance_total = 0;
    int index_origin = -1;
    for (int i = 0; i < tour->size(); i++) {
        if(tour->at(i) == origin){
            index_origin = i;
            break;
        }
    }
    if( index_origin != -1 ){
        int index_destiny = -1;
        for (int i = 0; i < tour->size(); i++) {
            int index_search = (i + index_origin) % tour->size();//Posição atual, % faz indice ter referencia circular
            int nex_position = (i + index_origin + 1) % tour->size();//Seguindo o círculo, qual a proxima posição
            if( tour->at(index_search) != destiny ){//Enquanto não encontrar o destino, soma a distancia até o próximo nó
                distance_total += distance->at( tour->at(index_search) )[ tour->at(nex_position) ].valor;
            }else{ //Se encontrou para
                index_destiny = index_search;
                break;
            }
        }
        if( index_destiny == -1 ){//Se retornar -1 não encontrou nada, ERRO de inconsistência de busca
            cout<<"ERRO: TSP::getSubTourDistance(search) DEST "<<destiny<<" not found!"<<endl;
            exit (EXIT_FAILURE);
        }
    }else{//Se retornar -1 não encontrou nada, ERRO de inconsistência de busca
        cout<<"ERRO: TSP::getSubTourDistance(search) ORG "<<origin<<" not found!"<<endl;
        exit (EXIT_FAILURE);
    }
    return distance_total;
}

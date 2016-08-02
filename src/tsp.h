#ifndef TSP_H
#define TSP_H

#include <iostream>
#include <vector>
#include <qsopt.h>
#include "structs.h"
extern "C" {
	#include <concorde.h>
}

using namespace std;

class TSP {
    private:
        vector< vector<celula> > * distance;
        vector<int> * tour;
    public:
        TSP(vector< vector<celula> > * distance);
        vector<int> * solve( vector<int> tour );
        void solve( vector<int> * tour );
        double getCost( vector<int> * tour );
        void minTour( vector<int> * tour );
        double getSubTourDistance( vector<int> * tour, int origin, int destiny );
};
#endif

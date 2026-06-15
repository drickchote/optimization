#ifndef BOUNDEDPARETOSET_CPP
#define BOUNDEDPARETOSET_CPP

#include "ParetoSet.cpp"
#include <cstdio>
#include <string>


using namespace std;

class BoundedParetoSet : public ParetoSet {
	private:
	const static int MAXARCSIZE = 1000;
	// FILE *globalf;
	// bool existeFileGlobal;
    // string nomeglobalf;
	
	public:
	/* Complexidade: O(n) */
	BoundedParetoSet () {
		// existeFileGlobal = false;
	}
	~BoundedParetoSet () {
		// if (existeFileGlobal) {
		// 	fclose(globalf);
		// }
	}

    // void setNomeGlobalf (char * n) {
    //     nomeglobalf = n;
    // }

	bool adicionarSol(SolucaoEdgeSet *s) {
		ASS ( assert( confereGrid() ); )
		list<SolucaoEdgeSet *>::iterator maisPopuloso = sol.begin();
		int maiorPositionCount = -1;
		if (sol.size() > 0) maiorPositionCount = getPositionCount( **sol.begin() );
		
		/* percorre o vetor de solucoes e de valores e, caso exista solucao dominada, retira e retorna true. caso contrario, retorna false */
		list<SolucaoEdgeSet *>::iterator i = sol.begin();
		list< list<SolucaoEdgeSet *>::iterator > remover;
		while (i != sol.end()) {
			// se a solucao que vai entrar domina a solucao i que ja se encontra no conjunto
			if (*s >> **i) {
				remover.push_back(i);
			}
			// se a solucao que vai entrar nao domina a solucao i, procuramos a solucao que esta no
			// local mais populoso, caso o pareto esteja no tamanho maximo
			if (remover.size() == 0 && getSize()+1 > MAXARCSIZE) {
				int positionCountAtual = getPositionCount(**i);
				if (maiorPositionCount < positionCountAtual) {
					maiorPositionCount = positionCountAtual;
					maisPopuloso = i;
				}
			}
			
			if (**i >> *s || **i == *s)
				return false;
			i++;
		}

		// se a solucao que vai entrar nao domina nenhuma e o tamanho do conjunto pareto ja esta no maximo
		// (se nenhuma solucao vai sair do conjunto), remove a mais populosa
		if (remover.size() == 0 && getSize()+1 > MAXARCSIZE) {
			// //fprintf(stderr,"BoundedPareto: removendo a solucao mais populosa de pop = %d\n",maiorPositionCount);
			// if (!existeFileGlobal) {
			// 	// globalf = fopen (nomeglobalf.data(), "a");
			// 	existeFileGlobal = true;
			// }
			// else {
			// 	if (globalf == NULL)
			// 		cerr << "erro ao abrir arquivo!!" << endl;
			// 	else
			// 		printOnePoint (maisPopuloso, globalf);
			// }
			remover.push_back(maisPopuloso);
		}

		//fprintf(stderr,"getSize = %d %d\n",getSize(),sol.size());
			
		list< list<SolucaoEdgeSet *>::iterator >::iterator j = remover.begin();
		while (j != remover.end()) {
			// remove do grid
			g.removeGrid( calcularGridPos(***j) );
			// remove a frequencia das arestas
			removerFrequencia( **j );
			delete( **j );
			// remove do conjunto pareto
			sol.erase( *j );
			j++;
		}

		SolucaoEdgeSet *t = new SolucaoEdgeSet(s->nEdges,*(s->rg));
		*t = *s;
		// adiciona ao conjunto pareto
		sol.push_front( t );
		if (sol.size() > MAXARCSIZE) fprintf(stderr,"ERRO!\n");
		// adiciona a frequencia das arestas da solucao
		adicionarFrequencia( t );
		// adiciona ao grid
		g.addGrid( calcularGridPos(*t) );

		for(int k=0;k<NUMOBJETIVOS;k++) {
			rangeNovo[k].min = min(rangeNovo[k].min,t->getObj(k));
			rangeNovo[k].max = max(rangeNovo[k].max,t->getObj(k));
		}

		// se houve uma mudanca grande nos ranges (maior que 10% valor), atualizar o grid
		for (int k=0;k<NUMOBJETIVOS;k++) {
			if (fabs(rangeNovo[k].min-rangeAtual[k].min) > 0.1*rangeAtual[k].min || fabs(rangeNovo[k].max-rangeAtual[k].max) > 0.1*rangeAtual[k].max) {
				//fprintf(stderr,"Atualizando grid!\n");
				updateGrid();
				break;
			}
		}

		ASS ( assert( confereGrid() ); )
		return true;
	}
};

#endif
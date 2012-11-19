#include <boost/mpi.hpp>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/assume_abstract.hpp>

#include <boost/optional.hpp>

#include <stdexcept>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <assert.h>

#include <thread>
#include <mutex>
#include <atomic>
#include <queue>

#include <eo>
#include <ga.h>
#include <ga/make_ga.h>
#include <utils/eoFuncPtrStat.h>

#include <eval/oneMaxEval.h>

#include "t-mpi-common.h"

#include "make_ls_op.h"
#include "LocalSearch.h"
#include "make_checkpoint_ga.h"

using namespace std;
using namespace eo::mpi;
using namespace MPI;
using namespace boost::mpi;
using namespace boost;

template < typename Atom >
class SquareMatrix : public std::vector< Atom >, public eoObject, public eoPersistent
{
public:
    // Ctor : sets size
    SquareMatrix(unsigned _s = 0) : std::vector< Atom >(_s*_s), rSize(_s) {}

    /** simple accessor */
    const Atom& operator()(unsigned _i, unsigned _j) const
    {
	return this->operator[](_i*rSize + _j);
    }

    /** reference - to set values */
    Atom& operator()(unsigned _i, unsigned _j)
    {
	return this->operator[](_i*rSize + _j);
    }

    /** returns a vector value in the row _i */
    std::vector< Atom > operator()(unsigned _i)
    {
	typename std::vector< Atom >::iterator begin = this->begin();
	std::advance(begin, _i * rSize);
	typename std::vector< Atom >::iterator end = begin;
	std::advance(end, rSize);
	std::vector< Atom > vec(begin, end);
	return vec;
    }

    /** just in case */
    virtual void printOn(std::ostream & _os) const
    {
	unsigned index=0;
	for (unsigned i=0; i<rSize; i++)
	    {
		for (unsigned j=0; j<rSize; j++)
		    {
			_os << this->operator[](index++) << " " ;
		    }
		_os << std::endl;
	    }
	_os << std::endl;
    }

    virtual void readFrom(std::istream& /*_is*/)
    {
    }

    virtual std::string className() const {return "SquareMatrix";};

    size_t size() const {return rSize;}

private:
    unsigned rSize;              // row size (== number of columns!)
};

class MigrationMatrix : public SquareMatrix< double >
{
public:
    MigrationMatrix(unsigned s = 0) : SquareMatrix< double >(s) {}

    virtual void printOn(std::ostream & os) const
    {
	os << "Migration probability (in %) among islands" << std::endl;
	os << "\t" << "1/l";

	for (size_t i = 0; i < this->size() - 1; ++i)
	    {
		os << "\t" << i * 2 + 1 << "f";
	    }

	os << std::endl;

	for (size_t i = 0; i < this->size(); ++i)
	    {
		if (i == 0)
		    os << "1/l";
		else
		    os << i*2-1 << "f";

		for (size_t j = 0; j < this->size(); ++j)
		    {
			os << "\t" << floor((*this)(i,j) * 100) / 1000;
		    }

		os << std::endl;
	    }

	os << std::endl;
	os << "sum";

	double sum;

	for (size_t i = 0; i < this->size(); ++i)
	    {
		sum = 0;
		for (size_t j = 0; j < this->size(); ++j)
		    {
			sum = sum + (*this)(j,i);
		    }
		os << "\t" << floor(sum * 100) / 1000;
	    }
	os << std::endl;
    }
};

class InitMatrix : public eoUF< SquareMatrix<double>&, void >
{
public:
    InitMatrix(bool initG = false, double same = 90) : _initG(initG), _same(same) {}

    void operator()(SquareMatrix<double>& matrix)
    {
	double sum;

	for (size_t i = 0; i < matrix.size(); ++i)
	    {
		sum = 0;

		for (size_t j = 0; j < matrix.size(); ++j)
		    {
			if (i == j)
			    matrix(i,j) = _same * 10;
			else
			    {
				if (_initG)
				    matrix(i,j) = (1000 - _same * 10) / (matrix.size()-1);
				else
				    matrix(i,j) = rand();

				sum = sum + matrix(i,j);
			    }
		    }

		for (size_t j = 0; j < matrix.size(); ++j)
		    {
			if (i != j)
			    {
				if (sum == 0)
				    matrix(i,j) = 0;
				else
				    matrix(i,j) = matrix(i,j) / sum * (1000 - _same * 10);
			    }
		    }
	    }
    }

private:
    bool _initG;
    double _same;
};

template <typename EOT>
double binary_value(const EOT& _chrom)
{
    double sum = 0.0;
    for (unsigned i=0; i<_chrom.size(); i++)
	sum += _chrom[i];
    return sum;
}

typedef eoBit<double> EOT;

std::queue< std::pair< size_t, eoPop<EOT> > > imm;
std::queue< std::pair< size_t, eoPop<EOT> > > em;

void print_sum( eoPop<EOT>& pop )
{
    communicator world;
    int size = pop.size();
    eo::log << eo::progress << "F" << pop.size() << " "; eo::log.flush();
    int sum = -1;
    all_reduce( world, size, sum, std::plus<int>() );
    if ( 0 == world.rank() )
    	{
	    eo::log << eo::progress << "sum: " << sum << endl; eo::log.flush();
	}
}

template < typename T >
class FixedValue : public eoUpdater, public eoValueParam<T>
{
public:
    FixedValue( T value = 0, std::string label = "Value" ) : eoValueParam<T>(value, label) {}

    virtual void operator()() { /* nothing to do */ }
};

template < typename EOT > size_t getPopSize(const eoPop< EOT >& pop) { return pop.size(); }

class GetInputOutput : public eoUpdater, public eoValueParam<int>
{
public:
    GetInputOutput( int value, std::string label = "Value" ) : eoValueParam<int>(value, label) {}

    virtual void operator()()
    {
	value() = 0;
    }
};

class GetMigrationProbability : public eoUpdater, public eoValueParam<double>
{
public:
    GetMigrationProbability( const vector< double >& vecProba, const size_t isl, std::string label = "Value" )
	: eoValueParam<double>(0, label), _vecProba(vecProba), _isl(isl) {}

    virtual void operator()()
    {
	value() = _vecProba[_isl] / 10;
    }

private:
    const vector< double >& _vecProba;
    const size_t _isl;
};

class GetSumVectorProbability : public eoUpdater, public eoValueParam<double>
{
public:
    GetSumVectorProbability( const vector< double >& vecProba, std::string label = "Value" )
	: eoValueParam<double>(0, label), _vecProba(vecProba) {}

    virtual void operator()()
    {
	value() = accumulate(_vecProba.begin(), _vecProba.end(), 0) / 10.;
    }

private:
    const vector< double >& _vecProba;
};

template <class EOT>
class AverageDeltaFitnessStat : public eoStat<EOT, typename EOT::Fitness>
{
public :
    using eoStat<EOT, typename EOT::Fitness>::value;

    typedef typename EOT::Fitness Fitness;

    AverageDeltaFitnessStat(std::string _description = "Average Delta Fitness")
	: eoStat<EOT, Fitness>(Fitness(), _description) {}

    static Fitness sumFitness(double _sum, const EOT& _eot){
        _sum += _eot.fitness() / _eot.getLastFitness();
        return _sum;
    }

    AverageDeltaFitnessStat(double _value, std::string _desc) : eoStat<EOT, double>(_value, _desc) {}

    virtual void operator()(const eoPop<EOT>& _pop){
	doit(_pop, Fitness()); // specializations for scalar and std::vector
    }

    virtual std::string className(void) const { return "AverageDeltaFitnessStat"; }

private :
    // Default behavior
    template <class T>
    void doit(const eoPop<EOT>& _pop, T)
    {
        Fitness v = std::accumulate(_pop.begin(), _pop.end(), Fitness(0.0), AverageDeltaFitnessStat::sumFitness);
        value() = v / _pop.size();
    }

};

/****************************************
 * Define MPI tags used for the purpose *
 ****************************************/

enum { INDIVIDUALS=0, FEEDBACKS=1, END=2 };

int main(int argc, char *argv[])
{
    /******************************
     * Initialisation de MPI + EO *
     ******************************/

    environment env(argc, argv);
    communicator world;

    eoParser parser(argc, argv);
    eoState state;    // keeps all things allocated

    /*****************************
     * Definition des paramètres *
     *****************************/

    // a
    double alpha = parser.createParam(double(0.8), "alpha", "Alpha", 'a', "Islands Model").value();
    // b
    double beta = parser.createParam(double(0.99), "beta", "Beta", 'b', "Islands Model").value();
    // p
    size_t probaMin = parser.createParam(size_t(10), "probaMin", "Minimum probability to stay in the same island", 'p', "Islands Model").value();
    // d
    size_t probaSame = parser.createParam(size_t(100), "probaSame", "Probability for an individual to stay in the same island", 'd', "Islands Model").value();
    // r
    size_t reward = parser.createParam(size_t(2), "reward", "reward", 'r', "Islands Model").value();
    size_t penalty = parser.createParam(size_t(1), "penalty", "penalty", 0, "Islands Model").value();
    // I
    bool initG = parser.createParam(bool(true), "initG", "initG", 'I', "Islands Model").value();

    bool printBest = parser.createParam(false, "printBestStat", "Print Best/avg/stdev every gen.", '\0', "Output").value();

    /****************************
     * Il faut au moins 4 nœuds *
     ****************************/

    const int ALL = world.size();
    const int RANK = world.rank();

    if ( ALL < 4 )
    	{
    	    if ( 0 == RANK )
    		{
    		    cerr << "Needs at least 4 processes to be launched!" << endl;
    		}
    	    return 0;
    	}

    /*********************************
     * Déclaration des composants EO *
     *********************************/

    unsigned chromSize = parser.getORcreateParam(unsigned(1500), "chromSize", "The length of the bitstrings", 'n',"Problem").value();
    eoInit<EOT>& init = make_genotype(parser, state, EOT(), 0);

    oneMaxEval<EOT> mainEval;
    eoEvalFuncCounter<EOT> eval(mainEval);

    unsigned popSize = parser.getORcreateParam(unsigned(100), "popSize", "Population Size", 'P', "Evolution Engine").value();
    eoPop<EOT>& pop = make_pop(parser, state, init);

    parser.getORcreateParam(unsigned(10000), "maxGen", "Maximum number of generations () = none)",'G',"Stopping criterion");
    double targetFitness = parser.getORcreateParam(double(chromSize), "targetFitness", "Stop when fitness reaches",'T', "Stopping criterion").value();
    eoContinue<EOT>& term = make_continue(parser, state, eval);

    eoCheckPoint<EOT>& checkpoint = state.storeFunctor( new eoCheckPoint<EOT>( term ) );

    eoGenContinue<EOT> algo_cont( 1 ); // just for algo

    const int GALL = world.size();
    const int GRANK = world.rank();

    /****************************************
     * Distribution des opérateurs aux iles *
     ****************************************/

    unsigned nbMove = parser.getORcreateParam(unsigned(1), "nbMove", "Number of move allowed for local search",'N',"Local search").value();
    my::LocalSearch<EOT> seqOp(nbMove);

    eoMonOp<EOT>* ptMon = NULL;
    if ( GRANK == 0 )
	{
	    eo::log << eo::logging << GRANK << ": bitflip ";
	    ptMon = new eoBitMutation<EOT>( 1, true );
	}
    else
	{
	    eo::log << eo::logging << GRANK << ": kflip(" << (GRANK-1) * 2 + 1 << ") ";
	    ptMon = new eoDetBitFlip<EOT>( (GRANK-1) * 2 + 1 );
	}
    eo::log << eo::logging << endl;
    eo::log.flush();
    state.storeFunctor(ptMon);
    seqOp.add(*ptMon, 1);

    /********************************
     * Initialize generic algorithm *
     ********************************/

    string comment = "Selection: DetTour(T), StochTour(t), Roulette, Ranking(p,e) or Sequential(ordered/unordered)";
    parser.getORcreateParam(eoParamParamType("DetTour(20)"), "selection", comment, 'S', "Evolution Engine");
    parser.getORcreateParam(eoParamParamType("Plus"), "replacement", "Replacement: Comma, Plus or EPTour(T), SSGAWorst, SSGADet(T), SSGAStoch(t)", 'R', "Evolution Engine");

    eoAlgo<EOT>& ga = make_algo_scalar(parser, state, eval, algo_cont, seqOp);

    /**************
     * EO routine *
     **************/

    make_parallel(parser);
    make_verbose(parser);
    make_help(parser);

    /******************************************************************************
     * Création de la matrice de transition et distribution aux iles des vecteurs *
     ******************************************************************************/

    MigrationMatrix probabilities( GALL );
    InitMatrix initmatrix( initG, probaSame );
    vector< double > vecProba( GALL );

    if ( 0 == GRANK )
    	{
	    initmatrix( probabilities );
    	    cout << probabilities;
	    vecProba = probabilities(GRANK);

	    for (int i = 1; i < GALL; ++i)
		{
		    world.send( i, 0, probabilities(i) );
		}
    	}
    else
	{
	    world.recv( 0, 0, vecProba );
	}

    /***************************************
     * Déclaration des opérateurs de stats *
     ***************************************/

    std::ostringstream ss;
    ss << "monitor.csv." << GRANK;
    eoFileMonitor& fileMonitor = state.storeFunctor( new eoFileMonitor( ss.str(), " ", false, true ) );
    checkpoint.add(fileMonitor);

    eoStdoutMonitor* stdMonitor = NULL;
    if (printBest)
	{
	    stdMonitor = new eoStdoutMonitor("\t", 8);
	    state.storeFunctor( stdMonitor );
	    checkpoint.add(*stdMonitor);
	}

    FixedValue<unsigned int>& islNum = state.storeFunctor( new FixedValue<unsigned int>( GRANK, "Island" ) );
    fileMonitor.add(islNum);
    if (printBest) { stdMonitor->add(islNum); }

    eoTimeCounter& tCounter = state.storeFunctor( new eoTimeCounter );
    checkpoint.add(tCounter);
    fileMonitor.add(tCounter);
    if (printBest) { stdMonitor->add(tCounter); }

    eoGenCounter& genCounter = state.storeFunctor( new eoGenCounter( 0, "migration" ) );
    checkpoint.add(genCounter);
    fileMonitor.add(genCounter);
    if (printBest) { stdMonitor->add(genCounter); }

    std::ostringstream ss_size;
    ss_size << "nb_individual_isl" << GRANK;
    eoFuncPtrStat<EOT, size_t>& popSizeStat = makeFuncPtrStat( getPopSize<EOT>, state, ss_size.str() );
    checkpoint.add(popSizeStat);
    fileMonitor.add(popSizeStat);
    if (printBest) { stdMonitor->add(popSizeStat); }

    std::ostringstream ss_avg;
    ss_avg << "avg_ones_isl" << GRANK;
    eoAverageStat<EOT>& avg = state.storeFunctor( new eoAverageStat<EOT>( ss_avg.str() ) );
    checkpoint.add(avg);
    fileMonitor.add(avg);
    if (printBest) { stdMonitor->add(avg); }

    std::ostringstream ss_delta;
    ss_delta << "delta_avg_ones_isl" << GRANK;
    AverageDeltaFitnessStat<EOT>& avg_delta = state.storeFunctor( new AverageDeltaFitnessStat<EOT>( ss_delta.str() ) );
    checkpoint.add(avg_delta);
    fileMonitor.add(avg_delta);
    if (printBest) { stdMonitor->add(avg_delta); }

    std::ostringstream ss_best;
    ss_best << "best_value_isl" << GRANK;
    eoBestFitnessStat<EOT>& best = state.storeFunctor( new eoBestFitnessStat<EOT>( ss_best.str() ) );
    // checkpoint.add(best);
    fileMonitor.add(best);
    if (printBest) { stdMonitor->add(best); }

    std::ostringstream ss_input;
    ss_input << "nb_input_ind_isl" << GRANK;
    GetInputOutput& input = state.storeFunctor( new GetInputOutput( 0, ss_input.str() ) );
    fileMonitor.add(input);
    if (printBest) { stdMonitor->add(input); }

    std::ostringstream ss_output;
    ss_output << "nb_output_ind_isl" << GRANK;
    GetInputOutput& output = state.storeFunctor( new GetInputOutput( 0, ss_output.str() ) );
    fileMonitor.add(output);
    if (printBest) { stdMonitor->add(output); }

    for (int i = 0; i < vecProba.size(); ++i)
	{
	    std::ostringstream ss;
	    ss << "P" << GRANK << "to" << i;
	    GetMigrationProbability& migProba = state.storeFunctor( new GetMigrationProbability( vecProba, i, ss.str() ) );
	    checkpoint.add(migProba);
	    fileMonitor.add(migProba);
	    if (printBest) { stdMonitor->add(migProba); }
	}

    std::ostringstream ss_sum;
    ss_sum << "P" << GRANK << "to*";
    GetSumVectorProbability& sumProba = state.storeFunctor( new GetSumVectorProbability( vecProba, ss_sum.str() ) );
    checkpoint.add(sumProba);
    fileMonitor.add(sumProba);
    if (printBest) { stdMonitor->add(sumProba); }

    vector<GetInputOutput*> outputsPerIsl;
    for (size_t i = 0; i < GALL; ++i)
	{
	    std::ostringstream ss;
	    ss << "nb_migrants_isl" << GRANK << "to" << i;
	    GetInputOutput* out = new GetInputOutput( 0, ss.str() );
	    outputsPerIsl.push_back( out );
	    state.storeFunctor( out );
	    fileMonitor.add(*out);
	    if (printBest) { stdMonitor->add(*out); }
	}

    /******************************************
     * Get the population size of all islands *
     ******************************************/

    world.barrier();
    print_sum(pop);

    /***************
     * Rock & Roll *
     ***************/

    apply<EOT>(eval, pop);

    for (int i = 0; i < pop.size(); ++i)
	{
	    pop[i].addIsland(GRANK);
	}

    vector< typename EOT::Fitness > vecAvg(GALL, 0);
    vector< typename EOT::Fitness > vecFeedbacks(GALL, 0);

    bool check = true;
    bool reached = false;

    do
	{

	    best(pop);
	    all_reduce(world, best.value(), best.value(), maximum<double>());

	    /**************************************************
	     * Initialize a few variables for each generation *
	     **************************************************/

	    eval.value(0);
	    input.value(0);
	    output.value(0);

	    vector< request > reqs;

	    /****************************
	     * Update transition vector *
	     ****************************/
	    {
		int best = -1;
		typename EOT::Fitness max = 0;

		for (size_t i = 0; i < GALL; ++i)
		    {
			if (vecFeedbacks[i] > max)
			    {
				best = i;
				max = vecFeedbacks[i];
			    }
		    }

		//computation of epsilon vector (norm is 1)
	    	double sum = 0;

		vector< double > epsilon( GALL );

		for ( size_t k = 0; k < GALL; ++k )
		    {
			epsilon[k] = rng.rand() % 100;
			sum += epsilon[k];
		    }

		for ( size_t k = 0; k < GALL; ++k )
		    {
			epsilon[k] = sum ? epsilon[k] / sum : 0;
		    }

		/*******************************************************************************
		 * Si p_i^t est le vecteur de migration de ile numéro i au temps t             *
		 * alors on peut faire un update "baysien" comme ça:                           *
		 *                                                                             *
		 * p_i^{t+1} =  b  *  ( a * p_i^t + (1 - a) * select )  +  (1 - b) * epsilon   *
		 *                                                                             *
		 * où:                                                                         *
		 * - a et b sont les coefficients qui reglent respectivement l'exploitation et *
		 * l'exploration                                                               *
		 * - select est le vecteur avec que des 0 sauf pour l'ile qui semble la plus   *
		 * prometteuse où cela vaut 1.                                                 *
		 * - epsilon un vecteur de 'bruit' (de norme 1) avec des coefficient           *
		 * uniforme entre 0.0 et 1.0                                                   *
		 *******************************************************************************/

		if (best < 0) //aucune ile n'améliore donc on rééquilibre
		    {
			for ( size_t i = 0; i < GALL; ++i )
			    {
				vecProba[i] = beta * vecProba[i] + (1 - beta) * 1000 * epsilon[i];
			    }
		    }
		else
		    {
			for (size_t i = 0; i < GALL; ++i)
			    {
				if ( i == best )
				    {
					vecProba[i] = beta * ( alpha * vecProba[i] + (1 - alpha) * 1000 ) + (1 - beta) * 1000 * epsilon[i];
				    }
				else
				    {
					vecProba[i] = beta * ( alpha * vecProba[i] ) + (1 - beta) * 1000 * epsilon[i];
				    }
			    }
		    }
	    }

	    /***********************************
	     * Send individuals to all islands *
	     ***********************************/
	    {
		vector< eoPop<EOT> > pops( GALL );

		/*************
		 * Selection *
		 *************/

		for (size_t i = 0; i < pop.size(); ++i)
		    {
			double s = 0;
			int r = rng.rand() % 1000 + 1;

			size_t j;
			for ( j = 0; j < GALL && r > s; ++j )
			    {
				s += vecProba[j];
			    }
			--j;

			pops[j].push_back(pop[i]);
		    }

		for (size_t i = 0; i < GALL; ++i)
		    {
			if (i == GRANK) { continue; }
			outputsPerIsl[i]->value( pops[i].size() );
			output.value( output.value() + pops[i].size() );
		    }

		pop.clear();

		if ( best.value() < targetFitness )
		    {
			for ( size_t i = 0; i < GALL; ++i )
			    {
				reqs.push_back( world.isend( i, INDIVIDUALS, pops[i] ) );
			    }
		    }
	    }

	    vector< eoPop<EOT> > pops( GALL );

	    /****************************************
	     * Receive individuals from all islands *
	     ****************************************/
	    {
		for (size_t i = 0; i < GALL; ++i)
		    {
			if ( best.value() >= targetFitness )
			    {
				if (i == GRANK) { continue; }
			    }
			reqs.push_back( world.irecv( i, INDIVIDUALS, pops[i] ) );
		    }
	    }

	    /****************************
	     * Process all MPI requests *
	     ****************************/

	    wait_all( reqs.begin(), reqs.end() );
	    reqs.clear();

	    /*********************
	     * Update population *
	     *********************/
	    {
		for (size_t i = 0; i < GALL; ++i)
		    {
			eoPop<EOT>& newpop = pops[i];
			for (size_t j = 0; j < newpop.size(); ++j)
			    {
				pop.push_back( newpop[j] );
			    }
			if (i != GRANK)
			    {
				input.value( input.value() + newpop.size() );
			    }
		    }
	    }

	    /******************************************************************
	     * Memorize last fitness and island of population before evolving *
	     ******************************************************************/

	    for (int i = 0; i < pop.size(); ++i)
		{
		    pop[i].addFitness();
		    pop[i].addIsland(GRANK);
		}

	    /**********
	     * Evolve *
	     **********/

	    run_ea(ga, pop);
	    algo_cont.reset();

	    /************************************************
	     * Send feedbacks back to all islands (ANALYSE) *
	     ************************************************/

	    for (size_t i = 0; i < GALL; ++i)
		{
		    vector<typename EOT::Fitness> sums(GALL, 0);
		    vector<int> nbs(GALL, 0);
		    for (int i = 0; i < pop.size(); ++i)
			{
			    sums[pop[i].getLastIsland()] += pop[i].fitness() - pop[i].getLastFitness();
			    ++nbs[pop[i].getLastIsland()];
			}
		    if ( best.value() < targetFitness )
			{
			    reqs.push_back( world.isend( i, FEEDBACKS, sums[i] / nbs[i] ) );
			}
		}

	    /**************************************
	     * Receive feedbacks from all islands *
	     **************************************/

	    for (size_t i = 0; i < GALL; ++i)
		{
		    if ( best.value() >= targetFitness )
			{
			    if (i == GRANK) { continue; }
			}
		    reqs.push_back( world.irecv( i, FEEDBACKS, vecFeedbacks[i] ) );
		}

	    /****************************
	     * Process all MPI requests *
	     ****************************/

	    wait_all( reqs.begin(), reqs.end() );
	    reqs.clear();

	    /*******************************
	     * Stopping criteria reached ? *
	     *******************************/

	    // if ( best.value() >= targetFitness )
	    // 	{
	    // 	    world.abort(0);
	    // 	}

	    // reached = ( best.value() >= targetFitness );

	    // for (size_t i = 0; i < GALL; ++i)
	    // 	{
	    // 	    reqs.push_back( world.isend( i, TEST, reached ) );
	    // 	}

	    // for (size_t i = 0; i < GALL; ++i)
	    // 	{
	    // 	    bool isReached = false;
	    // 	    reqs.push_back( world.irecv( i, TEST, isReached ) );
	    // 	    reached |= isReached;
	    // 	}

	    // wait_all( reqs.begin(), reqs.end() );
	    // reqs.clear();

	    // cout << "r" << reached << " "; cout.flush();
	}
    while ( checkpoint(pop) );

    /*************************************************************************
     * MAJ de la matrice de transition et récupération des vecteurs des iles *
     *************************************************************************/

    world.barrier();
    if ( GRANK > 0 )
    	{
	    world.send( 0, 0, vecProba );
	}
    else
	{
	    for (int i = 1; i < GALL; ++i)
		{
		    vector<double> proba(GALL);
		    world.recv( i, 0, proba );
		    for (int j = 0; j < proba.size(); ++j)
			{
			    probabilities(i,j) = proba[j];
			}
		}
	    for (int j = 0; j < vecProba.size(); ++j)
		{
		    probabilities(0,j) = vecProba[j];
		}

	    cout << probabilities;
	    cout.flush();
	}

    /******************************************
     * Get the population size of all islands *
     ******************************************/

    world.barrier();
    print_sum(pop);

    /*********
     * DEBUG *
     *********/

    world.barrier();
    eo::log << eo::progress;
    copy(vecAvg.begin(), vecAvg.end(), ostream_iterator<typename EOT::Fitness>(eo::log, " "));
    eo::log << endl; eo::log.flush();

    return 0;
}

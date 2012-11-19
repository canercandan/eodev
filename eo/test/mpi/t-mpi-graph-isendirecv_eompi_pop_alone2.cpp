#include <stdexcept>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <mpi.h>

#include <thread>
#include <mutex>
#include <atomic>
#include <queue>

#include <eompi.h>
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
	// std::copy(vec.begin(), vec.end(), std::ostream_iterator<Atom>(std::cout, " "));
	// std::cout << std::endl;
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
    InitMatrix(int initG = 0, double same = 90) : _initG(initG), _same(same) {}

    void operator()(SquareMatrix<double>& matrix)
    {
	double sum;

	for (size_t i = 0; i < matrix.size(); ++i)
	    {
		sum = 0;

		for (size_t j = 0; j < matrix.size(); ++j)
		    {
			if (i == j)
			    {
				matrix(i,j) = _same * 10;
			    }
			else
			    {
				if (_initG == 0)
				    matrix(i,j) = rand();
				else
				    matrix(i,j) = (1000 - _same * 10) / (matrix.size() - 1);
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
    int _initG;
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

int get_sum( eoPop<EOT>& pop, mpi::communicator& topology )
{
    int size = pop.size();
    int sum = -1;
    for (int i = 0; i < topology.size(); ++i)
	{
	    MPI_Reduce( &size, &sum, 1, MPI_INT, MPI_SUM, i, topology.mpi_comm() );
	}
    return sum;
}

void print_sum( eoPop<EOT>& pop, mpi::communicator& topology )
{
    cout << "F" << pop.size() << " "; cout.flush();
    int sum = get_sum(pop, topology);
    if ( 0 == topology.rank() )
    	{
	    cout << "sum: " << sum << endl; cout.flush();
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

// template < typename EOT >
// class GetOutput : public eoUF< const eoPop< EOT >&, int >
// {
// public:
//     GetOutput( int size ) : _size(size) {}

//     int operator()( const eoPop< EOT >& pop )
//     {
// 	int dist = _size - pop.size();
// 	_size = pop.size();
// 	return dist > 0 ? dist : 0;
//     }

// private:
//     int _size;
// };

// template < typename EOT >
// class GetInput : public eoUF< const eoPop< EOT >&, int >
// {
// public:
//     GetInput( int size ) : _size(size) {}

//     int operator()( const eoPop< EOT >& pop )
//     {
// 	int dist = pop.size() - _size;
// 	_size = pop.size();
// 	return dist > 0 ? dist : 0;
//     }

// private:
//     int _size;
// };


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

int main(int argc, char *argv[])
{
    /******************************
     * Initialisation de MPI + EO *
     ******************************/

    // Node::init( ac, av, MPI_THREAD_MULTIPLE );
    Node::init( argc, argv/*, MPI_THREAD_SERIALIZED*/ );
    eoParser parser(argc, argv);
    eoState state;    // keeps all things allocated

    /*****************************
     * Definition des paramètres *
     *****************************/

    // a
    double alpha = parser.createParam(double(0.9), "alpha", "Alpha", 'a', "Islands Model").value();
    // b
    double beta = parser.createParam(double(0.9), "beta", "Beta", 'b', "Islands Model").value();
    // p
    size_t probaMin = parser.createParam(size_t(10), "probaMin", "Minimum probability to stay in the same island", 'p', "Islands Model").value();
    // d
    size_t probaSame = parser.createParam(size_t(90), "probaSame", "Probability for an individual to stay in the same island", 'd', "Islands Model").value();
    // r
    size_t reward = parser.createParam(size_t(2), "reward", "reward", 'r', "Islands Model").value();
    size_t penalty = parser.createParam(size_t(1), "penalty", "penalty", 0, "Islands Model").value();

    /****************************
     * Il faut au moins 4 nœuds *
     ****************************/

    const int ALL = Node::comm().size();
    const int RANK = Node::comm().rank();

    // mpi_debug();

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

    // eoUniformGenerator<typename EOT::AtomType> gen(0, 2);
    // eoInitFixedLength<EOT> init(dimSize, gen);

    eoInit<EOT>& init = make_genotype(parser, state, EOT(), 0);

    // eoEvalFuncPtr<EOT> eval(binary_value);

    oneMaxEval<EOT> mainEval;
    eoEvalFuncCounter<EOT> eval(mainEval);

    // eoGenOp<EOT>& op = my::make_ls_op(parser, state, init);

    eoPop<EOT>& pop = make_pop(parser, state, init);

    // eoValueParam<unsigned long>& maxEvalParam = parser.getORcreateParam((unsigned long)0, "maxEval", "Maximum number of evaluations (0 = none)", 'E', "Stopping criterion");
    eoValueParam<unsigned>& maxGenParam = parser.getORcreateParam(unsigned(50), "maxGen", "Maximum number of generations () = none)",'G',"Stopping criterion");
    eoContinue<EOT>& term = make_continue(parser, state, eval);

    // eoCheckPoint<EOT>& checkpoint = my::make_checkpoint(parser, state, eval, term);

    eoCheckPoint<EOT>& checkpoint = state.storeFunctor( new eoCheckPoint<EOT>( term ) );

    eoGenContinue<EOT> cont( 1 ); // just for algo

    // mpi::ring topology;
    mpi::complete topology;

    const size_t GALL = topology.size();
    const size_t GRANK = topology.rank();

    if ( 0 == GRANK )
	{
	    // topology.print();
	    topology.test();
	}

    /****************************************
     * Distribution des opérateurs aux iles *
     ****************************************/

    my::LocalSearch<EOT> seqOp(100);

    eoMonOp<EOT>* ptMon = NULL;
    if ( GRANK == GALL - 1 )
	{
	    cout << "bitflip ";
	    ptMon = new eoBitMutation<EOT>( 1, true );
	}
    else
	{
	    cout << "kflip(" << GRANK * 2 + 1 << ") ";
	    ptMon = new eoDetBitFlip<EOT>( GRANK * 2 + 1 );
	}
    cout << endl;
    cout.flush();
    state.storeFunctor(ptMon);
    seqOp.add(*ptMon, 1);

    /***********************************
     * Déclaration de l'algo générique *
     ***********************************/

    string comment = "Selection: DetTour(T), StochTour(t), Roulette, Ranking(p,e) or Sequential(ordered/unordered)";
    eoValueParam<eoParamParamType>& selectionParam = parser.getORcreateParam(eoParamParamType("DetTour(20)"), "selection", comment, 'S', "Evolution Engine");
    eoValueParam<eoParamParamType>& replacementParam = parser.getORcreateParam(eoParamParamType("Plus"), "replacement", "Replacement: Comma, Plus or EPTour(T), SSGAWorst, SSGADet(T), SSGAStoch(t)", 'R', "Evolution Engine");

    eoAlgo<EOT>& ga = make_algo_scalar(parser, state, eval, cont, seqOp);

    /**************
     * Routine EO *
     **************/

    make_parallel(parser);
    make_verbose(parser);
    make_help(parser);

    topology.barrier();

    apply<EOT>(eval, pop);

    Intracomm comm = COMM_WORLD;

    size_t size = topology.neighbors_count();

    vector<int> to = topology.to();
    vector<int> from = topology.from();

    cout << to.size() << " " << from.size() << endl;

    topology.barrier();

    /******************************************************************************
     * Création de la matrice de transition et distribution aux iles des vecteurs *
     ******************************************************************************/

    MigrationMatrix probabilities( GALL );
    InitMatrix initmatrix;
    vector< double > vecProba( GALL );

    if ( 0 == GRANK )
    	{
	    initmatrix( probabilities );
    	    cout << probabilities;
	    vecProba = probabilities(GRANK);

	    for (int i = 1; i < GALL; ++i)
		{
		    MPI_Send( probabilities(i).data(), GALL, MPI_DOUBLE, i, 0, topology.mpi_comm() );
		}
    	}
    else
	{
	    MPI_Status stat;
	    MPI_Recv( vecProba.data(), GALL, MPI_DOUBLE, 0, 0, topology.mpi_comm(), &stat );
	}

    topology.barrier();

    print_sum(pop, topology);

    topology.barrier();

    size_t popSize = pop.size();

    vector< typename EOT::Fitness > vecEffectiveness(from.size(), 0);
    vector< int > vecMigrations(from.size(), 0);

    /***************************************
     * Déclaration des opérateurs de stats *
     ***************************************/

    std::ostringstream ss;
    ss << "monitor.csv." << GRANK;
    eoFileMonitor& monitor = state.storeFunctor( new eoFileMonitor( ss.str(), " ", false, true ) );
    checkpoint.add(monitor);

    FixedValue<unsigned int>& islNum = state.storeFunctor( new FixedValue<unsigned int>( GRANK, "Island" ) );
    monitor.add(islNum);

    eoGenCounter& genCounter = state.storeFunctor( new eoGenCounter( 0, "migration" ) );
    checkpoint.add(genCounter);
    monitor.add(genCounter);

    eoTimeCounter& tCounter = state.storeFunctor( new eoTimeCounter );
    checkpoint.add(tCounter);
    monitor.add(tCounter);

    eoFuncPtrStat<EOT, size_t>& popSizeStat = makeFuncPtrStat( getPopSize<EOT>, state, "nb_individual" );
    checkpoint.add(popSizeStat);
    monitor.add(popSizeStat);

    eoAverageStat<EOT>& avg = state.storeFunctor( new eoAverageStat<EOT>( "avg_ones" ) );
    checkpoint.add(avg);
    monitor.add(avg);

    // delta
    AverageDeltaFitnessStat<EOT>& avg_delta = state.storeFunctor( new AverageDeltaFitnessStat<EOT>( "delta_avg_ones" ) );
    checkpoint.add(avg_delta);
    monitor.add(avg_delta);

    eoBestFitnessStat<EOT>& best = state.storeFunctor( new eoBestFitnessStat<EOT>( "best_value" ) );
    checkpoint.add(best);
    monitor.add(best);

    // GetInput<EOT>& inputFunc = state.storeFunctor( new GetInput<EOT>( pop.size() ) );
    // eoFunctorStat<EOT,int>& inputStat = makeFunctorStat( inputFunc, state, "input" );
    // checkpoint.add(inputStat);
    // monitor.add(inputStat);

    // GetOutput<EOT>& outputFunc = state.storeFunctor( new GetOutput<EOT>( pop.size() ) );
    // eoFunctorStat<EOT, int>& outputStat = makeFunctorStat( outputFunc, state, "output" );
    // checkpoint.add(outputStat);
    // monitor.add(outputStat);

    GetInputOutput& input = state.storeFunctor( new GetInputOutput( 0, "input" ) );
    // checkpoint.add(input);
    monitor.add(input);

    GetInputOutput& output = state.storeFunctor( new GetInputOutput( 0, "output" ) );
    // checkpoint.add(output);
    monitor.add(output);

    for (int i = 0; i < vecProba.size(); ++i)
	{
	    std::ostringstream ss;
	    ss << "P" << GRANK << "to" << i;
	    GetMigrationProbability& migProba = state.storeFunctor( new GetMigrationProbability( vecProba, i, ss.str() ) );
	    checkpoint.add(migProba);
	    monitor.add(migProba);
	}

    {
	std::ostringstream ss;
	ss << "P" << GRANK << "to*";
	GetSumVectorProbability& sumProba = state.storeFunctor( new GetSumVectorProbability( vecProba, ss.str() ) );
	checkpoint.add(sumProba);
	monitor.add(sumProba);
    }

    vector<GetInputOutput*> outputsPerIsl;
    for (size_t i = 0; i < to.size(); ++i)
	{
	    std::ostringstream ss;
	    ss << "nb_migrants_isl" << GRANK << "to" << i;
	    GetInputOutput* out = new GetInputOutput( 0, ss.str() );
	    outputsPerIsl.push_back( out );
	    state.storeFunctor( out );
	    monitor.add(*out);
	}

    /***************
     * Rock & Roll *
     ***************/

    for (int i = 0; i < pop.size(); ++i)
	{
	    pop[i].addIsland(GRANK);
	}

    do
	{
	    eval.value(0);
	    input.value(0);
	    output.value(0);

	    for (int i = 0; i < pop.size(); ++i)
		{
		    pop[i].addFitness();
		}

	    // mpi::requests reqs;

	    // update transition
	    {
		int best = -1;
		typename EOT::Fitness max = 0;

		for (size_t i = 0; i < from.size(); ++i)
		    {
			typename EOT::Fitness mean = ( (typename EOT::Fitness)vecMigrations[i]/(typename EOT::Fitness)vecEffectiveness[i] );
			if (vecMigrations[i] && mean > max)
			    {
				best = i;
				max = mean;
			    }
		    }

		//computation of epsilon vector (norm is 1)
	    	double sum = 0;

		vector< double > epsilon( to.size() );

		for ( size_t k = 0; k < to.size(); ++k )
		    {
			epsilon[k] = rng.rand() % 100;
			sum += epsilon[k];
		    }

		for ( size_t k = 0; k < to.size(); ++k )
		    {
			epsilon[k] = sum ? epsilon[k] / sum : 0;
		    }

		// Si p_i^t est le vecteur de migration de ile numéro i au temps t
		// alors on peut faire un update "baysien" comme ça:
		//
		// p_i^{t+1} =  b  *  ( a * p_i^t + (1 - a) * select )  +  (1 - b) * epsilon
		//
		// où:
		// - a et b sont les coefficients qui reglent respectivement
		// l'exploitation et l'exploration
		// - select est le vecteur avec que des 0 sauf pour l'ile qui semble la
		// plus prometteuse où cela vaut 1.
		// - epsilon un vecteur de 'bruit' (de norme 1) avec des coefficient
		// uniforme entre 0.0 et 1.0

		if (best < 0) //aucune ile n'améliore donc on rééquilibre
		    {
			for ( size_t i = 0; i < to.size(); ++i )
			    {
				vecProba[i] = beta * vecProba[i] + (1 - beta) * 1000 * epsilon[i];
			    }
		    }
		else
		    {
			for (size_t i = 0; i < to.size(); ++i)
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

	    // send
	    {
		vector< eoPop<EOT> > pops( to.size(), eoPop<EOT>() );

		for (size_t i = 0; i < pop.size(); ++i)
		    {
			double s = 0;
			int r = rng.rand() % 1000 + 1;

			size_t j;
			for ( j = 0; j < to.size() && r > s; ++j )
			    {
				s += vecProba[j];
			    }
			--j;

			pops[j].push_back(pop[i]);
		    }

		for (size_t i = 0; i < to.size(); ++i)
		    {
			if (i == GRANK) { continue; }
			outputsPerIsl[i]->value( pops[i].size() );
			output.value( output.value() + pops[i].size() );
		    }

		pop.clear();

		for ( size_t i = 0; i < to.size(); ++i )
		    {
			eoserial::Object* obj = pops[i].pack();
			std::stringstream ss;
			obj->print( ss );
			delete obj;

			const std::string& str = ss.str();

			int size = str.size() + 1;
			int pos = 0;
			std::vector< char > buf( size, 0 );

			MPI_Pack( &size, 1, MPI_INT, buf.data(), buf.size(), &pos, topology.mpi_comm() );

			buf.resize(size + pos);

			MPI_Pack( (char*)str.c_str(), size, MPI_CHAR, buf.data(), buf.size(), &pos, topology.mpi_comm() );

			buf.resize(pos + 1000);
			MPI_Pack( &(vecEffectiveness[i]), 1, MPI_DOUBLE, buf.data(), buf.size(), &pos, topology.mpi_comm() );
			MPI_Pack( &(vecMigrations[i]), 1, MPI_INT, buf.data(), buf.size(), &pos, topology.mpi_comm() );
			buf.resize(pos);

			// MPI_Request req;
			// MPI_Isend( buf.data(), pos, MPI_PACKED, to[i], 0, topology.mpi_comm(), &req );
			// reqs.push_back(req);

			MPI_Send( buf.data(), pos, MPI_PACKED, to[i], 0, topology.mpi_comm() );
		    }
	    }

	    // // recv
	    // {
	    // 	// TODO: estimer la taille max
	    // 	size_t bufSize = sizeof(int) + (popSize*from.size()) + sizeof(typename EOT::Fitness) + sizeof(int);
	    // 	vector< vector< char > > from_buf( from.size(), vector< char >( bufSize * 200, 0 ) );

	    // 	for (size_t i = 0; i < from.size(); ++i)
	    // 	    {
	    // 		MPI_Request req;
	    // 		MPI_Irecv( from_buf[i].data(), from_buf[i].size(), MPI_PACKED, from[i], 0, topology.mpi_comm(), &req);
	    // 		reqs.push_back(req);
	    // 	    }

	    // 	vector<MPI_Status> status(reqs.size());
	    // 	MPI_Waitall(reqs.size(), reqs.data(), status.data());

	    // 	for (size_t i = 0; i < from.size(); ++i)
	    // 	    {
	    // 		int pos = 0;

	    // 		int size = -1;
	    // 		MPI_Unpack( from_buf[i].data(), from_buf[i].size(), &pos, &size, 1, MPI_INT, topology.mpi_comm() );
	    // 		string str(size, 0);
	    // 		MPI_Unpack( from_buf[i].data(), from_buf[i].size(), &pos, (char*)str.data(), size, MPI_CHAR, topology.mpi_comm() );

	    // 		eoserial::Object* obj = eoserial::Parser::parse( str );
	    // 		eoPop<EOT> newpop;
	    // 		newpop.unpack( obj );
	    // 		for (size_t j = 0; j < newpop.size(); ++j)
	    // 		    {
	    // 			pop.push_back( newpop[j] );
	    // 		    }
	    // 		delete obj;

	    // 		MPI_Unpack( from_buf[i].data(), from_buf[i].size(), &pos, &(vecEffectiveness[i]), 1, MPI_DOUBLE, topology.mpi_comm() );
	    // 		MPI_Unpack( from_buf[i].data(), from_buf[i].size(), &pos, &(vecMigrations[i]), 1, MPI_INT, topology.mpi_comm() );

	    // 		if (i != GRANK)
	    // 		    {
	    // 			input.value( input.value() + newpop.size() );
	    // 		    }
	    // 	    }
	    // }

	    // vector<MPI_Status> status(reqs.size());
	    // MPI_Waitall(reqs.size(), reqs.data(), status.data());

	    // recv
	    {
		for (size_t i = 0; i < from.size(); ++i)
		    {
			vector< char > buf( 10000, 0 );

			MPI_Status stat;
			MPI_Recv( buf.data(), buf.size(), MPI_PACKED, from[i], 0, topology.mpi_comm(), &stat);

			int pos = 0;

			int size = -1;
			MPI_Unpack( buf.data(), buf.size(), &pos, &size, 1, MPI_INT, topology.mpi_comm() );
			string str(size, 0);
			MPI_Unpack( buf.data(), buf.size(), &pos, (char*)str.data(), size, MPI_CHAR, topology.mpi_comm() );

			eoserial::Object* obj = eoserial::Parser::parse( str );
			eoPop<EOT> newpop;
			newpop.unpack( obj );
			for (size_t j = 0; j < newpop.size(); ++j)
			    {
				pop.push_back( newpop[j] );
			    }
			delete obj;

			MPI_Unpack( buf.data(), buf.size(), &pos, &(vecEffectiveness[i]), 1, MPI_DOUBLE, topology.mpi_comm() );
			MPI_Unpack( buf.data(), buf.size(), &pos, &(vecMigrations[i]), 1, MPI_INT, topology.mpi_comm() );

			if (i != GRANK)
			    {
				input.value( input.value() + newpop.size() );
			    }
		    }
	    }

	    // learn
	    {
	    }

	    if ( pop.size() )
		{
		    run_ea(ga, pop);
		}

	    cont.reset();

	    // analyse
	    {
		for (int i = 0; i < pop.size(); ++i)
		    {
			vecEffectiveness[pop[i].getLastIsland()] += pop[i].fitness() - pop[i].getLastFitness();
			++vecMigrations[pop[i].getLastIsland()];
		    }

		for (size_t i = 0; i < pop.size(); ++i)
		    {
			pop[i].addIsland(GRANK);
		    }
	    }

	    cout << "m "; cout.flush();

	}
    while ( checkpoint(pop) );

    topology.barrier();

    /*************************************************************************
     * MAJ de la matrice de transition et récupération des vecteurs des iles *
     *************************************************************************/

    if ( GRANK > 0 )
    	{
	    MPI_Send( vecProba.data(), GALL, MPI_DOUBLE, 0, 0, topology.mpi_comm() );
	}
    else
	{
	    for (int i = 1; i < GALL; ++i)
		{
		    vector<double> proba(GALL);
		    MPI_Status stat;
		    MPI_Recv( proba.data(), GALL, MPI_DOUBLE, i, 0, topology.mpi_comm(), &stat );
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

    topology.barrier();
    print_sum(pop, topology);

    topology.barrier();

    // if ( 0 == GRANK )
    // 	{
    // cout << "GRANK: " << GRANK << endl;
    for (int i = 0; i < pop.size(); ++i)
	{
	    // pop[i].printLastIslands();
	    cout << pop[i].getLastIslandsSize() << " ";
	    cout.flush();
	}
    cout << endl;
    cout.flush();
    // }

    // if ( 0 == GRANK )
    // 	{
    // 	    for (int i = 0; i < pop.size(); ++i)
    // 	    	{
    // 	    	    const vector<typename EOT::Fitness>& last = pop[i].getLastFitnesses();
    // 	    	    for (int j = 0; j < last.size(); ++j)
    // 	    		{
    // 	    		    cout << last[j] << " ";
    // 	    		}
    // 	    	    cout << endl;
    // 	    	}
    // 	}

    topology.barrier();

    copy(vecEffectiveness.begin(), vecEffectiveness.end(), ostream_iterator<typename EOT::Fitness>(cout, " "));
    cout << endl;

    // copy(vecProba.begin(), vecProba.end(), ostream_iterator<typename EOT::Fitness>(cout, " "));
    // cout << endl;

    return 0;
}

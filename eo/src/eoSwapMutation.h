// -*- mode: c++; c-indent-level: 4; c++-member-init-indent: 8; comment-column: 35; -*-

//-----------------------------------------------------------------------------
// eoSwapMutation.h
// (c) GeNeura Team, 2000 - EEAAX 2000 - Maarten Keijzer 2000
// (c) INRIA Futurs - Dolphin Team - Thomas Legrand 2007
/*
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

  Contact: todos@geneura.ugr.es, http://geneura.ugr.es
  thomas.legrand@lifl.fr
  Marc.Schoenauer@polytechnique.fr
  mak@dhi.dk
  caner.candan@univ-angers.fr
*/
//-----------------------------------------------------------------------------

#ifndef eoSwapMutation_h
#define eoSwapMutation_h

//-----------------------------------------------------------------------------


/**
 * Swap two components of a chromosome.
 *
 * @ingroup Variators
 */
template<class Chrom> class eoSwapMutation: public eoMonOp<Chrom>
{
public:

    /// CTor
    eoSwapMutation(const unsigned _howManySwaps=1): howManySwaps(_howManySwaps)
    {
        // consistency check
        if(howManySwaps < 1)
	    throw std::runtime_error("Invalid number of swaps in eoSwapMutation");
    }

    /// The class name.
    virtual std::string className() const { return "eoSwapMutation"; }

    /**
     * Swap two components of the given chromosome.
     * @param chrom The cromosome which is going to be changed.
     */
    bool operator()(Chrom& chrom)
    {
	unsigned i, j;

	for(unsigned int k = 0; k < howManySwaps; ++k)
	    {
		// generate two different indices
		i=eo::rng.random(chrom.size());
		do { j = eo::rng.random(chrom.size()); } while (i == j);

		// swap
		std::swap(chrom[i],chrom[j]);
	    }
	return true;
    }

private:
    unsigned int howManySwaps;
};

/**
 * Swap two components of a chromosome with the guarantee to have one improvement.
 *
 * @ingroup Variators
 */
template<typename EOT>
class eoFirstImprovementSwapMutation : public eoMonOp<EOT>
{
public:
    eoFirstImprovementSwapMutation(eoEvalFunc<EOT>& eval) : _eval(eval) {}

    /// The class name.
    virtual std::string className() const { return "eoFirstImprovementSwapMutation"; }

    /**
     * Swap two components of the given chromosome.
     * @param chrom The cromosome which is going to be changed.
     */
    bool operator()(EOT& sol)
    {
	for (size_t k = 0; k < sol.size()-1; ++k)
	    {
		EOT candidate = sol;
		candidate.invalidate();
		_op(candidate);
		_eval(candidate);
		if ( candidate.fitness() > sol.fitness() )
		    {
			sol = candidate;
			return true;
		    }
	    }
	return false;
    }

private:
    eoSwapMutation<EOT> _op;
    eoEvalFunc<EOT>& _eval;
};

/**
 * Swap two components of a chromosome while the best improvement wasnt been reached.
 *
 * @ingroup Variators
 */
template<typename EOT>
class eoRelativeBestImprovementSwapMutation : public eoMonOp<EOT>
{
public:
    /// ctor
    eoRelativeBestImprovementSwapMutation(eoEvalFunc<EOT>& eval) : _eval(eval) {}

    /// The class name.
    virtual std::string className() const { return "eoRelativeBestImprovementSwapMutation"; }

    /**
     * Swap two components of the given chromosome.
     * @param chrom The cromosome which is going to be changed.
     */
    bool operator()(EOT& sol)
    {
	// keep a best solution
	EOT best = sol;

	// select two indices from the initial solution
	size_t i, j;
	i = eo::rng.random(sol.size());
	do { j = eo::rng.random(sol.size()); } while (i == j);

	for (size_t k = 0; k < sol.size()-1; ++k)
	    {
		// create a candidate solution
		EOT candidate = sol;
		candidate.invalidate();

		// swap
		std::swap(candidate[i], candidate[j]);

		// evaluate
		_eval(candidate);

		// if the candidate is better than best solution we replace best by candidate
		if ( candidate.fitness() > best.fitness() )
		    {
			best = candidate;
		    }

		// increment j in order to swap with the other indexes
		do { j = (j+1) % candidate.size(); } while (i == j);
	    }

	// if the best solution is better than the initial one, we replace the solution by best
	if ( best.fitness() > sol.fitness() )
	    {
		sol = best;
		return true;
	    }

	return false;
    }

private:
    eoEvalFunc<EOT>& _eval;
};

/**
 * Swap two components of a chromosome while the best improvement wasnt been reached.
 *
 * @ingroup Variators
 */
template<typename EOT>
class eoBestImprovementSwapMutation : public eoMonOp<EOT>
{
public:
    /// ctor
    eoBestImprovementSwapMutation(eoEvalFunc<EOT>& eval) : _eval(eval) {}

    /// The class name.
    virtual std::string className() const { return "eoBestImprovementSwapMutation"; }

    /**
     * Swap two components of the given chromosome.
     * @param chrom The cromosome which is going to be changed.
     */
    bool operator()(EOT& sol)
    {
	// keep a best solution
	EOT best = sol;

	// select two indices from the initial solution
	size_t i, j;
	i = eo::rng.random(sol.size());

	for (size_t k = 0; k < sol.size()-1; ++k)
	    {
		do { j = eo::rng.random(sol.size()); } while (i == j);

		for (size_t l = 0; l < sol.size()-1; ++l)
		    {
			// create a candidate solution
			EOT candidate = sol;
			candidate.invalidate();

			// swap
			std::swap(candidate[i], candidate[j]);

			// evaluate
			_eval(candidate);

			// if the candidate is better than best solution we replace best by candidate
			if ( candidate.fitness() > best.fitness() )
			    {
				best = candidate;
			    }

			// increment j in order to swap with the other indexes
			do { j = (j+1) % sol.size(); } while (i == j);
		    }

		// increment i in order to swap with the other indexes
		i = (i+1) % sol.size();
	    }

	// if the best solution is better than the initial one, we replace the solution by best
	if ( best.fitness() > sol.fitness() )
	    {
		sol = best;
		return true;
	    }

	return false;
    }

private:
    eoEvalFunc<EOT>& _eval;
};

/** @example t-eoSwapMutation.cpp
 */

//-----------------------------------------------------------------------------
#endif

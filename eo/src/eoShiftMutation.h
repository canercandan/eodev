// -*- mode: c++; c-indent-level: 4; c++-member-init-indent: 8; comment-column: 35; -*-

//-----------------------------------------------------------------------------
// eoShiftMutation.h
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
 */
//-----------------------------------------------------------------------------

#ifndef eoShiftMutation_h
#define eoShiftMutation_h

//-----------------------------------------------------------------------------


/**
 * Shift two components of a chromosome.
 *
 * @ingroup Variators
 */
template<class EOT> class eoShiftMutation: public eoMonOp<EOT>
{
 public:

  typedef typename EOT::AtomType GeneType;

  /// CTor
  eoShiftMutation(){}


  /// The class name.
  virtual std::string className() const { return "eoShiftMutation"; }


  /**
   * Shift two components of the given eoosome.
   * @param _eo The cromosome which is going to be changed.
   */
  bool operator()(EOT& _eo)
    {

      unsigned i, j, from, to;
      GeneType tmp;

      // generate two different indices
      i=eo::rng.random(_eo.size());
      do j = eo::rng.random(_eo.size()); while (i == j);

      // indexes
      from=std::min(i,j);
      to=std::max(i,j);

      // keep the first component to change
      tmp=_eo[to];

      // shift
      for(unsigned int k=to ; k > from ; k--)
	  {
	      _eo[k]=_eo[k-1];
	  }

      // shift the first component
      _eo[from]=tmp;

      return true;
    }

};

template <typename EOT>
void shift(EOT& _eo, size_t i, size_t j)
{
    unsigned from, to;
    typename EOT::AtomType tmp;

    // indexes
    from=std::min(i,j);
    to=std::max(i,j);

    // keep the first component to change
    tmp=_eo[to];

    // shift
    for(unsigned int k=to ; k > from ; k--)
	{
	    _eo[k]=_eo[k-1];
	}

    // shift the first component
    _eo[from]=tmp;
}

/**
 * Shift two components of a chromosome with the guarantee to have one improvement.
 *
 * @ingroup Variators
 */
template<typename EOT>
class eoFirstImprovementShiftMutation : public eoMonOp<EOT>
{
public:
    eoFirstImprovementShiftMutation(eoEvalFunc<EOT>& eval) : _eval(eval) {}

    /// The class name.
    virtual std::string className() const { return "eoFirstImprovementShiftMutation"; }

    /**
     * Shift two components of the given chromosome.
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
    eoShiftMutation<EOT> _op;
    eoEvalFunc<EOT>& _eval;
};

/**
 * Shift two components of a chromosome while the best improvement wasnt been reached.
 *
 * @ingroup Variators
 */
template<typename EOT>
class eoRelativeBestImprovementShiftMutation : public eoMonOp<EOT>
{
public:
    /// ctor
    eoRelativeBestImprovementShiftMutation(eoEvalFunc<EOT>& eval) : _eval(eval) {}

    /// The class name.
    virtual std::string className() const { return "eoRelativeBestImprovementShiftMutation"; }

    /**
     * Shift two components of the given chromosome.
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

		// shift
		shift(candidate, i, j);

		// evaluate
		_eval(candidate);

		// if the candidate is better than best solution we replace best by candidate
		if ( candidate.fitness() > best.fitness() )
		    {
			best = candidate;
		    }

		// increment j in order to shift with the other indexes
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
 * Shift two components of a chromosome while the best improvement wasnt been reached.
 *
 * @ingroup Variators
 */
template<typename EOT>
class eoBestImprovementShiftMutation : public eoMonOp<EOT>
{
public:
    /// ctor
    eoBestImprovementShiftMutation(eoEvalFunc<EOT>& eval) : _eval(eval) {}

    /// The class name.
    virtual std::string className() const { return "eoBestImprovementShiftMutation"; }

    /**
     * Shift two components of the given chromosome.
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

			// shift
			shift(candidate, i, j);

			// evaluate
			_eval(candidate);

			// if the candidate is better than best solution we replace best by candidate
			if ( candidate.fitness() > best.fitness() )
			    {
				best = candidate;
			    }

			// increment j in order to shift with the other indexes
			do { j = (j+1) % sol.size(); } while (i == j);
		    }

		// increment i in order to shift with the other indexes
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

/** @example t-eoShiftMutation.cpp
 */


//-----------------------------------------------------------------------------
#endif

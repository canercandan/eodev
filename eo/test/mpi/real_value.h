#ifndef real_value_h
#define real_value_h

// EVAL
//-----------------------------------------------------------------------------
// a simple fitness function that computes the euclidian norm of a real vector
//    @param _indi A real-valued individual

template <typename Indi>
double real_value(const Indi & _indi)
{
    double sum = 0;
    for (unsigned i = 0; i < _indi.size(); i++)
        sum += _indi[i]*_indi[i];
    return (-sum);            // maximizing only
}

#endif

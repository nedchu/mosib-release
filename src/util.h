#ifndef util_h
#define util_h

#include "global.h"


double Jaccard(const VI & x, const VI & y);

template <class _InIt1, class _InIt2>
VI get_intersection(const _InIt1& x, const _InIt2& y)
{
	VI ans;
	std::set_intersection(x.begin(), x.end(),
		y.begin(), y.end(),
		back_inserter(ans));
	return ans;
}

VI get_intersection(const VI &x, VI::const_iterator y_b, VI::const_iterator y_e);
double get_duration(hclock::time_point s, hclock::time_point t);

#endif /* util_h */
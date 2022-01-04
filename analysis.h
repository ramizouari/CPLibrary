#ifndef __ANALYSIS_H__
#define __ANALYSIS_H__
#include "topology.h"

template<typename R,typename ...StructureMetaData>
R exp_n(R u, int n, StructureMetaData ...meta_info)
{
	R r(0,meta_info...),w(1, meta_info...);
	for (int i = 0; i < n; i++, w *= u / R(i, meta_info...))
		r += w;
	return r;
}

template<typename R,typename Norm, typename ...StructureMetaData>
R exp(R u,const Norm & N=Norm(),real eps=1e-6, StructureMetaData ...meta_info)
{
	R r1(0, meta_info...),r2(1,meta_info...), w=u;
	int k = 1;
	while (N.metric(r1,r2) > eps)
	{
		r1 = r2;
		r2 += w;
		k++;
		w *= u / R(k, meta_info...);
	}
	return r2;
}

template<typename R, typename ...StructureMetaData>
R log_n(R u, int n, StructureMetaData ...meta_info)
{
	R s= u - R(1, meta_info...),r(0,meta_info...), w(1, meta_info...);
	for (int i = 1; i < n; i++, w *= s)
		r += w/R(i%2==1?i:-i,meta_info...);
	return r;
}

template<typename R, typename Norm, typename ...StructureMetaData>
R log(R u, const Norm& N = Norm(), real eps = 1e-6, StructureMetaData ...meta_info)
{
	R s=u-R(1,meta_info...),r1(0, meta_info...), r2=s, w = s*s;
	int k = 2;
	while (N.metric(r1, r2) > eps)
	{
		r1 = r2;
		r2 += w / R(k % 2 == 1 ? k : -k, meta_info...);
		k++;
		w *= s;
	}
	return r2;
}

#endif
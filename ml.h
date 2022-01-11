//
// Created by ramizouari on 01/12/2021.
//

#ifndef __ML_H__
#define __ML_H__
#include "abstract_algebra.h"
#include "optimisation.h"
#include "linear_algebra.h"

/*
* Basic Machine Learning Support
*/
class ml_model
{
	public:
	virtual ml_model& fit(const d_matrix<real> &X, const d_vector<real> &y)=0;
	virtual d_vector<real> predict(const d_matrix<real>& X) const=0;
	virtual real score(const d_matrix<real>& X, const d_vector<real>& y) const = 0;
};

class linear_regression :public ml_model
{
	d_vector<real> w;
public:
	ml_model& fit(const d_matrix<real>& X, const d_vector<real>& y) override
	{
		w = (X.T() * X).solve(X.T() * y);
		return *this;
	}
	d_vector<real> predict(const d_matrix<real>& X) const override
	{
		return X * w;
	}

	real score(const d_matrix<real>& X, const d_vector<real>& y) const override
	{
		auto y_pred = X * w;
		real err = 0;
		for (auto [p, s] : zip(y_pred, y))
			err += pow(p - s, 2);
		return err;
	}
};

class logistic_regression :public ml_model
{
	d_vector<real> w;
public:
	int limit=2000;
	ml_model& fit(const d_matrix<real>& X, const d_vector<real>& y) override
	{
		derivator<d_vector<real>, real, d_vector<real>> D;
		barzilai_borwein_gradient_descent<d_vector<real>, L2_inner_product<real, d_vector<real>>>
			GD(D,1e-3);
		int k = 0;
		w = GD.argmin([&X,&y,&k](const d_vector<real>& u)->std::pair<real,d_vector<real>>
			{
				k++;
				d_vector<real> y_pred = pointwise_function([](auto s)
					{
						return 1 / (1 + std::exp(-s));
					}, X * u),residual=y-y_pred,df{v_shape{(int)X.col_dim()}};
				real err = 0;
					for (auto [p, q] : zip(y, y_pred))
						err += p == 0 ? -std::log(1-q):-std::log(q);
					err /= X.row_dim();
					for (auto [R, r] : zip(X, residual))
						for (auto [g, s] : zip(df, R))
							g -= s * r/X.row_dim();
				return std::make_pair(std::move(err),std::move(df));
			},d_vector<real>{v_shape{(int)X.col_dim()}},limit);
		return *this;
	}

	d_vector<real> boundary_function(const d_matrix<real>& X) const
	{
		return X * w;
	}

	d_vector<real> decision_function(const d_matrix<real>& X) const
	{
		return pointwise_function([](auto s)
			{
				return 1 / (1 + std::exp(-s));
			}, boundary_function(X));
	}

	d_vector<real> predict(const d_matrix<real>& X) const
	{
		auto y = decision_function(X);
		for (auto& s : y)
			s = std::round(s);
		return y;
	}

	real score(const d_matrix<real>& X, const d_vector<real>& y) const override
	{
		auto y_pred = predict(X);
		int a = 0;
		for (auto [s, p] : zip(y, y_pred))
			if (s == p)
				a++;
		return a / (real)y.dim();
	}

	real error(const d_matrix<real>& X, const d_vector<real>& y) const
	{
		auto y_pred = decision_function(X);
		real err = 0;
		for (auto [p, q] : zip(y, y_pred))
			err += p == 0 ? -std::log(1 - q) : -std::log(q);
		return err/y.dim();
	}
};


auto xlogy(auto x, auto y)
{
	if (x == 0)
		return 0;
	return x * std::log(y);
}

class multilogistic_regression :public ml_model
{
	d_matrix<real> W;
public:
	int limit = 2000;
	ml_model& fit(const d_matrix<real>& X, const d_vector<real>& y) override
	{
		d_vector<real> w;
		int C = 0;
		for (auto s : y)
			C = std::max(C, (int)s);
		C++;
		d_matrix<real> Y(0,m_shape{(int)y.dim(),C});
		for (auto [s,R] : zip(y,Y))
			R[(int)s] = 1;
		derivator<d_vector<real>, real, d_vector<real>> D;
		barzilai_borwein_gradient_descent<d_vector<real>, L2_inner_product<real, d_vector<real>>>
			GD(D, 1e-3);
		int k = 0;
		w = GD.argmin([&X, &y, &k,C](const d_vector<real>& u)->real
			{
				k++;
				d_matrix<real> U(0, m_shape{ (int)X.col_dim(),C});
				for (int i = 0; i < X.col_dim(); i++) for (int j = 0; j < C; j++)
					U[i][j] = u[i * C + j];
				d_matrix<real> Y_pred=X*U;
				for (auto& y_pred : Y_pred)
				{
					for (auto& s : y_pred)
						s = std::exp(s);
					real r = std::reduce(y_pred.begin(), y_pred.end());
					for (auto& s : y_pred)
						s /= r;
				}
				d_vector<real> y_pred = pointwise_function([](auto s)
					{
						return 1 / (1 + std::exp(-s));
					}, X * u), residual = y - y_pred, df{ v_shape{(int)X.col_dim()*C} };
					real err = 0;
					for (auto [p, Q] : zip(y, Y_pred))
						err -= std::log(Q[(int)p]);
					err /= X.row_dim();

					return err;
			}, d_vector<real>{v_shape{ (int)X.col_dim()*C }}, limit);
		W = d_matrix<real>(0,m_shape{ (int)X.col_dim(),C });
		for (int i = 0; i < X.col_dim(); i++)
			for (int j = 0; j < C; j++)
				W[i][j] = w[i * C + j];
		return *this;
	}

	d_matrix<real> boundary_function(const d_matrix<real>& X) const
	{
		return X * W;
	}

	d_matrix<real> decision_function(const d_matrix<real>& X) const
	{
		auto Y_pred = boundary_function(X);
		for (auto& y_pred : Y_pred)
		{
			for (auto& s : y_pred)
				s = std::exp(s);
			real r = std::reduce(y_pred.begin(), y_pred.end());
			for (auto& s : y_pred)
				s /= r;
		}
		return Y_pred;
	}

	d_vector<real> predict(const d_matrix<real>& X) const
	{
		auto Y = decision_function(X);
		d_vector<real> y_pred(v_shape{(int)X.row_dim()});
		int C = Y.col_dim();
		for (auto [y,s] : zip(Y,y_pred))
		{
			int k = 0;
			real M = -1;
			for (int i = 0; i < C; i++)
				if (y[i] > M)
				{
					k = i;
					M = y[i];
				}
			s = k;
		}
		return y_pred;
	}

	real score(const d_matrix<real>& X, const d_vector<real>& y) const override
	{
		auto y_pred = predict(X);
		int a = 0;
		for (auto [s, p] : zip(y, y_pred))
			if (s == p)
				a++;
		return a / (real)y.dim();
	}

	real error(const d_matrix<real>& X, const d_vector<real>& y) const
	{
		auto Y_pred = decision_function(X);
		real err = 0;
		for (auto [p, Q] : zip(y, Y_pred))
			err -= std::log(Q[(int)p]);
		err /= X.row_dim();
		return err;
	}
};
#endif //ACPC_PREPARATION_ABSTRACT_ALGEBRA_H

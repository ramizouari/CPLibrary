/*
 * A general implementation of the 2-Phase Simplex algorithm
 * */

#ifndef CPLIBRARY_SIMPLEX2_H
#define CPLIBRARY_SIMPLEX2_H
#include <vector>
#include <algorithm>
#include <concepts>
#include <span>
#include <stdexcept>

namespace cp::topology
{

    template<std::floating_point Float>
    constexpr Float inf = std::numeric_limits<Float>::infinity();

    enum Comparison : int
    {
        LessEqual=-1, Equal, GreaterEqual,
        LE=LessEqual , EQ=Equal, GE=GreaterEqual
    };

    enum ObjectiveDirection : bool
    {
        Min, Max
    };

    template<std::floating_point Float>
    struct LinearConstraint
    {
        LinearConstraint() = default;
        LinearConstraint(const std::vector<int> &v, const std::vector<Float> & coeff, Float target, Comparison cmp = EQ):
                    variables(v),coefficients(coeff),target(target),direction(cmp)
        {

        }

        explicit LinearConstraint(const std::vector<Float> & coeff, Float target, Comparison cmp = EQ) : target(target),coefficients(coeff),direction(cmp)
        {
            auto n =coeff.size();
            for(int i=0;i<n;i++)
                variables.push_back(i);
        }

        std::vector<int> variables;
        std::vector<Float> coefficients;
        Float target{};
        Comparison direction=EQ;
        void addVariable(int v, Float c){
            variables.push_back(v);
            coefficients.push_back(c);
        }
        void setTarget(Float t)
        {
            target= t;
        }
    };

    template<std::floating_point Float>
    struct LinearObjective
    {
        LinearObjective() = default;
        LinearObjective(const std::vector<int> &v, const std::vector<Float> coeff): variables(v),coefficients(coeff){}
        explicit LinearObjective(const std::vector<Float> coeff): coefficients(coeff)
        {
            auto n = coeff.size();
            for(int i=0;i<n;i++)
                variables.push_back(i);
        }

        std::vector<int> variables;
        std::vector<Float> coefficients;
    };

    template<std::floating_point Float,std::integral Int>
    Float max(const std::vector<Float> & Z, Int * pos = nullptr)
    {
        auto it = std::max_element(Z.begin(),Z.end());
        if(pos)
            *pos = std::distance(Z.begin(),it);
        return *it;
    }

    template<std::floating_point Float,std::integral Int>
    Float min(const std::vector<Float> & Z, Int * pos = nullptr)
    {
        auto it = std::min_element(Z.begin(),Z.end());
        if(pos)
            *pos = std::distance(Z.begin(),it);
        return *it;
    }


    enum State
    {
        Infeasible, Unbounded,Optimal
    };

    template<std::floating_point Float>
    void truncate(Float& u, Float eps)
    {
        if(std::abs(u) <= eps)
            u=0;
    }

    template<std::floating_point Float>
    struct SimplexTable
    {
        int n,m;
        std::vector<std::vector<Float>> A;
        std::vector<Float> b,Z;
        Float &W;
        SimplexTable(const std::vector<std::vector<Float>>& A, const std::vector<Float> &b, const std::vector<Float> & Z, Float &W):A(A),b(b),Z(Z),W(W)
        {
            m=b.size();
            n=Z.size()-m;
        }

        void pivot(int p,int q,Float eps)
        {
            for(int j=0;j<n+m;j++) if(j!=q)
                    A[p][j] /= A[p][q];
            b[p]/=A[p][q];
            A[p][q]=1;
            for(int i=0;i<m;i++) if(i!=p)
            {
                auto r=A[i][q];
                for(int j=0;j<n+m;j++)
                    A[i][j] -= A[p][j] * r;
                A[i][q]=0;
                b[i] -= r * b[p];
            }
            auto r = Z[q];
            for(int j=0;j<n+m;j++)
                Z[j] -= r * A[p][j];
            for(int i=0;i<m;i++) for(int j=0;j<n+m;j++)
                truncate(A[i][j],eps);
            for(int i=0;i<m;i++)
                truncate(b[i],eps);
            for(int i=0;i<n+m;i++)
                truncate(Z[i],eps);
            W+=r*b[p];
            Z[q] = 0;
        }

        std::vector<bool> extractSolution(std::span<Float> solution,Float eps) const
        {
            std::vector<bool> basics(n+m);
            for(int j=0;j<n+m;j++)
            {
                int nonZeros=0;
                int k=0;
                for(int i=0;i<m && nonZeros < 2;i++) if(std::abs(A[i][j]) > eps)
                    {
                        nonZeros++;
                        k = i;
                    }
                basics[j] = nonZeros == 1;
                //Basic variable
                if(basics[j])
                    solution[j] = b[k] / A[k][j];
                else
                    solution[j]=0;
            }
            return basics;
        }

        void extractPrimalSolution(std::span<Float> solution,Float eps) const
        {
            for(int j=0;j<n;j++)
            {
                int nonZeros=0;
                int k=0;
                for(int i=0;i<m && nonZeros < 2;i++) if(std::abs(A[i][j]) > eps)
                {
                    nonZeros++;
                    k = i;
                }

                //Basic variable
                if(nonZeros == 1)
                    solution[j] = b[k] / A[k][j];
                else
                    solution[j]=0;
            }
        }

        void extractDualSolution(std::span<Float> solution) const
        {
            for(int j=0;j<m;j++)
                solution[j]= -Z[n+j];
        }
    };


    template<std::floating_point Float>
    struct Simplex
    {
        std::vector<LinearConstraint<Float>> constraints;
        LinearObjective<Float> objective;
        ObjectiveDirection direction=Max;
        int nbrVariables,dualOffset;
        std::vector<Float> primal,dual;
        Float optimal{};
        State state;
        Float eps;
        explicit Simplex(Float eps=1e-6): eps(eps){}


        void addLinearConstraint(const LinearConstraint<Float> & C)
        {
            if(C.direction!=EQ)
                constraints.push_back(C);
            else
            {
                constraints.push_back(LinearConstraint<Float>(C.variables,C.coefficients,C.target,LE));
                auto negCoeffs = C.coefficients;
                for(auto &c:negCoeffs)
                    c=-c;
                constraints.push_back(LinearConstraint<Float>(C.variables,negCoeffs,-C.target,LE));

            }
        }

        void setLinearObjective(const LinearObjective<Float>& z)
        {
            objective = z;
        }

        void setDirection(ObjectiveDirection dir)
        {
            direction=dir;
        }

        State optimize()
        {
            auto tableOpt= prepare();
            if(!tableOpt)
                return state=Infeasible;
            auto & table=*tableOpt;
            auto & [n,m,A,b,Z,_] = table;
            int q;
            while(opt(Z,&q))
            {
                int p=-1;
                std::vector<Float> c(m);
                for(int i=0;i<m;i++)
                    c[i] = b[i] / A[i][q];
                for(int i=0;i<m;i++) if(A[i][q] > 0 && c[i] >= 0 && (p==-1 || c[i] < c[p]))
                    p=i;
                if(p==-1) return state=Unbounded;
                table.pivot(p,q,eps);
            }
            table.extractPrimalSolution(primal,eps);
            table.extractDualSolution(dual);
            return state=Optimal;
        }

        std::span<const Float> primalSolution() const
        {
            return primal;
        }

        std::span<const Float> dualSolution() const
        {
            return dual;
        }

    private:

        template<std::integral I>
        Float opt(const std::vector<Float> & Z, I *pos)
        {
            if(direction == Max) return max(Z,pos) > 0;
            else return min(Z,pos) < 0;
        }


        std::optional<SimplexTable<Float>> prepare()
        {
            updateVarCount();
            std::vector<Float> X;
            auto m = constraints.size();
            auto n = nbrVariables;
            std::vector<std::vector<Float>> A(m,std::vector<Float>(n+m));
            std::vector<Float> b(m),Z(n+m);
            bool originFeasible=true;
            for(int i=0;i<m;i++)
            {
                auto &[vars,coeffs,target,dir] = constraints[i];
                auto r = constraints[i].variables.size();
                for(int j=0;j<r;j++)
                    A[i][vars[j]] = coeffs[j];
                switch (dir)
                {
                    case LE:
                        A[i][n+i]=1;
                        if(target <0) originFeasible=false;
                        b[i]=target;
                        break;
                    case GE:
                        if(target>0) originFeasible=false;
                        b[i]=-target;
                        for(int j=0;j<n;j++) A[i][j]=-A[i][j];
                        A[i][n+i]=1;
                        break;
                    default:
                        throw std::runtime_error("Equal constraints should be preprocessed");
                }
            }
            auto & [vars,coeffs] = objective;
            auto r=vars.size();
            for(int i=0;i<r;i++)
                Z[vars[i]] = coeffs[i];
            primal.resize(n);
            dual.resize(m);
            if(originFeasible)
                return SimplexTable<Float>{A,b,Z,optimal};
            else
            {
                Float W=0,t=0;
                SimplexTable<Float> table{A,b,Z,W};
                ++table.n;
                for(auto &a:table.A)
                    a.push_back(-1);
                std::fill(table.Z.begin(),table.Z.end(),0);
                int k=0;
                // Find an initial pivot of the auxiliary LP
                for(int i=0;i<m;i++) if(table.b[i] < t)
                {
                    k=i;
                    t=table.b[i];
                }
                table.Z.push_back(1);
                table.pivot(k,n+m,eps);
                // Apply simplex to the auxiliary LP
                int q;
                while(min(table.Z,&q) <0)
                {
                    int p=-1;
                    std::vector<Float> c(m);
                    for(int i=0;i<m;i++)
                        c[i] = table.b[i] / table.A[i][q];
                    for(int i=0;i<m;i++) if(table.A[i][q] > 0 && c[i] >= 0 && (p==-1 || c[i] < c[p]))
                        p=i;
                    if(p==-1) throw std::runtime_error("Auxiliary LP cannot be unbounded");
                    table.pivot(p,q,eps);
                }
                std::vector<Float> sol(n+m+1);
                auto basics=table.extractSolution(sol,eps);
                // If the auxiliary objective is positive, then the original LP is infeasible
                if(W > eps)
                    return std::nullopt;
                // If the artificial variable is basic, apply a degenerate pivot
                if(basics[n+m])
                {
                    int p;
                    for(p=0;p<m; p++) if(std::abs(table.A[p][n+m]) > eps)
                            break;
                    for(q=0;q < n+m;q++) if(std::abs(table.A[p][q]) > eps && !basics[q])
                            break;
                    table.pivot(p,q,eps);
                    basics=table.extractSolution(sol,eps);
                }
                // Transform the original LP to the basis found by the auxiliary LP
                for(int j=0;j<table.n+table.m;j++) if(basics[j]) for(int i=0;i<table.m;i++) if(std::abs(table.A[i][j]) > eps)
                {
                    for(int k=0;k<n+m;k++) if(k!=j)
                        Z[k]-=table.A[i][k]*Z[j] / table.A[i][j];
                    optimal+=table.b[i]*Z[j] / table.A[i][j];
                    if(j<n+m) Z[j]=0;
                    break;
                }
                for(auto & a:table.A) a.pop_back();
                return SimplexTable<Float>{table.A,table.b,Z,optimal};
            }
        }

        void updateVarCount()
        {
            auto it = std::max_element(objective.variables.begin(),objective.variables.end());
            nbrVariables = it==objective.variables.end()? 0 : *it+1;
            for(const auto & [V,_1,_2,_3] : constraints)
            {
                auto maxCount = std::max_element(V.begin(),V.end());
                nbrVariables=std::max<int>(nbrVariables,maxCount == V.end()? 0 : *maxCount+1);
            }
        }
    };
}


#endif //CPLIBRARY_SIMPLEX2_H

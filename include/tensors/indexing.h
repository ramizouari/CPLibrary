//
// Created by ramizouari on 07/01/2026.
//

#ifndef CPLIBRARY_INDEXING_H
#define CPLIBRARY_INDEXING_H
#include <generator>
#include <ranges>

namespace cp {
    namespace indexing {
        template<typename Range>
        concept sized_random_access = std::ranges::random_access_range<Range> && std::ranges::sized_range<Range>;

        constexpr std::size_t extent_concat_dim(std::size_t extent1,std::size_t extent2) {
            if (extent1 == dynamic_extent || extent2 == dynamic_extent) return dynamic_extent;
            return extent1+extent2;
        }

        template<std::size_t extent>
        struct index_vector : std::array<std::size_t,extent>
        {
            using std::array<std::size_t,extent>::array;

            auto &operator+=(const index_vector &other) {
                for (int i=0;i<extent;++i) (*this)[i]+=other[i];
                return *this;
            }

            auto operator+(const index_vector &other) const {
                index_vector result(*this);
                return result+=other;
            }

            auto& operator-=(const index_vector &other) {
                for (int i=0;i<extent;++i) (*this)[i]-=other[i];
                return *this;
            }

            auto operator-(const index_vector &other) const {
                index_vector result(*this);
                return result-=other;
            }

            auto& operator*=(const index_vector &other) {
                for (int i=0;i<extent;++i) (*this)[i]*=other[i];
                return *this;
            }

            auto operator*(const index_vector &other) const {
                index_vector result(*this);
                return result*=other;
            }

            auto& operator|=(const index_vector &other) const {
                for (int i=0;i<extent;++i) (*this)[i]|=other[i];
                return *this;
            }

            auto operator|(const index_vector &other) const {
                index_vector result(*this);
                return result|=other;
            }

            auto operator&=(const index_vector &other) const {
                for (int i=0;i<extent;++i) (*this)[i]&=other[i];
                return *this;
            }

            auto operator&(const index_vector &other) const {
                index_vector result(*this);
                return result&=other;
            }

            auto operator~() const {
                index_vector result(*this);
                for (auto &r:result) r=~r;
                return result;
            }

            auto& operator<<=(std::size_t n)  {
                for (auto &r:*this) r<<=n;
                return *this;
            }
            auto operator<<(std::size_t n) const {
                index_vector result(*this);
                return result<<=n;
            }

            auto operator>>=(std::size_t n)  {
                for (auto &r:*this) r>>=n;
                return *this;
            }
            auto operator>>(std::size_t n) const {
                index_vector result(*this);
                return result>>=n;
            }
        };

        template<std::size_t ext1, std::size_t ext2>
        index_vector<extent_concat_dim(ext1,ext2)> concat(const index_vector<ext1> &A, const index_vector<ext2> &B)
        {
            constexpr auto ext12=extent_concat_dim(ext1,ext2);
            index_vector<ext12> result;
            if constexpr (ext12 == dynamic_extent)
                result.resize(A.size()+B.size());
            std::copy(A.begin(),A.end(),result.begin());
            std::copy(B.begin(),B.end(),result.begin()+result.size());
            return result;
        }

        template<sized_random_access P,sized_random_access Q>
        bool increment(P &X, const Q &shape)
        {
            for(int i=0;i<shape.size();i++)
            {
                ++X[i];
                if(X[i]<shape[i]) return false;
                X[i]=0;
            }
            return true;
        }

        template<sized_random_access P,sized_random_access Q>
        bool rincrement(P &X, const Q &shape)
        {
            for(int i=shape.size()-1;i>=0;--i)
            {
                ++X[i];
                if(X[i]<shape[i]) return false;
                X[i]=0;
            }
            return true;
        }

        template<sized_random_access P,sized_random_access Q>
        bool decrement(P &X, const Q &shape)
        {
            for(int i=shape.size()-1;i>=0;--i)
            {
                if (X[i] > 0) {
                    --X[i];
                    return false;
                }
                X[i]=shape[i]-1;
            }
            return true;
        }


        template<sized_random_access P,sized_random_access Q>
        bool rdecrement(P &X, const Q &shape)
        {
            for(int i=0;i<shape.size();i++)
            {
                if (X[i] > 0) {
                    --X[i];
                    return false;
                }
                X[i]=shape[i]-1;
            }
            return true;
        }

        enum class iterate_direction {
            right_up,right_down,left_up,left_down
        };

        inline bool is_up_direction(iterate_direction dir) {
            return dir==iterate_direction::right_up || dir==iterate_direction::left_up;
        }


        template<sized_random_access Q>
        class multi_index_iterator {
            const Q shape;
            Q index;
            const iterate_direction iterate_dir;
            bool is_end=false;
        public:
            explicit multi_index_iterator(
                const Q shape,
                const iterate_direction iterate_dir = iterate_direction::right_up
                ):
                shape(std::move(shape)), index(shape), iterate_dir(iterate_dir)
            {
                if (is_up_direction(iterate_dir))
                    std::fill(index.begin(),index.end(),0);
                else for (auto &i:index) --i;
            }

            explicit multi_index_iterator(
                const Q shape,
                Q index,
                const iterate_direction iterate_dir = iterate_direction::right_up
                ):
                shape(std::move(shape)), index(std::move(index)), iterate_dir(iterate_dir)
            {
            }

            explicit multi_index_iterator(
                const Q shape,
                bool is_end,
                const iterate_direction iterate_dir = iterate_direction::right_up
                ):
                multi_index_iterator(shape,iterate_dir)
            {
                this->is_end=is_end;
            }

            auto& operator*() const {
                return index;
            }

            auto& operator*() {
                return index;
            }

            auto& operator->() {
                return index;
            }

            const auto& operator->() const {
                return index;
            }

            bool operator==(const multi_index_iterator& other) const {
                return is_end == other.is_end && (is_end || index == other.index);
            }

            multi_index_iterator& operator++() {
                switch (iterate_dir) {
                    case iterate_direction::left_up:
                        is_end=increment(index,shape);
                        break;
                    case iterate_direction::left_down:
                        is_end=decrement(index,shape);
                        break;
                    case iterate_direction::right_up:
                        is_end=rincrement(index,shape);
                        break;
                    case iterate_direction::right_down:
                        is_end=rdecrement(index,shape);
                }
                return *this;
            }

            multi_index_iterator operator++(int) {
                auto tmp=*this;
                ++*this;
                return tmp;
            }
        };

        template<sized_random_access Q>
        struct iterated_multi_index {
            Q shape;
            iterate_direction direction = iterate_direction::right_up;
            auto begin() const {return multi_index_iterator<Q>(shape,direction);}
            auto end() const {return multi_index_iterator<Q>(shape,true,direction);}

        };

        template<sized_random_access Q>
        std::generator<Q> iterate_multi_index_right_up(const Q &shape)
        {
            Q I=shape;
            std::fill(I.begin(),I.end(),0);
            bool is_end = I.size() == 0;
            while (!is_end) {
                co_yield I;
                is_end = rincrement(I,shape);
            }
        }

        template<sized_random_access Q>
        std::generator<Q> iterate_multi_index_left_up(const Q &shape)
        {
            Q I=shape;
            std::fill(I.begin(),I.end(),0);
            bool is_end = I.size() == 0;
            while (!is_end) {
                co_yield I;
                is_end = increment(I,shape);
            }
        }
        template<sized_random_access Q>
        std::generator<Q> iterate_multi_index_right_down(const Q &shape)
        {
            Q I(shape);
            for (auto &i:I) --i;
            bool is_end = I.size() == 0;
            while (!is_end) {
                co_yield I;
                is_end = rdecrement(I,shape);
            }
        }

        template<sized_random_access Q>
        std::generator<Q> iterate_multi_index_left_down(const Q &shape)
        {
            Q I(shape);
            for (auto &i:I) --i;
            bool is_end = I.size() == 0;
            while (!is_end) {
                co_yield I;
                is_end = decrement(I,shape);
            }
        }
    }

    using namespace indexing;
    namespace tensors {
        using namespace indexing;
    }
}

#endif //CPLIBRARY_INDEXING_H
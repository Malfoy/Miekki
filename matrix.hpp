#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <memory>
#include <algorithm>
#include "common.hpp"

namespace matrix {

using idx_t = int; // Signed indices allows better loop optimizations

struct slice {
    idx_t offset;
    idx_t size;
};

struct strided_slice {
    idx_t offset;
    idx_t size;
    idx_t stride;
};


template<typename T> struct slicevec;

/// Reference to a subarray
template<typename T>
struct refvec {
    using idx_t = int;
    using scalar_t = T;
    using ref_t = scalar_t&;
    using const_ref_t = scalar_t&; // Constness is contained in T if needed
    using iterator = scalar_t*;
    using const_iterator = scalar_t*;

    refvec(scalar_t* data ,idx_t size) : _data(data), _size(size) {}
    refvec(refvec&&) = default;
    refvec(const refvec&) = default;

    scalar_t* data() const { return _data; }
    iterator begin() const { return _data; }
    iterator end() const { return _data + _size; }

    ref_t operator[](idx_t i) const {
        assume(i >= 0, "negative indice");
        assume(i < _size, "out of bound indice");
        return _data[i];
    }

    refvec<T> operator[](slice s) {
        assume(s.offset >= 0, "negative offset");
        assume(s.offset + s.size <= _size, "Out of bound slice");
        return refvec<T>(_data + s.offset, s.size);
    }

    slicevec<T> operator[](strided_slice s) {
        assume(s.offset >= 0, "negative offset");
        idx_t min = s.offset;
        idx_t max = s.offset + (s.size-1)*s.stride;
        assume(min >= 0 && min < _size && max >= 0 && max < _size, "Slice out of bound");
        return slicevec<T>(_data, s);
    }

    idx_t size() const { return _size; }

    scalar_t* _data;
    idx_t _size;
};


/// Reference to a subarray with strides
template<typename T>
struct slicevec : protected refvec<T> {
    using scalar_t = T;
    using ref_t = scalar_t&;

    using refvec<T>::size;

    // Refereced memory size: this->_size * this->_stride

    slicevec(scalar_t* data , strided_slice s) : refvec<T>(data + s.offset, s.size), _stride(s.stride) {}
    slicevec(slicevec&&) = default;
    slicevec(const slicevec&) = default;

    idx_t stride() const { return _stride; }

    ref_t operator[](idx_t i) const {
        assume(i >= 0, "negative indice");
        assume(i < this->_size, "out of bound indice");
        return this->_data[i*_stride];
    }

    slicevec<T> operator[](slice s) {
        assume(s.offset >= 0, "negative offset");
        assume(s.offset + s.size <= this->_size, "Out of bound slice");
        return slicevec<T>(this->_data, {s.offset * _stride, s.size, _stride});
    }

    slicevec<T> operator[](strided_slice s) {
        assume(s.offset >= 0, "negative offset");
        idx_t min = s.offset;
        idx_t max = s.offset + (s.size-1)*s.stride;
        assume(min >= 0 && min < this->_size && max >= 0 && max < this->_size, "Slice out of bound");
        return slicevec<T>(this->_data, {s.offset * _stride, s.size, s.stride * _stride});
    }


    struct iterator {
        scalar_t* ptr;
        idx_t stride;

        iterator& operator++() { ptr += stride; return *this; }
        bool operator<(const iterator& other) { return ptr < other.ptr; }
        bool operator!=(const iterator& other) { return ptr != other.ptr; }
        ref_t operator*() const { return *ptr; }
    };

    using const_iterator = iterator;

    iterator begin() const { return {this->data(), _stride}; }
    iterator end() const { return {this->data() + this->size() * _stride, _stride}; }

    idx_t _stride;
};


/// Own a memory range of T values
template<typename T>
struct vector {
    using idx_t = int;
    using scalar_t = T;
    using ref_t = scalar_t&;
    using const_ref_t = const scalar_t&;
    using iterator = scalar_t*;
    using const_iterator = const scalar_t*;
    using unique_arr_t = std::unique_ptr<scalar_t[]>;
    typedef scalar_t (&array_ref)[];
    typedef const scalar_t (&const_array_ref)[];

    vector(idx_t size) : _data(new scalar_t[size_t(size)]), _size(size) {}
    vector(idx_t size, scalar_t v) : _data(new scalar_t[size_t(size)]{v}), _size(size) {}
    vector(unique_arr_t&& data, idx_t size) : _data(std::move(data)), _size(size) {}
    vector(iterator* data, idx_t size) : vector(size) {
        std::copy(data, data + _size);
    }
    vector(vector&&) = default;
    vector& operator=(vector&&) = default;
    explicit vector(const vector& from) : vector(from._data.get(), from._size) {}

    array_ref data() { return _data.get(); }
    const_array_ref data() const { return _data.get(); }
    idx_t size() const {return _size; }

    unique_arr_t as_unique_ptr() && {
        return std::move(_data);
    }

    operator refvec<T>() { return {_data.get(), _size}; }
    operator refvec<const T>() const { return {_data.get(), _size}; }

    template<typename I> auto operator[](I i) -> decltype(this->operator refvec<T>()[i])
    { return this->operator refvec<T>()[i]; }
    template<typename I> auto operator[](I i) const -> decltype(this->operator refvec<const T>()[i])
    { return *this->operator refvec<const T>()[i]; }

    iterator begin() { return _data.get(); }
    iterator end() { return _data.get() + _size; }
    const_iterator begin() const { return _data.get(); }
    const_iterator end() const { return _data.get() + _size; }

    unique_arr_t _data;
    idx_t _size;
};

template<typename base>
struct matrix_t : public base {
    using typename base::idx_t;
    using typename base::scalar_t;

    matrix_t(idx_t r, idx_t c) : base(r*c), _rows(r), _cols(c) {}
    matrix_t(idx_t r, idx_t c, scalar_t v) : base(r*c, v), _rows(r), _cols(c) {}
//    template<typename T, typename=decltype(base(std::declval<T>()), idx_t{})>
//    matrix(T&& x, idx_t r, idx_t c) : base(std::forward<T>(x), r * c), _rows(r), _cols(c) {}

    idx_t cols() const { return _cols; }
    idx_t rows() const { return this->size() / _cols; }

    auto row(idx_t i) -> decltype((*this)[slice{}]) { return (*this)[slice{_cols * i, _cols}]; }
    auto row(idx_t i) const -> decltype((*this)[slice{}]) { return (*this)[slice{_cols * i, _cols}]; }

    auto col(idx_t i) -> decltype((*this)[strided_slice{}]) { return (*this)[strided_slice{i, _rows, _cols}]; }
    auto col(idx_t i) const -> decltype((*this)[strided_slice{}]) { return (*this)[strided_slice{i, _rows, _cols}]; }

    auto operator() (idx_t r, idx_t c) -> decltype((*this)[idx_t{}]) { return static_cast<refvec<scalar_t>>(*this)[r*_cols + c]; }
    auto operator() (idx_t r, idx_t c) const -> decltype((*this)[idx_t{}]) { return this[r*_cols + c]; }

    idx_t _rows, _cols;
};

template<typename T> using matrix = matrix_t<vector<T>>;
template<typename T> using matrix_ref = matrix_t<refvec<T>>;

} /* namespace matrix */

#endif // MATRIX_HPP

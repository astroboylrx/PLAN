//
//  global.hpp
//  PLAN: PLantesimal ANalyzer
//
//  Created by Rixin Li on 4/26/16.
//  Copyright © 2016 Rixin Li. All rights reserved.
//

/*! \file global.hpp
 *  \brief provide library headers, I/O-related class and utilities */

#ifndef global_hpp
#define global_hpp

// Include C libraries first
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cassert>
#include <algorithm>
// Include C++ libraries
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <bitset>
#include <iomanip>
#include <numeric>
#include <array>
#include <type_traits>
// Include other libraries
#include <unistd.h>
#include <getopt.h>

#define MPI_ON // Comment out this line before committing!!
#ifdef MPI_ON // "ifdef" options are defined during compilation
#include "mpi.h"
#endif // MPI_ON

// Check c++11 support
#if __cplusplus <= 199711L
#error This program wants a C++11 compliant compiler (option -DOLDCPP which supports old compilers has been abandoned).
#endif // __cplusplus

/***********************************/
/********** SmallVec Part **********/
/***********************************/

/*
 * N.B.: template classes need to have the method definitions inside the header file, in order to let the linker work. Alternatively, you can put definitions in other file and include after declaration in the header file.
 */

/*! \struct template <bool, class T, class U> __SelectIf_base
 *  \brief a template for struct __SelectIF to accept a bool as condition */
template <bool, class T, class U>
struct __SelectIf {};

/*! \struct template <class T, class U> __SelectIf_true
 *  \brief a template for struct __SelectIF to select T type if bool is true */
template <class T, class U>
struct __SelectIf<true, T, U> { typedef T type; };

/*! \struct template <class T, class U> __SelectIf_false
 *  \brief a template for struct __SelectIF to select U type if bool is false */
template <class T, class U>
struct __SelectIf<false, T, U> { typedef U type; };

/*! \struct template <class T, class U> PromoteNumeric
 * this struct controls the type promotion, this struct nests many levels of selections. Read comments/explanations from inside and notice that comments start with new lines are different with comments after statements */
template <class T, class U>
struct PromoteNumeric {
    typedef typename __SelectIf<
    //if T and U are both integers or both non-integers
    std::numeric_limits<T>::is_integer == std::numeric_limits<U>::is_integer, // middle, served as bool for outermost
    //then pick the larger type
    typename __SelectIf<(sizeof(T) > sizeof(U)), T,
    //else if they are equal
    typename __SelectIf<(sizeof(T) == sizeof(U)),
    //pick the one which is unsigned
    typename __SelectIf<std::numeric_limits<T>::is_signed, U, T>::type, // this is the innermost layer
    //else select U as bigger
    U
    >::type // this is the second innermost layer
    >::type, // middle, served as T for outermost, nests other layers
    //else pick the one which is not integer
    typename __SelectIf<std::numeric_limits<T>::is_integer, U, T>::type // middle, served as U for outermost
    >::type type; // outermost layer
};

/*! \class template <class T, int D> SmallVec
 *  \brief define a vector class and all related arithmetic computations
 *  \tparam T data type of elements of this vector
 *  \tparam D dimension of this vector */
template <class T, int D>
class SmallVec {
private:
    
public:
    /*! \var T data[D]
     *  \brief data of elements of this vector */
    T data[D];
    
    /*! \fn SmallVec() : data{0} {}
     *  \brief the most basic constructor, list initialization with {0} */
    SmallVec() : data{0} {}
    
    /*! \fn template <class U> SmallVec(const SmallVec<U, D> &initializer)
     *  \brief implicit copy constructor, read explanations below */
    template <class U>
    SmallVec(const SmallVec<U, D> &initializer) {
        for (int i = 0; i != D; i++) {
            data[i] = static_cast<T>(initializer.data[i]);
        }
    }
    
    /*
     * If explicit keyword is applied to a copy constructor, it means that object of that class can't be copied when being passed tofunctions or when being returned from function - (this type of copying is called implicit copying). So statement "SmallVec<T, D> new_vec = old_vec;" will cause error during compilation. Passing SmallVec as an function argument by value or make it the return type of a function will also cause error, like "Func(old_vec)" or "SmallVec<T, D> Func()". Only explicit copying, i.e., "SmallVec<T, D> new_vec (old_vec);" is allowed. And only being passed to a function (or being returned from a funtion) by reference/pointer is allowed.
     */
    
    /*! \fn template <class U, typename std::enable_if<!(std::is_pointer<U>::value), int>::type = 0> explicit SmallVec(const U& scalar)
     *  \brief overloading constructor to broadcast a scalar, e.g., unit_vec = SmallVec<double, 3>(1.0). For keyword "explicit", read explanations below */
    template <class U, typename std::enable_if<!(std::is_pointer<U>::value), int>::type = 0>
    explicit SmallVec(const U& scalar) {
        for (int i = 0; i != D; i++) {
            data[i] = static_cast<T>(scalar);
        }
    }
    
    /*! \fn template <class U, typename std::enable_if<!(std::is_pointer<U>::value), int>::type = 0> explicit SmallVec(const U (&vec)[D])
     *  \brief overloading constructor to copy an old-fashion array */
    template <class U, typename std::enable_if<!(std::is_pointer<U>::value), int>::type = 0>
    explicit SmallVec(const U (&vec)[D]) { // remember how to reference an array
        for (int i = 0; i != D; i++) {
            data[i] = static_cast<T>(vec[i]);
        }
    }
    
    /*! \fn template <class U, size_t E, typename std::enable_if<!(std::is_pointer<U>::value), int>::type = 0, typename std::enable_if<E != D, int>::type = 0> explicit SmallVec(const U (&vec)[E])
     *  \brief overloading constructor to notify user a wrong old-fashion array is input */
    template <class U, size_t E, typename std::enable_if<!(std::is_pointer<U>::value), int>::type = 0, typename std::enable_if<E != D, int>::type = 0>
    explicit SmallVec(const U (&vec)[E]) { // remember how to reference an array
        std::cerr << "Error: Use an array with wrong size to construct class SmallVec." << std::endl;
        exit(4); // wrong function argument
    }
    
    /*! \fn template <class U, typename std::enable_if<!(std::is_pointer<U>::value), int>::type = 0> explicit SmallVec(const std::array<U, D> &vec)
     *  \brief overloading constructor to copy an STL array */
    template <class U, typename std::enable_if<!(std::is_pointer<U>::value), int>::type = 0>
    explicit SmallVec(const std::array<U, D> &vec) {
        for (int i = 0; i != D; i++) {
            data[i] = static_cast<T>(vec[i]);
        }
    }
    
    /*! \fn template <class U, size_t E, typename std::enable_if<!(std::is_pointer<U>::value), int>::type = 0, typename std::enable_if<E != D, int>::type = 0> explicit SmallVec(const std::array<U, E> &vec)
     *  \brief overloading constructor to notify user a wrong STL array is input */
    template <class U, size_t E, typename std::enable_if<!(std::is_pointer<U>::value), int>::type = 0, typename std::enable_if<E != D, int>::type = 0>
    explicit SmallVec(std::array<U, E> const &vec) {
        std::cerr << "Error: Use an array with wrong size to construct class SmallVec." << std::endl;
        exit(4); // wrong function argument
    }
    
    template <class U, typename std::enable_if<std::is_pointer<U>::value, int>::type = 0>
    SmallVec(U arg) {
        std::cerr << "Error: Please don't provide a pointer to construct class SmallVec. Whether it is a pointer to scalar or a pointer to array remains unknown. " << std::endl;
        exit(4); // wrong function argument
    }
    
    /*
     * Why do we use the keyword "explicit"? (reference: http://www.geeksforgeeks.org/g-fact-93/ )
     * In C++, if a class has a constructor which can be called with a single argument, then this constructor becomes conversion constructor because such a constructor allows conversion of the single argument to the class being constructed.
     * We can avoid such implicit conversions as these may lead to unexpected results. We can make the constructor explicit with the help of explicit keyword.
     */
    
    /*! \fn template <typename... Tail> SmallVec(typename std::enable_if<sizeof...(Tail)+1==D, T>::type head, Tail... tail) : data{ head, T(tail)... }
     *  \brief overloading constructor for initialization in the form of e.g., SmallVec<double, 3>(3, 4.5, 5) */
    template <typename... Tail>
    SmallVec(typename std::enable_if<sizeof...(Tail)+1==D, T>::type head, Tail... tail) : data{head, T(tail)...} {}
    /*
     * C++ has a special parameter type, ellipsis, that can be used to pass a varying number of arguments. It is called the ellipsis operator. Try to undersatnd it by recalling the usage of printf() in C where you can input any amount of arguments. Also, sizeof...() is used to count the number of arguments, which is different with sizeof()(reference: http://www.cplusplus.com/articles/EhvU7k9E/ )
     * Then "data{ head, T(tail)... }" is called an initialization list. Usually it is for passing arguments to the constructor of a parent class.
     */
    
    /*! \fn T operator[] ( const size_t i ) const
     *  \brief allow access with usual vector notation */
    T operator[] (const size_t i) const {
        assert(i < D);
        return *(data+i);
    }
    /*
     * Here "const" after a function declaration means that the function is not allowed to change any class members (except ones that are marked "mutable").
     */
    
    /*! \fn T& operator [] ( const size_t i )
     *  \brief allow access with usual vector notation */
    T& operator [] (const size_t i) {
        assert(i < D);
        return *(data+i);
    }
    
    /*! \fn const SmallVec& operator +()
     *  \brief unary operations (sign): make +SmallVec return itself */
    const SmallVec& operator +() {
        return *this;
    }
    
    /*! \fn SmallVec operator -()
     *  \brief unary operations (sign): make -SmallVec return -1*data */
    SmallVec operator -() {
        SmallVec<T, D> result;
        for (int i = 0; i != D; i++) {
            result[i] = -data[i];
        }
        return result;
    }
    
    /*! \fn template <class U> SmallVec<T, D>& operator = (const SmallVec<U, D>& rhs)
     *  \brief assignment operator = */
    template <class U>
    SmallVec<T, D>& operator = (const SmallVec<U, D>& rhs) {
        for (int i = 0; i != D; i++) {
            data[i] = static_cast<T>(rhs.data[i]);
        }
        return *this;
    }
    
    /*! \fn template <class U> SmallVec<T, D>& operator += (const SmallVec<U, D>& rhs)
     *  \brief compound assignment operator += */
    template <class U>
    SmallVec<T, D>& operator += (const SmallVec<U, D>& rhs) {
        for (int i = 0; i != D; i++) {
            data[i] += rhs.data[i];
        }
        return *this;
    }
    
    /*! \fn template <class U> SmallVec<T, D>& operator -= (const SmallVec<U, D>& rhs)
     *  \brief compound assignment operator -= */
    template <class U>
    SmallVec<T, D>& operator -= (const SmallVec<U, D>& rhs) {
        for (int i = 0; i != D; i++) {
            data[i] -= rhs.data[i];
        }
        return *this;
    }
    
    /*! \fn template <class U> SmallVec<T, D>& operator *= (const U rhs)
     *  \brief compound assignment operator *= */
    template <class U>
    SmallVec<T, D>& operator *= (const U rhs) {
        for (int i = 0; i != D; i++) {
            data[i] *= rhs;
        }
        return *this;
    }
    
    /*! \fn template <class U> SmallVec<T, D>& operator /= (const U rhs)
     *  \brief compound assignment operator /= */
    template <class U>
    SmallVec<T, D>& operator /= (const U rhs) {
        for (int i = 0; i != D; i++) {
            data[i] /= rhs;
        }
        return *this;
    }
    
    /*! \fn typename PromoteNumeric<T, double>::type Norm() const
     *  \brief calculate the norm of this vector, |x|, for at least double precision */
    typename PromoteNumeric<T, double>::type Norm() const {
        typename PromoteNumeric<T, double>::type sum = 0;
        for (int i = 0; i != D; i++) {
            sum += data[i] * data[i];
        }
        return sqrt(sum);
    }
    
    /*! \fn typename PromoteNumeric<T, double>::type Norm2() const
     *  \brief calculate the norm^2 of this vector, |x|^2, for at least double precision */
    typename PromoteNumeric<T, double>::type Norm2() const {
        typename PromoteNumeric<T, double>::type sum = 0;
        for (int i = 0; i != D; i++) {
            sum += data[i] * data[i];
        }
        return sum;
    }
    
    /*! \fn T MaxElement() const
     *  \brief return the maximum element of this vector */
    inline T MaxElement() const {
        return *std::max_element(data, data+D);
    }
    
    /*! \fn T MinElement() const
     *  \brief return the minimum element of this vector */
    inline T MinElement() const {
        return *std::min_element(data, data+D);
    }
    
    /*! \fn template<class U> typename PromoteNumeric<T, U>::type Dot(const SmallVec<U, D>& rhs) const
     *  \brief calculate the dot product with another vector rhs (with the same dimension) */
    template<class U>
    typename PromoteNumeric<T, U>::type Dot(const SmallVec<U, D>& rhs) const {
        typename PromoteNumeric<T, U>::type tmp = data[0] * rhs[0];
        for (int i = 1; i < D; i++) {
            tmp += data[i] * rhs[i];
        }
        return tmp;
    }
    
    /*! \fn template <class U> inline SmallVec<typename PromoteNumeric<T, U>::type, 3> Cross(const SmallVec<U, 3>& rhs) const
     *  \brief calculate cross product of this_vector x rhs in 3-dimension-space */
    template <class U>
    inline SmallVec<typename PromoteNumeric<T, U>::type, 3> Cross(const SmallVec<U, 3>& rhs) const {
        SmallVec<typename PromoteNumeric<T, U>::type, 3> tmp;
        tmp[0] = data[1] * rhs[2] - data[2] * rhs[1];
        tmp[1] = data[2] * rhs[0] - data[0] * rhs[2];
        tmp[2] = data[0] * rhs[1] - data[1] * rhs[0];
        return tmp;
    }
    
    /*! \fn template <class U, class V> bool AbsClose(const SmallVec<U, D>& rhs, const V epsilon) const
     *  \brief determine if the absolute difference between this_vector and rhs is less than epsilon in each element */
    template <class U, class V>
    bool AbsClose(const SmallVec<U, D>& rhs, const V epsilon) const {
        SmallVec<typename PromoteNumeric<T, U>::type, D> diff;
        diff = *this - rhs;
        bool val = true;
        for (int i = 0; i != D; i++) {
            val = val && fabs(diff[i] < epsilon);
        }
        return val;
    }
    
    /*! \fn template <class U, class V> int relclose(const SmallVec<U, D>& rhs, const V epsilon) const
     *  \brief determinte (the relative difference between this_vector and rhs) divided by (the average length of them) is less than epsilon in each element */
    template <class U, class V>
    int RelClose(const SmallVec<U, D>& rhs, const V epsilon) const {
        SmallVec<typename PromoteNumeric<T, U>::type, D> sum, diff;
        for (int i = 0; i != D; i++) {
            sum[i] = fabs(data[i]) + fabs(rhs[i]);
        }
        diff = *this - rhs;
        bool val = true;
        for (int i = 0; i != D; i++) {
            val = val && ((2 * fabs(diff[i]) / sum[i]) < epsilon);
        } // if sum[i] = 0, we will get nan
        return val;
    }
    
    /*! \fn bool IsFinite() const
     *  \brief determine if every element is finite */
    bool IsFinite() const {
        bool val = true;
        for (int i = 0; i != D; i++) {
            val = val && std::isfinite(data[i]);
        }
        return val;
    }
    
    /*! \fn template <class U> bool operator == (const SmallVec<U, D>& rhs) const
     *  \brief overloading operator == */
    template <class U>
    bool operator == (const SmallVec<U, D>& rhs) const {
        bool val = true;
        for (int i = 0; i != D; i++) {
            val = val && (data[i] == rhs[i]);
        }
        return val;
    }
    
    /*! \fn template <class U> bool operator != (const SmallVec<U, D>& rhs) const
     *  \brief overloading operator != */
    template <class U>
    bool operator != (const SmallVec<U, D>& rhs) const {
        return !(*this == rhs);
    }
    
    /*! \fn template <class U> inline bool operator < (const SmallVec<U, D>& rhs) const
     *  \brief overloading operator <, this is no point to compare each element. We compare the norm of vectors. */
    template <class U>
    inline bool operator < (const SmallVec<U, D>& rhs) const {
        return (this->Norm() < rhs.Norm());
    }
    
    /*! \fn template <class U> bool operator <= (const SmallVec<U, D>& rhs) const
     *  \brief overloading operator <= */
    template <class U>
    bool operator <= (const SmallVec<U, D>& rhs) const {
        if ((*this < rhs) || (*this == rhs)) return true;
        return false;
    }
    
    /*! \fn template <class U> inline bool operator > (const SmallVec<U, D>& rhs) const
     *  \brief overloading operator >, this is no point to compare each element. We compare the norm of vectors. */
    template <class U>
    inline bool operator > (const SmallVec<U, D>& rhs) const {
        return (this->Norm() > rhs.Norm());
    }
    
    /*! \fn template <class U> bool operator >= (const SmallVec<U, D>& rhs) const
     *  \brief overloading operator >= */
    template <class U>
    bool operator >= (const SmallVec<U, D>& rhs) const {
        if((*this>rhs) || (*this==rhs)) return true;
        return false;
    }
    
    /*! \fn SmallVec<T, D> SetZeros()
     *  \brief reset data to {0} */
    SmallVec<T, D> SetZeros() {
        for (int i = 0; i != D; i++) {
            data[i] = 0;
        }
        return *this;
    }
    
    /*! \fn template <class U, class V> bool InRange(U low, V high)
     *  \brief return true if all elements are on [a, b] */
    template <class U, class V>
    bool InRange(U low, V high) {
        bool val = true;
        for (int i = 0; i != D; i++) {
            val = val && (data[i] >= low && data[i] <= high);
        }
        return val;
    }
    
    /*! \fn template <class U, class V> bool InRange(SmallVec<U, D> low, SmallVec<V, D> high)
     *  \brief overloading InRange by taking SmallVec-type arguments instead of arithmetic type */
    template <class U, class V>
    bool InRange(SmallVec<U, D> low, SmallVec<V, D> high) {
        bool val = true;
        for (int i = 0; i != D; i++) {
            val = val && (data[i] >= low[i] && data[i] <= high[i]);
        }
        return val;
    }
    
    /*! \fn friend std::ostream& operator << (std::ostream& stream, const SmallVec<T, D>& vec)
     *  \brief overloading ostream operator, keep the format settings (outside this function) from being destroyed by the non-numeric characters output */
    friend std::ostream& operator << (std::ostream& stream, const SmallVec<T, D>& vec) {
        std::streamsize tmp_width = stream.width(); // this function return std::streamsize type, which depends on machine
        std::streamsize tmp_precision = stream.precision();
        char tmp_fill = stream.fill();
        std::ios::fmtflags tmp_flags = stream.flags();  // format flags like "scientific" and "left" and "showpoint"
        
        stream << std::setw(1);
        stream << "(";
        for (int i = 0; i < D-1; i++) {
            stream.flags(tmp_flags);
            stream << std::setfill(tmp_fill) << std::setprecision(static_cast<int>(tmp_precision)) << std::setw(static_cast<int>(tmp_width));
            stream << vec.data[i];
            stream << ", ";
        }
        stream.flags(tmp_flags);
        stream << std::setfill(tmp_fill) << std::setprecision(static_cast<int>(tmp_precision)) << std::setw(static_cast<int>(tmp_width));
        stream << vec.data[D-1];
        stream << ")";
        return stream;
    }
    
    /*
     * Overload a binary operator can be done either as a non-member function or a member function.
     * For non-member function, all arguments are passed explicitly as the operator's operands. The lhs is the first argument and rhs is the second. To avoid the overhead of making a copy, operands are usually passed as reference. Take the addition as an example, "c = a + b" is equilavent to "c = operator+(a, b)".
     * For member function, it is called by an object of the class. Since the object can be accessed in the function, it is in fact passed implicitly to the member function, hence overloading binary operator as a member function requires 1 (only rhs) or 3 arguments, otherwise the compiler gives warning. While calling the operator, the lhs is the object that calls the operator, the rhs is the very explicit argument (say only 1 argument). Take the addition as an example, "c = a + b" is equilavent to "c = a.operator+(b);"
     * To avoid the complaint of compiler about member funtion with 2 arguments which overloads a binary operator, apply "friend" before the function declaration does the trick (e.g., operator << above). However, such an operator will have access to private part of this class.
     */
    
};

/*! \fn template <class T, class U, int D> inline SmallVec<typename PromoteNumeric<T, U>::type, D> operator + (const SmallVec<T, D>& lhs, const SmallVec<U, D>& rhs)
 *  \brief overloading binary operator + for class SmallVec */
template <class T, class U, int D>
inline SmallVec<typename PromoteNumeric<T, U>::type, D>
operator + (const SmallVec<T, D>& lhs, const SmallVec<U, D>& rhs) {
    SmallVec<typename PromoteNumeric<T, U>::type, D> tmp;
    for (int i = 0; i != D; i++) {
        tmp.data[i] = lhs.data[i] + rhs.data[i];
    }
    return tmp;
}

/*! \fn template <class T, class U, int D> inline SmallVec<typename PromoteNumeric<T, U>::type, D> operator - (const SmallVec<T, D>& lhs, const SmallVec<U, D>& rhs)
 *  \brief overloading binary operator - for class SmallVec */
template <class T, class U, int D>
inline SmallVec<typename PromoteNumeric<T, U>::type, D>
operator - (const SmallVec<T, D>& lhs, const SmallVec<U, D>& rhs) {
    SmallVec<typename PromoteNumeric<T, U>::type, D> tmp;
    for (int i = 0; i != D; i++) {
        tmp.data[i] = lhs.data[i] - rhs.data[i];
    }
    return tmp;
}

/*! \fn template <class T, class U, int D> inline SmallVec<typename PromoteNumeric<T, U>::type, D> operator * (const SmallVec<T, D>& lhs, const U rhs)
 *  \brief overloading binary operator * for class SmallVec (times scalar) */
template <class T, class U, int D>
inline SmallVec<typename PromoteNumeric<T, U>::type, D>
operator * (const SmallVec<T, D>& lhs, const U rhs) {
    SmallVec<typename PromoteNumeric<T, U>::type, D> tmp;
    for (int i = 0; i <D ; i++) {
        tmp.data[i] = lhs.data[i] * rhs;
    }
    return tmp;
}

/*! \fn template <class T, class U, int D> inline SmallVec<typename PromoteNumeric<T, U>::type, D> operator * (const T lhs, const SmallVec<U, D>& rhs)
 *  \brief overloading binary operator * for (scalar times) class SmallVec */
template <class T, class U, int D>
inline SmallVec<typename PromoteNumeric<T, U>::type, D>
operator * (const T lhs, const SmallVec<U, D>& rhs) {
    SmallVec<typename PromoteNumeric<T, U>::type, D> tmp;
    for (int i = 0; i != D; i++) {
        tmp.data[i] = lhs * rhs.data[i];
    }
    return tmp;
}

/*! \fn template <class T, class U, int D> inline SmallVec<typename PromoteNumeric<T, U>::type, D> operator / (const SmallVec<T, D>& lhs, const U rhs)
 *  \brief overloading binary operator / for class SmallVec (divided by scalar) */
template <class T, class U, int D>
inline SmallVec<typename PromoteNumeric<T, U>::type, D>
operator / (const SmallVec<T, D>& lhs, const U rhs) {
    SmallVec<typename PromoteNumeric<T, U>::type, D> tmp;
    for (int i = 0; i != D; i++) {
        tmp.data[i] = lhs.data[i] / rhs;
    }
    return tmp;
}

/********************************************************/
/********** Definitions of Numerical Varialbes **********/
/********************************************************/

/*! \var constexpr int dim {3}
 *  \brief dimension of simulation */
constexpr int dim {3};

/*! \class NumericalParameters
 *  \brief default numerical parameters, can be read from file */
class NumericalParameters {
private:
    
public:    
    /*! \var SmallVec<double, dim> box_center {SmallVec<double, dim>(0.0)}
     *  \brief default center position of the box */
    SmallVec<double, dim> box_center {SmallVec<double, dim>(0.0)};
    
    /*! \var SmallVec<double, dim> box_length {SmallVec<double, dim>(0.2)}
     *  \brief default side length of the box */
    SmallVec<double, dim> box_length {SmallVec<double, dim>(0.2)};
    
    /*! \var SmallVec<double, dim> box_min {SmallVec<double, dim>(-0.1)}
     *  \brief default minimum coordinates for the box */
    SmallVec<double, dim> box_min {SmallVec<double, dim>(-0.1)};
    
    /*! \var SmallVec<double, dim> box_max {SmallVec<double, dim>(0.1)}
     *  \brief default maximum coordinates for the box */
    SmallVec<double, dim> box_max {SmallVec<double, dim>(0.1)};
    
    /*! \var SmallVec<double, dim> box_half_width {SmallVec<double, dim>(0.1)}
     *  \brief half the default box side length */
    SmallVec<double, dim> box_half_width {SmallVec<double, dim>(0.1)};
    
    /*! \var double max_half_width
     *  \brief maximum half width of box
     *  This is trick for building tree for non-cubic box. But later we can implement tree that use exact half width */
    double max_half_width {0.1};
    
    /*! \var double ghost_zone_ratio = 0.025
     *  \brief width of the ghost zone
     *  Since planetesimal's size usually < 0.025H, so we adopt this figure. */
    double ghost_zone_width {0.025};
    
    /*! \var double q
     *  \brief shearing */
    double q {1.5};
    
    /*! \var double Omega
     *  \brief time unit */
    double Omega {1.0};
    
    /*! \var double shear_speed {q * Omega * box_length[0]}
     *  \brief q * Omega * Lx, shear distance per unit time */
    double shear_speed {q * Omega * box_length[0]};
    
    /*! \fn NumericalParameters()
     *  \brief constructor */
    NumericalParameters();
    
    /*! \fn void CalculateNewParameters()
     *  \brief calculate new parameters based on box shape */
    void CalculateNewParameters();
    
    /*! \fn void ReadNumericalParameters(std::string filename)
     *  \brief read numerical parameters from file */
    void ReadNumericalParameters(std::string filename);
    
};

/**********************************************/
/********** Basic_IO_Operations Part **********/
/**********************************************/

/*! \class IO_FileName
 *  \brief contains all I/O-related file path & names
 */
class IO_FileName {
private:
    
public:
    /*! \var std::string data_file_dir
     *  \brief path of directory that contains data files */
    std::string data_file_dir;
    
    /*! \var std::string data_file_basename
     *  \brief base name for data files */
    std::string data_file_basename;
    
    /*! \var std::string data_file_postname
     *  \brief post name for data files */
    std::string data_file_postname;
    
    /*! \var std::string output_file_path
     *  \brief output file for basic_analyses result */
    std::string output_file_path;
    
    /*! \var std::string input_const_path
     *  \brief input file for constant data */
    std::string input_const_path;
    
    /*! \var std::vector<std::string> lis_data_file_name
     *  \brief file names construced (for particle's lis files) */
    std::vector<std::string> lis_data_file_name;
    
#ifdef SMR_ON
    /*! \var std::string data_level
     *  \brief set level for Static-Mesh-Refinement Run */
    std::string data_level;
    
    /*! \var std::string data_domain
     *  \brief set level for Static-Mesh-Refinement Run */
    std::string data_domain;
#endif // SMR_ON
    
};

/*! \class PhysicalQuantities
 *  \brief all physical quantities that might needed to be calculated */
class PhysicalQuantities {
public:
    /*! \var float time;
     *  \brief simulation time */
    float time;
    
    /*! \var float dt;
     *  \brief simulation time step */
    float dt;
    
    /*! \var std::vector<double> particle_scale_height;
     *  \brief particle scale height for all particle sizes */
    std::vector<double> particle_scale_height;
    
    /*! \var double max_particle_density
     *  \brief maximum particle density: $\rho_p$ */
    float max_particle_density {0.0};
};

/*! \enum OutputLevel
 *  \brief served as output_level in Output() */
enum OutputLevel {
    // put two underscore at first to avoid possbile name conflict
    __normal_output = 0,        /*!< only basic output */
    __more_output,              /*!< print more info for debugging */
    __even_more_output          /*!< give all possible info */
};

/*! \enum MPI_Level
 *  \brief served as mpi_level in Output() */
enum MPI_Level {
    // put two underscore at first to avoid possbile name conflict
    __master_only = 0,          /*!< only master processor prints */
    __all_processors            /*!< all processors speak */
};

/*! \class IO_Flags
 *  \brief all possible flags used in executation */
class IO_Flags {
public:
    /*! \var int debug_flag
     *  \brief set this flag to make log_level = 2 */
    int debug_flag {0};
    
    /*! \var int verbose_flag
     *  \brief set this flag to make log_level = 1 */
    int verbose_flag {0};
    
    /*! \var int combined_flag
     *  \brief set this flag to read combined data */
    int combined_flag {0};
    
    /*! \var int find_clumps_flag
     *  \brief set this flag to find clumps & planetesimals */
    int find_clumps_flag {0};
    
    /*! \var int basic_analyses_flag
     *  \brief set this flag to perform basic data analysis */
    int basic_analyses_flag {0};
    
    /*! \var int help_flag
     *  \brief set this flag to print out usage information */
    int help_flag {0};
};

/*! \class Basic_IO_Operations
 *  \brief Handle command line and most of I/O functions together with MPI_Wrapper
 *  PS: reading method for lis files is in ParticleSet class */
class Basic_IO_Operations {
private:
    /*! \var int ostream_level
     *  \brief set level for different amount of output */
    int ostream_level {__normal_output};
    
public:
    /*! \var NumericalParameters numerical_parameters
     *  \brief default numerical parameters, can be read from file */
    NumericalParameters numerical_parameters;
    
    /*! \var IO_FileName file_name
     *  \brief all I/O-related file path & names */
    IO_FileName file_name;
    
    /*! \var std::vector<PhysicalQuantities> physical_quantities;
     *  \brief all physical quantities that might needed to be calculated */
    std::vector<PhysicalQuantities> physical_quantities;
    
    /*! \var IO_Flags flags
     *  \brief all possible flags used in executation */
    IO_Flags flags;
    
    /*! \var std::ostringstream out_content
     *  \brief this stores the content to output */
    std::ostringstream out_content;
    
    /*! \var std::ostringostream log_info
     *  \brief this stores the lof information */
    std::ostringstream log_info;
    
    /*! \var std::ostringstream error_message
     *  \brief this stores any temporary error message */
    std::ostringstream error_message;
    
    /*! \var int start_num, end_num, interval
     *  \brief start/end number and interval fo the entir file loop */
    int start_num, end_num, interval;
    
    /*! \var int num_file
     *  \brief the number of files in total */
    int num_file;
    
    /*! \var int num_cpu
     *  \brief the number of processors used in simulation */
    int num_cpu;
    
    /*! \var int width {15}
     *  \brief set default width of one data unit during output */
    int width {15};
    
    /*! \var int column {4}
     *  \brief data column in basic_analyses result file */
    int column {2};
    
    /*! \fn Basic_IO_Operations()
     *  \brief constructor */
    Basic_IO_Operations();
    
    /*! \fn int Initialization(int argc, const char * argv[])
     *  \brief initialization */
    int Initialize(int argc, const char * argv[]);
    
    /*! \fn void Output(std::ostream &stream, std::ostringstream &content, const OutputLevel &output_level, const MPI_Level &mpi_level)
     *  \brief handle output by log level & MPI status */
    void Output(std::ostream &stream, std::ostringstream &content, const OutputLevel &output_level, const MPI_Level &mpi_level);
    
    /*! \fn int PrintUsage(const char *program_name)
     *  \brief print usage */
    void PrintUsage(const char *program_name);
    
    /*! \fn void PrintStars(std::ostream &stream, const OutputLevel &output_level)
     *  \brief print 80 * symbols as a divider line */
    void PrintStars(std::ostream &stream, const OutputLevel &output_level);
    
    /*! \fn void GenerateFilenames()
     *  \brief generate the name of data files for processing */
    void GenerateFilenames();
    
    /*! \fn inline void Reset(std::ostringstream &content)
     *  \brief reset the content of ostringstream */
    inline void Reset(std::ostringstream &content) {
        content.str(std::string());
        content.clear();
        // The most elegant way to clear state and empty content is below, but it sometimes needs c++14 since different compilers have different implementations
        //std::ostringstream().swap(content);
    }
    
    /*! \fn ~Basic_IO_Operations()
     *  \brief destructor */
    ~Basic_IO_Operations();
};

/*! \var FileOperation *progIO
 *  \brief Handle command line and most of I/O functions together with MPI_Wrapper */
extern Basic_IO_Operations *progIO;

/**********************************/
/********** Utility Part **********/
/**********************************/

/*! \class MPI_Wrapper
 *  \brief served as wrappers of MPI routines plus related variables
 *  Note that even if running as a serial program, this still works since MPI_Wrapper takes care of it. The goal of this encapsulation is to reduce the use of (annoying) "#ifdef MPI_ON" in main function. */
class MPI_Wrapper {
private:
    
public:
    /*! \var int num_proc
     *  \brief number of processors */
    int num_proc;
    
    /*! \var int rank, master
     *  \brief rank of this cpu / master cpu */
    int myrank, master;
    
    /*! \var int loop_begin, loop_end, loop_step
     *  \brief begin/end/step of the file loop handled by this cpu */
    int loop_begin, loop_end, loop_step;
    
#ifdef MPI_ON
    /*! \var MPI_Comm world
     *  \brief a wrapper of MPI_COMM_WORLD */
    MPI_Comm world;
    
    /*! \var MPI_Status status
     *  \brief MPI status */
    MPI_Status status;
    
    /*! \var MPI_Offset offset
     *  \brief stream offset used during parallel file I/O */
    MPI_Offset offset {0};
    
    /*! \var MPI_Offset header_offset
     *  \brief stream offset used during parallel file I/O */
    MPI_Offset header_offset {0};
    
    /*! \var using file_obj = MPI_File
     *  \brief define a type for opening files */
    using file_obj = MPI_File;
    
#else // MPI_ON
    /*! \var using file_obj = std::ofstream
     *  \brief define a type for opening files */
    using file_obj = std::ofstream;
    
#endif // MPI_ON
    
    /*! \var file_obj result_file
     *  \brief default file to output result */
    file_obj result_file;
    
    /*! \fn MPI_Wrapper()
     *  \brief constructor */
    MPI_Wrapper();
    
    /*! \fn void Initialization(int argc, const char * argv[])
     *  \brief MPI initializaion */
    void Initialization(int argc, const char * argv[]);
    
    /*! \fn void Determine_Loop(int num_file)
     *  \brief determine the begin/end/step for file loop */
    void DetermineLoop(int num_file);
    
    /*! \fn int Barrier()
     *  \brief a wrapper of MPI_Barrier */
    void Barrier();
    
    /*! \fn int Finalize()
     *  \brief a wrapper of MPI_Finalize */
    void Finalize();
    
    /*! \fn std::string RankInfo()
     *  \brief return a string contains "Processor myrank: " */
    std::string RankInfo();
    
    /*! \fn void OpenFile(file_obj &__file, std::string filename)
     *  \brief open file */
    void OpenFile(file_obj &__file, std::string file_name);
    
    /*! \fn void WriteSingleFile(file_obj &__file, std::ostringstream &content)
     *  \brief all processor write into a file, use with cautions -> read the assumptions in descriptions
     *  When MPI_ON is on, this function assumes that you only write file header if specifying __master_only, and it assumes that every processor are writing the same amout of chunk into the file every time */
    void WriteSingleFile(file_obj &__file, std::ostringstream &content, const MPI_Level &mpi_level);
    
    /*! \fn void WriteItsOwnFile(file_obj &__file, std::ostringstream &content)
     *  \brief all processor write to its own file */
    void WriteItsOwnFile(file_obj &__file, std::ostringstream &content);
    
    /*! \fn void CloseFile(file_obj &__file)
     *  \brief close the file */
    void CloseFile(file_obj &__file);
    
    /*! \fn ~MPI_Wrapper()
     *  \brief destructor */
    ~MPI_Wrapper();
    
};

/*! \var extern MPI_Wrapper *mpi
 *  \brief declaration of gloal MPI wrapper */
extern MPI_Wrapper *mpi;

/*! \class Timer
 *  \brief served as timer */
class Timer {
private:
    /*! \var double begin_time
     *  \brief the time when this class is initiated */
    double begin_time;
    
    /*! \var bool timer_on_flag
     *  \biref a flag indicate if this timer is on or off */
    bool timer_on_flag;
    
    /*! \var double skip_time
     *  \brief if timer is resumed, [stop->resume] is skipped */
    double skip_time;
    
    /*! \fn double GetCurrentTime()
     *  \brief get current time */
    double GetCurrentTime();
    
    /*! \var double end_time
     *  \brief the time when we stop this timer */
    double stop_time;
    
public:
    /*! \var std::vector<double> lap_time;
     *  \brief lap time, just like stopwatch */
    std::vector<double> lap_time;
    
    /*! \fn Timer()
     *  \brief constructor */
    Timer();
    
    /*! \fn void StartTimer()
     *  \brief start the timer */
    void StartTimer();
    
    /*! \fn int void()
     *  \brief record a lap time and return index */
    int Lap();
    
    /*! \fn void StopTimer()
     *  \brief stop the timer */
    void StopTimer();
    
    /*! \fn void ClearTimer()
     *  \brief reset the timer */
    void ClearTimer();
    
    /*! \fn void ResumeTimer()
     *  \brief resume the timer */
    void ResumeTimer();
    
    /*! \fn double GiveTime()
     *  \brief give the time (on) or [start->stop] (off) */
    double GiveTime();
    
    /*! \fn ~Timer()
     *  \brief destructor */
    ~Timer();
};

/*! \var extern std::vector<Timer> timer
 *  \brief declaration of global timer */
extern std::vector<Timer> timer;

/*! \enum TimerTypeIndex
 *  \brief served as index of the timer */
enum TimeTypeIndex {
    // put two underscore at first to avoid possbile name conflict
    __total_elapse_time = 0,      /*!< total elapsed time */
    __waiting_time,               /*!< time for waiting, used when MPI is on */
    __tmp_used_timer,             /*!< calculate time for temporary purpose */
    __time_type_count             /*!< the number of time type */
};

/*! \fn template <typename T> T MaxOf(const T &a, const T &b)
 *  \brief return the larger one of a and b */
template <typename T>
T MaxOf(const T &a, const T &b) {
    return std::max(a, b);
}

/*! \fn template <typename T, typename... Args> T MaxOf(const T &a, const T &b, Args... args)
 *  \brief return the maximum one in a list, downgrade to MaxOf(a, b) */
template <typename T, typename... Args>
T MaxOf(const T &a, const T &b, Args... args) {
    return MaxOf(std::max(a, b), args...);
}

/*! \fn template <typename T> T MinOf(const T &a, const T &b)
 *  \brief return the smaller one of a and b */
template <typename T>
T MinOf(const T &a, const T &b) {
    return std::min(a, b);
}

/*! \fn template <typename T, typename... Args> T MinOf(const T &a, const T &b, Args... args)
 *  \brief return the minimum one in a list, downgrade to MaxOf(a, b) */
template <typename T, typename... Args>
T MinOf(const T &a, const T &b, Args... args) {
    return MaxOf(std::min(a, b), args...);
}

#endif /* global_hpp */

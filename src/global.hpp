//
//  global.hpp
//  PLAN: PLanetesimal ANalyzer
//
//  Created by Rixin Li on 4/26/16.
//  Copyright © 2016 Rixin Li. All rights reserved.
//

/*! \file global.hpp
 *  \brief provide library headers, I/O-related class and utilities */

#ifndef global_hpp
#define global_hpp

// Include C libraries first
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
// Include C++ libraries
#include <algorithm>
#include <array>
#include <bitset>
#include <chrono> // todo: use this instead of ctime in Timer class
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <functional>
// Include other libraries
#include "boost/endian/conversion.hpp"  // Boost Endian library
#include "boost/multi_array.hpp"        // Boost MultiArray library
#include "boost/dynamic_bitset.hpp"     // Boost Dynamic Bitset library
#include "boost/algorithm/string/trim.hpp"
#include <getopt.h>
#include <unistd.h>

//#define MPI_ON // Comment out this line before committing!!
#ifdef MPI_ON // "ifdef" options are defined during compilation
#include "mpi.h"
#endif // MPI_ON

// Check c++11 support
#if __cplusplus <= 199711L
#error This program wants a C++11 compliant compiler (option -DOLDCPP which supports old compilers has been abandoned).
#endif // __cplusplus

// Disable assert() for the production version (we may want to merge Debug flag with this)
//#define NDEBUG

/***********************************/
/********** SmallVec Part **********/
/***********************************/

/*
 * N.B.: template classes need to have the method definitions inside the header file, in order to let the linker work. Alternatively, you can put definitions in other file and include after declaration in the header file.
 */

/*! \class template <bool, class T, class U> struct __SelectIf_base
 *  \brief a template for struct __SelectIF to accept a bool as condition */
template <bool, class T, class U>
struct __SelectIf {};

/*! \class template <class T, class U> struct __SelectIf_true
 *  \brief a template for struct __SelectIF to select T type if bool is true */
template <class T, class U>
struct __SelectIf<true, T, U> { typedef T type; };

/*! \class template <class T, class U> struct __SelectIf_false
 *  \brief a template for struct __SelectIF to select U type if bool is false */
template <class T, class U>
struct __SelectIf<false, T, U> { typedef U type; };

/*! \class template <class T, class U> struct PromoteNumeric
 *  \brief this struct controls the type promotion, this struct nests many levels of selections. Read comments/explanations from inside and notice that comments start with new lines are different with comments after statements */
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

    ////////////////////////////////////
    /////////// Constructors ///////////
    ////////////////////////////////////
    
    /*! \fn SmallVec() : data{0}
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
     * If explicit keyword is applied to a copy constructor, it means that object of that class can't be copied when being passed to functions or when being returned from function - (this type of copying is called implicit copying). So statement "SmallVec<T, D> new_vec = old_vec;" will cause error during compilation. Passing SmallVec as an function argument by value or make it the return type of a function will also cause error, like "Func(old_vec)" or "SmallVec<T, D> Func()". Only explicit copying, i.e., "SmallVec<T, D> new_vec (old_vec);" is allowed. And only being passed to a function (or being returned from a functionn) by reference/pointer is allowed.
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
     * C++ has a special parameter type, ellipsis, that can be used to pass a varying number of arguments. It is called the ellipsis operator. Try to understand it by recalling the usage of printf() in C where you can input any amount of arguments. Also, sizeof...() is used to count the number of arguments, which is different with sizeof()(reference: http://www.cplusplus.com/articles/EhvU7k9E/ )
     * Then "data{ head, T(tail)... }" is called an initialization list. Usually it is for passing arguments to the constructor of a parent class.
     */

    /////////////////////////////////
    /////////// Operators ///////////
    /////////////////////////////////

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
     *  \brief compound assignment operator *= with scalar */
    template <class U>
    SmallVec<T, D>& operator *= (const U rhs) {
        for (int i = 0; i != D; i++) {
            data[i] *= rhs;
        }
        return *this;
    }

    /*! \fn template <class U> SmallVec<T, D>& operator *= (const SmallVec<U, D>& rhs)
 *  \brief compound assignment operator *= with class SmallVec */
    template <class U>
    SmallVec<T, D>& operator *= (const SmallVec<U, D>& rhs) {
        for (int i = 0; i != D; i++) {
            data[i] *= rhs[i];
        }
        return *this;
    }
    
    /*! \fn template <class U> SmallVec<T, D>& operator /= (const U rhs)
     *  \brief compound assignment operator /= with scalar */
    template <class U>
    SmallVec<T, D>& operator /= (const U rhs) {
        for (int i = 0; i != D; i++) {
            data[i] /= rhs;
        }
        return *this;
    }

    /*! \fn template <class U> SmallVec<T, D>& operator /= (const SmallVec<U, D>& rhs)
*  \brief compound assignment operator /= with class SmallVec */
    template <class U>
    SmallVec<T, D>& operator /= (const SmallVec<U, D>& rhs) {
        for (int i = 0; i != D; i++) {
            data[i] /= rhs[i];
        }
        return *this;
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
     *  \brief overloading operator <, we compare the norm of vectors. */
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
     *  \brief overloading operator >, we compare the norm of vectors. */
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

    ////////////////////////////////////////
    /////////// Simple Iterators ///////////
    ////////////////////////////////////////

    // not implemented yet

    ////////////////////////////////////////
    /////////// Member Functions ///////////
    ////////////////////////////////////////
    
    /*! \fn typename PromoteNumeric<T, double>::type Norm() const
     *  \brief calculate the norm of this vector, |x|, for at least double precision */
    typename PromoteNumeric<T, double>::type Norm() const {
        typename PromoteNumeric<T, double>::type sum = 0;
        for (int i = 0; i != D; i++) {
            sum += data[i] * data[i];
        }
        return std::sqrt(sum);
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
    
    /*! \fn template <class U> inline SmallVec<typename PromoteNumeric<T, U>::type, 3> ParaMultiply(const SmallVec<U, 3>& rhs) const
     *  \brief give a new SmallVec, where data[i] = this[i] * rhs[i]
     *  Note: the multiply operator '*' has been overloaded to do this */
    template <class U>
    inline SmallVec<typename PromoteNumeric<T, U>::type, 3> ParaMultiply(const SmallVec<U, 3>& rhs) const {
        SmallVec<typename PromoteNumeric<T, U>::type, 3> tmp;
        tmp[0] = data[0] * rhs[0];
        tmp[1] = data[1] * rhs[1];
        tmp[2] = data[2] * rhs[2];
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
            val = val && std::fabs(diff[i] < epsilon);
        }
        return val;
    }
    
    /*! \fn template <class U, class V> int relclose(const SmallVec<U, D>& rhs, const V epsilon) const
     *  \brief determine (the relative difference between this_vector and rhs) divided by (the average length of them) is less than epsilon in each element */
    template <class U, class V>
    int RelClose(const SmallVec<U, D>& rhs, const V epsilon) const {
        SmallVec<typename PromoteNumeric<T, U>::type, D> sum, diff;
        for (int i = 0; i != D; i++) {
            sum[i] = std::fabs(data[i]) + std::fabs(rhs[i]);
        }
        diff = *this - rhs;
        bool val = true;
        for (int i = 0; i != D; i++) {
            val = val && ((2 * std::fabs(diff[i]) / sum[i]) < epsilon);
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
    
    /*! \fn SmallVec<T, D> SetZeros()
     *  \brief reset data to {0} */
    SmallVec<T, D> SetZeros() {
        for (int i = 0; i != D; i++) {
            data[i] = 0;
        }
        return *this;
    }

    /*! \fn SmallVec<T, D> AbsSelf()
     *  \brief use abs() to each element */
    SmallVec<T, D> AbsSelf() {
        for (int i = 0; i != D; i++) {
            data[i] = std::abs(data[i]);
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
    

    /*
     * Overload a binary operator can be done either as a non-member function or a member function.
     * For non-member function, all arguments are passed explicitly as the operator's operands. The lhs is the first argument and rhs is the second. To avoid the overhead of making a copy, operands are usually passed as reference. Take the addition as an example, "c = a + b" is equivalent to "c = operator+(a, b)".
     * For member function, it is called by an object of the class. Since the object can be accessed in the function, it is in fact passed implicitly to the member function, hence overloading binary operator as a member function requires 1 (only rhs) or 3 arguments, otherwise the compiler gives warning. While calling the operator, the lhs is the object that calls the operator, the rhs is the very explicit argument (say only 1 argument). Take the addition as an example, "c = a + b" is equivalent to "c = a.operator+(b);"
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

/*! \fn template <class T, class U, int D> inline SmallVec<typename PromoteNumeric<T, U>::type, D> operator * (const SmallVec<T, D>& lhs, const SmallVec<U, D>& rhs)
 *  \brief overloading binary operator * for class SmallVec times class SmallVec */
template <class T, class U, int D>
inline SmallVec<typename PromoteNumeric<T, U>::type, D>
operator * (const SmallVec<T, D>& lhs, const SmallVec<U, D>& rhs) {
    SmallVec<typename PromoteNumeric<T, U>::type, D> tmp;
    for (int i = 0; i != D; i++) {
        tmp.data[i] = lhs.data[i] * rhs.data[i];
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

/*! \fn template <class T, class U, int D> inline SmallVec<typename PromoteNumeric<T, U>::type, D> operator / (const SmallVec<T, D>& lhs, const SmallVec<U, D>& rhs)
 *  \brief overloading binary operator / for class SmallVec divided by class SmallVec */
template <class T, class U, int D>
inline SmallVec<typename PromoteNumeric<T, U>::type, D>
operator / (const SmallVec<T, D>& lhs, const SmallVec<U, D>& rhs) {
    SmallVec<typename PromoteNumeric<T, U>::type, D> tmp;
    for (int i = 0; i != D; i++) {
        tmp.data[i] = lhs.data[i] / rhs.data[i];
    }
    return tmp;
}

template<class T, class U, int D>
bool SmallVecLessEq (const SmallVec<T, D>& lhs, const SmallVec<U, D>& rhs) {
    bool val = true;
    for (int i = 0; i != D; i++) {
        val = val && (lhs[i] <= rhs[i]);
    }
    return val;
};

template<class T, class U, int D>
bool SmallVecGreatEq (const SmallVec<T, D>& lhs, const SmallVec<U, D>& rhs) {
    bool val = true;
    for (int i = 0; i != D; i++) {
        val = val && (lhs[i] >= rhs[i]);
    }
    return val;
};


/********************************************************/
/********** Definitions of Numerical Varialbes **********/
/********************************************************/

/*! \var constexpr int dim {3}
 *  \brief dimension of simulation */
constexpr int dim {3};

/*! \class NumericalParameters
 *  \brief numerical parameters used in simulations, can be read from file */
class NumericalParameters {
private:
    
public:

    /***** The following parameters can be retrieved from *.lis files *****/

    /*! \var SmallVec<double, dim> box_center {SmallVec<double, dim>(0.0)}
     *  \brief center position of the box, default (0, 0, 0) */
    SmallVec<double, dim> box_center {SmallVec<double, dim>(0.0)};

    /*! \var SmallVec<double, dim> box_length {SmallVec<double, dim>(0.2)}
     *  \brief side length of the box, default (0.2, 0.2, 0.2) */
    SmallVec<double, dim> box_length {SmallVec<double, dim>(0.2)};

    /*! \var SmallVec<double, dim> box_min {SmallVec<double, dim>(-0.1)}
     *  \brief minimum coordinates for the box, default (-0.1, -0.1, -0.1) */
    SmallVec<double, dim> box_min {SmallVec<double, dim>(-0.1)};

    /*! \var SmallVec<double, dim> box_max {SmallVec<double, dim>(0.1)}
     *  \brief maximum coordinates for the box, default (0.1, 0.1, 0.1) */
    SmallVec<double, dim> box_max {SmallVec<double, dim>(0.1)};

    /*! \var SmallVec<double, dim> box_half_width {SmallVec<double, dim>(0.1)}
     *  \brief half the box side length, default (0.1, 0.1, 0.1) */
    SmallVec<double, dim> box_half_width {SmallVec<double, dim>(0.1)};

    /*! \var double max_half_width
     *  \brief maximum half width of box
     *  This is trick for building tree for non-cubic box. But later we can implement tree that use exact half width */
    double max_half_width {0.1};

    /***** The following parameters need input information and/or derivation *****/

    /*! \var SmallVec<unsigned int, dim> Nx {SmallVec<unsigned int, dim>(128)}
     *  \brief box resolution, default (128, 128, 128) */
    SmallVec<unsigned int, dim> box_resolution {SmallVec<unsigned int, dim>(128)};
    
    /*! \var SmallVec<double, dim> cell_length {SmallVec<double, dim>(0.0015625)}
     *  \brief cell side length, default (0.0015625, 0.0015625, 0.0015625) */
    SmallVec<double, dim> cell_length {SmallVec<double, dim>(0.0015625)};

    /*! \var double cell_volume
     *  \brief cell volume, default 0.0015625^3 */
    double cell_volume {0.0015625*0.0015625*0.0015625};

    /*! \var double ghost_zone_width {0.025}
     *  \brief width of the ghost zone
     *  Since planetesimal's size usually < 0.025H, we adopt this number. A smaller number may be chosen. */
    double ghost_zone_width {0.025};
    
    /*! \var double q
     *  \brief shearing */
    double q {1.5};
    
    /*! \var double Omega
     *  \brief inverse time unit */
    double Omega {1.0};

    /*! \var double Omega_squared
     *  \brief Omega^2 */
    double Omega_squared {Omega*Omega};

    /*! \var double rho_g0
     *  \brief density unit, midplane gas density */
    double rho_g0 {1.0};

    /*! \var double solid_to_gas_ratio
     *  \brief solid-to-gas ratio, usually 0.02 */
    double solid_to_gas_ratio {0.02};

    /*! \var std::vector<double> mass_fraction_per_species
     *  \brief the fraction of total mass that each particle species has */
    std::vector<double> mass_fraction_per_species {std::vector<double>(1, 1)};

    /*! \var std::vector<double> mass_per_particle
     *  \brief mass of one particle for each type */
    std::vector<double> mass_per_particle; // {std::vector<double>(1, 03.82481121006927567e-9)};

    /*! \var double mass_total_code_units
     *  \brief total particle mass in code unit */
    double mass_total_code_units {0.0020053};

    /*! \var const double PI
     *  \brief the constant PI */
    const double PI = 3.14159265358979323846;

    /*! \var const double four_PI
     *  \brief constant 4 * PI */
    const double four_PI = 4. * PI;

    /*! \var const double four_PI_over_three
     *  \brief constant 4 * PI / 3 */
    const double four_PI_over_three = PI * 4. / 3.;

    /*! \var double four_PI_G
     *  \brief the numerical value provided to Athena */
    double four_PI_G {0.05};

    /*! \var double G_tilde
     *  \brief the code unit indicating the strength of self-gravity */
    double G_tilde {four_PI_G * rho_g0 / Omega / Omega};

    /*! \var const double mass_ceres
     *  \brief the mass of asteroid Ceres (in kg) */
    const double mass_ceres = 9.4e20;

    /*! \var double mass_physical
     *  \brief the physical mass in the box assuming 3 AU and MMSN (in kg) */
    double mass_physical {(G_tilde/0.1)*1.3e24}; // For G=0.1, r=3AU, M_0 = rho_0*H^3 = 1.3e24

    /*! \var double grav_constant
     *  \brief the gravitational constant in code unit */
    double grav_constant {four_PI_G / four_PI};

    /*! \var double shear_speed {q * Omega * box_length[0]}
     *  \brief q * Omega * Lx, shear distance per unit time */
    double shear_speed {q * Omega * box_length[0]};

    /*! \var unsigned int num_neighbors_in_knn_search
     *  \brief determine how many (K) Nearest Neighbors that we need to search for density calculations */
    unsigned int num_neighbors_in_knn_search {64};

    /*! \var unsigned int num_neighbors_to_hop
     *  \brief how many neighbors to search while looking for the densest neighbor */
    unsigned int num_neighbors_to_hop {32};

    /*! \var unsigned int num_peaks
     *  \brief the maximum number of clumps/peaks to OUTPUT */
    unsigned int num_peaks {0};

    /*! \var double min_trusted_mass_code_unit
     *  \brief the minimum mass of a clump/peak to OUTPUT, determined by R_Hill = cell_length
     *  The Hill Radius is defined as R_Hill^3 = G Mass / (3 * Omega^2), where G = four_PI_G / four_PI.
     *  Substitute dx into R_Hill, we have: Mass_min / (rho_g0 H^3) = 3 (dx/H)^3 (Omega^2/(rho_g0 G)).
     *  For the fiducial run with 128^3 cells in (0.2H)^3 box, Mass_min = 2.87621e-6, i.e., 0.1434% of total mass.
     *  Such a Mass_min corresponds to ~24 particles for the subsample of 16384 single-species particles. */
    double min_trusted_mass_code_unit {2.87621e-6};

    /*! \var std::unordered_map<std::string, double> input_paras
     *  \brief take input and serve as dictionary for single value parameters */
    std::unordered_map<std::string, double> input_paras;

    /*! \fn NumericalParameters()
     *  \brief constructor */
    NumericalParameters();
    
    /*! \fn void CalculateNewParameters()
     *  \brief calculate new parameters based on box shape */
    void CalculateNewParameters();
    
    /*! \fn void ReadNumericalParameters(std::string filename)
     *  \brief read numerical parameters from Athena input file */
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
    
    /*! \var std::string max_rhop_vs_scale_file_path
     *  \brief output file for max_rhop at all scales */
    std::string max_rhop_vs_scale_file;

    /*! \var std::string mean_sigma_file
     *  \brief output file for <Sigma_g>_{yz}(t) and <Sigma_p>_{yz}(t) */
    std::string mean_sigma_file;

    /*! \var std::string input_const_path
     *  \brief input file for constant data */
    std::string input_const_path;
    
    /*! \var std::vector<std::string> lis_data_file_name
     *  \brief file names (for particle's lis files) */
    std::vector<std::string> lis_data_file_name;
    
    /*! \var std::vector<std::string> vtk_data_file_name
     *  \brief file names (for vtk files) */
    std::vector<std::string> vtk_data_file_name;

    /*! \var std::string planetesimals_file
     *  \brief output file for basic planetesimal-info */
    std::string planetesimals_file;
    
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
    /*! \var double time;
     *  \brief simulation time */
    double time;
    
    /*! \var double dt;
     *  \brief simulation time step */
    double dt;
    
    /*! \var std::vector<double> particle_scale_height;
     *  \brief particle scale height for all particle sizes */
    std::vector<double> particle_scale_height;
    
    /*! \var double max_particle_density
     *  \brief maximum particle density: $\rho_p$ */
    double max_particle_density {0.0};
    
    /*! \var std::vector<double> max_rhop_vs_scale
     *  \brief maximum sphere-based particle density at all length scales */
    std::vector<double> max_rhop_vs_scale;

    /*! \var double vertical_flux
     *  \brief vertical gas mass flux on vertical boundaries (per unit area) */
    double vertical_flux {0.0};

    /*! \var std::vector<double> mean_sigma
     *  \brief <Sigma_g>_{yz}(t) and <Sigma_p>_{yz}(t) */
    std::vector<double> mean_sigma;

    /*! \var double max_planetesimal_mass
     *  \brief the maximum mass in the planetesimals */
    double max_planetesimal_mass {0};

    /*! \var double total_planetesimal_mass
     *  \brief total mass in planetesimals */
    double total_planetesimal_mass {0};
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
 *  \brief all possible flags used in execution */
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
    
    /*! \var int density_vs_scale_flag
     *  \brief set this flag to calculate max(rho_p) at all scales */
    int density_vs_scale_flag {0};
    
    /*! \var int no_ghost_particle_flag
     *  \brief set this flag to prevent making ghost particles */
    int no_ghost_particle_flag {0};
    
    /*! \var int tmp_calculation_flag
     *  \brief set this flag to do some temporary calculations */
    int tmp_calculation_flag {0};
    
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
     *  \brief all possible flags used in execution */
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
    
    /*! \var int num_files
     *  \brief the number of files in total */
    int num_files;
    
    /*! \var int num_cpus
     *  \brief the number of processors used in simulation */
    int num_cpus;
    
    /*! \var int width {15}
     *  \brief set default width of one data unit during output */
    int width {15};
    
    /*! \var int column {4}
     *  \brief data column in basic_analyses result file
     *  Here we start with 2, time and max(rho_p) */
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

    /*! \fn std::string LocalTime()
     *  \brief return a string containing the date and time information */
    std::string LocalTime();
    
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
    /*! \var int num_processors
     *  \brief number of processors */
    int num_processors;
    
    /*! \var int rank, master
     *  \brief rank of this cpu / master cpu */
    int myrank, master;
    
    /*! \var int loop_begin, loop_end, loop_step
     *  \brief begin/end/step of the file loop handled by this cpu */
    int loop_begin, loop_end, loop_step;
    
#ifndef MPI_ON

#if defined(__GNUC__) && (__GNUC__ < 5) && (!__clang__)
    /* RL: std::ofstream should be move constructable but there is no move constructor declared before GCC 5.1
     * Thus, we use pointers to compromise this bug here.
     * refer to https://stackoverflow.com/questions/28775673/why-cant-i-move-stdofstream
     *      and https://gcc.gnu.org/bugzilla/show_bug.cgi?id=54316
     *
     * For Clang users, it should be fine if
     *     Apple LLVM version > 6.0 or LLVM version > 3.5
     * Since Mac OS uses clang as default, __clang__ = true(1)
     */

    /*! \alias using file_obj = std::ofstream *
     *  \brief define a type for opening files */
    using file_obj = std::ofstream *;
#else
    /*! \alias using file_obj = std::ofstream
     *  \brief define a type for opening files */
    using file_obj = std::ofstream;
#endif

#else // MPI_ON
    /*! \var MPI_Comm world
     *  \brief a wrapper of MPI_COMM_WORLD */
    MPI_Comm world;

    /*! \var MPI_Status status
     *  \brief MPI status */
    MPI_Status status;

    /*! \alias using file_obj = MPI_File
     *  \brief define a type for opening files */
    using file_obj = MPI_File;

    /*! \var MPI_Comm file_writing_communicator
     *  \brief MPI communicator used for writing files (to avoid Bcase error while num_files < num_processors) */
    MPI_Comm file_writing_communicator;

    /*! \var int file_writing_master
     *  \brief mark out the processor used to write file headers */
    int file_writing_master;

    /*! \var std::map<file_obj, MPI_Offset> offset
     *  \brief stream offset used during parallel file I/O */
    std::map<file_obj, MPI_Offset> offset;

    /*! \var std::map<file_obj, MPI_Offset> header_offset
     *  \brief stream offset used during parallel file I/O */
    std::map<file_obj, MPI_Offset> header_offset;

#endif // MPI_ON
    
    /*! \var std::vector<file_obj> result_files
     *  \brief files to output results */
    std::vector<file_obj> result_files;
    
    /*! \var std::map<std::string, std::vector<file_obj>::size_type> file_pos;
     *  \brief mapping file names into result_files vector
     *  I tried to use std::map<std::string, std::vector<file_obj>::iterator>, but ended up with segmentation fault on std::ofstream. It seems like using the iterator of file_obj is not a good idea. */
    std::map<std::string, std::vector<file_obj>::size_type> file_pos;
    
    /*! \fn MPI_Wrapper()
     *  \brief constructor */
    MPI_Wrapper();
    
    /*! \fn void Initialization(int argc, const char * argv[])
     *  \brief MPI initializaion */
    void Initialization(int argc, const char * argv[]);
    
    /*! \fn void Determine_Loop(int num_files)
     *  \brief determine the begin/end/step for file loop */
    void DetermineLoop(int num_files);
    
    /*! \fn int Barrier()
     *  \brief a wrapper of MPI_Barrier */
    void Barrier();
    
    /*! \fn int Finalize()
     *  \brief a wrapper of MPI_Finalize */
    void Finalize();
    
    /*! \fn std::string RankInfo()
     *  \brief return a string contains "Processor myrank: " */
    std::string RankInfo();
    
    /*! \fn void OpenFile(std::string filename)
     *  \brief open file */
    void OpenFile(std::string file_name);
    
    /*! \fn void WriteSingleFile(file_obj &__file, std::ostringstream &content)
     *  \brief all processor write into a file, use with cautions -> read the assumptions in descriptions
     *  When MPI_ON is on, this function assumes that you only write file header if specifying __master_only, and it assumes that every processor are writing the same amout of chunk into the file every time */
    void WriteSingleFile(file_obj &__file, std::ostringstream &content, const MPI_Level &mpi_level);
    
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
    
    /*! \fn double GiveTime(const unsigned int i)
     *  \brief overload: give the lap time */
    double GiveTime(const unsigned int i);
    
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

/*! \namespace sn
 *  \brief store simplified names */
namespace sn {
    
    /*! \alias using b_range = boost::multi_array_types::index_range
     *  \brief used to generate boost index range */
    using b_range = boost::multi_array_types::index_range;
    
    // RL: consider add more here

    using dvec = SmallVec<double, dim>;
}



#endif /* global_hpp */

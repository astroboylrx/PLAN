//
//  tree.hpp
//  PLATO: PLAneTesimal locatOr
//
//  Created by Rixin Li on 3/11/16.
//  Copyright © 2016 Rixin Li. All rights reserved.
//
//  The following classes are based on programs written by Dr. Philip Pinto during his course (ASTR596) in fall 2015

#ifndef tree_hpp
#define tree_hpp

#include "global.hpp"

/*
 * An important thing to remember before everything: template classes need to have the method definitions inside the header file.
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
 * this controls the type promotion, this struct nests many levels of selections. Read comments/explanations from inside and notice that comments start with new lines are different with comments after statements */
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
public:
    /*! \var T data[D]
     *  \brief data of elements of this vector */
    T data[D];
    
    /*! \fn SmallVec() : data{0} {}
     *  \brief the most basic constructor, list initialization with {0} */
    SmallVec() : data{0} {}
    
    /*! \fn template <class U> explicit SmallVec( const U& scalar)
     *  \brief overloading constructor to broadcast a scalar, e.g., unit_vec = SmallVec<double, 3>(1.0). For keyword "explicit", read explanations below */
    template <class U>
    explicit SmallVec( const U& scalar) {
        for(int i=0; i<D; i++) data[i] = static_cast<T>(scalar);
    }
    /*
     * Why do we use the keyword "explicit"? (reference: http://www.geeksforgeeks.org/g-fact-93/ )
     * In C++, if a class has a constructor which can be called with a single argument, then this constructor becomes conversion constructor because such a constructor allows conversion of the single argument to the class being constructed.
     * We can avoid such implicit conversions as these may lead to unexpected results. We can make the constructor explicit with the help of explicit keyword.
     */
    
    /*! \fn template <typename... Tail> SmallVec(typename std::enable_if<sizeof...(Tail)+1==D, T>::type head, Tail... tail) : data{ head, T(tail)... }
     *  \brief overloading constructor for initialization in the form of e.g., SmallVec<double, 3>(3, 4.5, 5) */
    template <typename... Tail>
    SmallVec(typename std::enable_if<sizeof...(Tail)+1==D, T>::type head, Tail... tail)
    : data{ head, T(tail)... } {}
    /*
     * C++ also has a special parameter type, ellipsis, that can be used to pass a varying number of arguments. It is called the ellipsis operator. Try to undersatnd it by recalling the usage of printf() in C where you can input any amount of arguments.
     */
    
    // allow access with usual vector notation
    T operator[] ( const size_t i ) const {
        assert( i < D );
        return *(data+i);
    }
    T& operator [] ( const size_t i ) {
        assert( i < D );
        return *(data+i);
    }
    
    //unary operations (sign)
    const SmallVec& operator +() {
        return *this;
    }
    
    SmallVec operator -() {
        SmallVec<T,D> result;
        for(int i=0; i<D; i++) result[i] = -data[i];
        return result;
    }
    
    // combinations of arithmetic functions and =
    template <class U>
    SmallVec<T,D>& operator += ( const SmallVec<U,D>& rhs ) {
        for(int i=0; i<D; i++) data[i] += rhs.data[i];
        return *this;
    }
    
    template <class U>
    SmallVec<T,D>& operator -= ( const SmallVec<U,D>& rhs ) {
        for(int i=0; i<D; i++) data[i] -= rhs.data[i];
        return *this;
    }
    
    template <class U>
    SmallVec<T,D>& operator *= ( const U rhs ) {
        for(int i=0; i<D; i++) data[i] *= rhs;
        return *this;
    }
    
    template <class U>
    SmallVec<T, D>& operator /= ( const U rhs ) {
        for(int i=0; i<D; i++) data[i] /= rhs;
        return *this;
    }
    
    // functions of components
    typename PromoteNumeric<T, double>::type norm() const {
        typename PromoteNumeric<T, double>::type sum = 0;
        for(int i=0; i<D; i++) sum += data[i]*data[i];
        return sqrt(sum);
    }
    
    typename PromoteNumeric<T, double>::type norm2() const {
        typename PromoteNumeric<T, double>::type sum = 0;
        for(int i=0; i<D; i++) sum += data[i]*data[i];
        return sum;
    }
    
    T maxcomponent() const {
        T maxc = data[0];
        for(int i=1; i<D; i++) maxc = (data[i]>maxc)?data[i]:maxc;
        return maxc;
    }
    T mincomponent() const {
        T minc = data[0];
        for(int i=1; i<D; i++) minc = (data[i]<minc)?data[i]:minc;
        return minc;
    }
    
    template<class U>
    typename PromoteNumeric<T, U>::type dot(const SmallVec<U, D>& rhs) const {
        typename PromoteNumeric<T, U>::type tmp = data[0]*rhs[0];
        for(int i=1; i<D; i++) tmp += data[i]*rhs[i];
        return tmp;
    }
    
    template <class U>
    inline SmallVec<typename PromoteNumeric<T, U>::type, 3>
    cross( const SmallVec<U, 3>& rhs ) const {
        SmallVec<typename PromoteNumeric<T, U>::type, 3> tmp;
        tmp[0] = data[1] * rhs[2] - data[2] * rhs[1];
        tmp[1] = data[2] * rhs[0] - data[0] * rhs[2];
        tmp[2] = data[0] * rhs[1] - data[1] * rhs[0];
        return tmp;
    }
    
    // abs. diff. between *this and rhs is less than epsilon in each component
    template <class U, class V>
    bool absclose( const SmallVec<U, D>& rhs, const V epsilon ) const {
        SmallVec<typename PromoteNumeric<T, U>::type, D> diff;
        diff = *this - rhs;
        bool val = true;
        for(int i=0; i<D; i++) val = val && fabs(diff[i] < epsilon);
        return val;
    }
    
    // rel. diff. between *this and rhs is less than epsilon in each component
    template <class U, class V>
    int relclose( const SmallVec<U, D>& rhs, const V epsilon ) const {
        SmallVec<typename PromoteNumeric<T, U>::type, D> sum, diff;
        for(int i=0; i<D; i++) sum[i] = fabs(data[i]) + fabs(rhs[i]);
        diff = *this - rhs;
        bool val = true;
        for(int i=0; i<D; i++) val = val && (( 2*fabs(diff[i]) / sum[i] ) < epsilon);
        return val;
    }
    
    bool is_finite() const {
        bool val = true;
        for(int i=0; i<D; i++) val = val && std::isfinite(data[i]);
        return val;
    }
    
    // relational operators
    template <class U>
    bool operator == ( const SmallVec<U, D>& rhs ) const {
        bool val = true;
        for(int i=0; i<D; i++) val = val && (data[i] == rhs[i]);
        return val;
    }
    
    template <class U>
    bool operator != ( const SmallVec<U, D>& rhs ) const {
        return !(*this == rhs);
    }
    
    // define < as in sorted order by dimension
    template <class U>
    bool operator < ( const SmallVec<U, D>& rhs ) const {
        if(data[0] < rhs[0]) return true;
        for(int i=1; i<D; i++)
            if( data[i-1]==rhs[i-1] && data[i] < rhs[i] ) return true;
        return false;
    }
    
    template <class U>
    bool operator <= ( const SmallVec<U, D>& rhs ) const {
        if(  (*this < rhs) || (*this==rhs) ) return true;
        return false;
    }
    
    // define > as in reverse sorted order by dimension
    template <class U>
    bool operator > ( const SmallVec<U, D>& rhs ) const {
        if(data[0] > rhs[0]) return true;
        for(int i=1; i<D; i++)
            if( data[i-1]==rhs[i-1] && data[i] > rhs[i] ) return true;
        return false;
    }
    
    template <class U>
    bool operator >= ( const SmallVec<U, D>& rhs ) const {
        if( (*this>rhs) || (*this==rhs) ) return true;
        return false;
    }
    
    SmallVec<T, D> zero() {
        for(int i=0; i<D; i++) data[i] = 0;
        return *this;
    }
    
    // return true if all components are on [a,b]
    template <class U, class V>
    bool inrange(U low, V hi) {
        bool val = true;
        for(int i=0; i<D; i++) val = val && (data[i] >= low && data[i] <= hi);
        return val;
    }
    
    template <class U, class V>
    bool inrange(SmallVec<U, D> low, SmallVec<V, D> hi) {
        bool val = true;
        for(int i=0; i<D; i++) val = val && (data[i] >= low[i] && data[i] <= hi[i]);
        return val;
    }
    
    // stream operator: keep the format settings from being destroyed by the
    // non-numeric characters output
    friend std::ostream& operator <<( std::ostream& o, const SmallVec<T, D>& v ) {
        std::streamsize tmpw = o.width();
        std::streamsize tmpp = o.precision();
        char tmps = o.fill();
        std::ios::fmtflags tmpf = o.flags();  // format flags like "scientific" and "left" and "showpoint"
        o << std::setw(1);
        o << "(";
        for(int i=0; i<D-1; i++) {
            o.flags(tmpf); o << std::setfill(tmps) << std::setprecision(tmpp) << std::setw(tmpw);
            o << v.data[i];
            o << ",";
        }
        o.flags(tmpf); o << std::setfill(tmps) << std::setprecision(tmpp) << std::setw(tmpw);
        o << v.data[D-1];
        o << ")";
        return o;
    }
    
};

// componentwise addition and subtraction
template <class T, class U, int D>
inline SmallVec<typename PromoteNumeric<T, U>::type, D>
operator + (const SmallVec<T, D>& lhs, const SmallVec<U, D>& rhs ) {
    SmallVec<typename PromoteNumeric<T, U>::type, D> tmp;
    for(int i=0; i<D; i++) tmp.data[i] = lhs.data[i] + rhs.data[i];
    return tmp;
}

template <class T, class U, int D>
inline SmallVec<typename PromoteNumeric<T, U>::type, D>
operator - (const SmallVec<T, D>& lhs, const SmallVec<U, D>& rhs ) {
    SmallVec<typename PromoteNumeric<T, U>::type, D> tmp;
    for(int i=0; i<D; i++) tmp.data[i] = lhs.data[i] - rhs.data[i];
    return tmp;
}

// left and right multiplication by a scalar
template <class T, class U, int D>
inline SmallVec<typename PromoteNumeric<T, U>::type, D>
operator * ( const SmallVec<T, D>& lhs, const U rhs ) {
    SmallVec<typename PromoteNumeric<T, U>::type, D> tmp;
    for(int i=0; i<D; i++) tmp.data[i] = lhs.data[i]*rhs;
    return tmp;
}

template <class T, class U, int D>
inline SmallVec<typename PromoteNumeric<T, U>::type, D>
operator * ( const T lhs, const SmallVec<U, D>& rhs ) {
    SmallVec<typename PromoteNumeric<T, U>::type, D> tmp;
    for(int i=0; i<D; i++) tmp.data[i] = lhs*rhs.data[i];
    return tmp;
}

// right division by a scalar
template <class T, class U, int D>
inline SmallVec<typename PromoteNumeric<T, U>::type, D>
operator / ( const SmallVec<T, D>& lhs, const U rhs ) {
    SmallVec<typename PromoteNumeric<T, U>::type, D> tmp;
    for(int i=0; i<D; i++) tmp.data[i] = lhs.data[i]/rhs;
    return tmp;
}














#endif /* tree_hpp */

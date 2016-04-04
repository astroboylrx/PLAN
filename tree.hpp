//
//  tree.hpp
//  PLATO: PLAneTesimal locatOr
//
//  Created by Rixin Li on 3/11/16.
//  Copyright © 2016 Rixin Li. All rights reserved.
//
//  The following code are based on programs written by Dr. Philip Pinto during his course (ASTR596) in fall 2015.
//  Descriptions and comments are written by Rixin Li.

#ifndef tree_hpp
#define tree_hpp

#include "global.hpp"

/***********************************/
/********** SmallVec Part **********/
/***********************************/

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

/*
 * N.B.: template classes need to have the method definitions inside the header file, in order to let the linker work.
 * Alternatively, you can put definitions in other file and include after declaration in the header file.
 */

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
    
    /*! \fn template <class U> explicit SmallVec( const U& scalar)
     *  \brief overloading constructor to broadcast a scalar, e.g., unit_vec = SmallVec<double, 3>(1.0). For keyword "explicit", read explanations below */
    template <class U>
    explicit SmallVec(const U& scalar) {
        for (int i = 0; i < D; i++) {
            data[i] = static_cast<T>(scalar);
        }
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
        for (int i = 0; i < D; i++) {
            result[i] = -data[i];
        }
        return result;
    }
    
    /*! \fn template <class U> SmallVec<T, D>& operator += (const SmallVec<U, D>& rhs)
     *  \brief compound assignment operator += */
    template <class U>
    SmallVec<T, D>& operator += (const SmallVec<U, D>& rhs) {
        for (int i = 0; i < D; i++) {
            data[i] += rhs.data[i];
        }
        return *this;
    }
    
    /*! \fn template <class U> SmallVec<T, D>& operator -= (const SmallVec<U, D>& rhs)
     *  \brief compound assignment operator -= */
    template <class U>
    SmallVec<T, D>& operator -= (const SmallVec<U, D>& rhs) {
        for (int i = 0; i < D; i++) {
            data[i] -= rhs.data[i];
        }
        return *this;
    }
    
    /*! \fn template <class U> SmallVec<T, D>& operator *= (const U rhs)
     *  \brief compound assignment operator *= */
    template <class U>
    SmallVec<T, D>& operator *= (const U rhs) {
        for (int i = 0; i < D; i++) {
            data[i] *= rhs;
        }
        return *this;
    }
    
    /*! \fn template <class U> SmallVec<T, D>& operator /= (const U rhs)
     *  \brief compound assignment operator /= */
    template <class U>
    SmallVec<T, D>& operator /= (const U rhs) {
        for (int i = 0; i < D; i++) {
            data[i] /= rhs;
        }
        return *this;
    }
    
    /*! \fn typename PromoteNumeric<T, double>::type Norm() const
     *  \brief calculate the norm of this vector, |x|, for at least double precision */
    typename PromoteNumeric<T, double>::type Norm() const {
        typename PromoteNumeric<T, double>::type sum = 0;
        for (int i = 0; i < D; i++) {
            sum += data[i] * data[i];
        }
        return sqrt(sum);
    }
    
    /*! \fn typename PromoteNumeric<T, double>::type Norm2() const
     *  \brief calculate the norm^2 of this vector, |x|^2, for at least double precision */
    typename PromoteNumeric<T, double>::type Norm2() const {
        typename PromoteNumeric<T, double>::type sum = 0;
        for (int i = 0; i < D; i++) {
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
        for (int i = 0; i < D; i++) {
            val = val && fabs(diff[i] < epsilon);
        }
        return val;
    }
    
    /*! \fn template <class U, class V> int relclose(const SmallVec<U, D>& rhs, const V epsilon) const
     *  \brief determinte (the relative difference between this_vector and rhs) divided by (the average length of them) is less than epsilon in each element */
    template <class U, class V>
    int RelClose(const SmallVec<U, D>& rhs, const V epsilon) const {
        SmallVec<typename PromoteNumeric<T, U>::type, D> sum, diff;
        for (int i = 0; i < D; i++) {
            sum[i] = fabs(data[i]) + fabs(rhs[i]);
        }
        diff = *this - rhs;
        bool val = true;
        for (int i = 0; i < D; i++) {
            val = val && ((2 * fabs(diff[i]) / sum[i]) < epsilon);
        } // if sum[i] = 0, we will get nan
        return val;
    }
    
    /*! \fn bool IsFinite() const
     *  \brief determine if every element is finite */
    bool IsFinite() const {
        bool val = true;
        for (int i = 0; i < D; i++) {
            val = val && std::isfinite(data[i]);
        }
        return val;
    }
    
    /*! \fn template <class U> bool operator == (const SmallVec<U, D>& rhs) const
     *  \brief overloading operator == */
    template <class U>
    bool operator == (const SmallVec<U, D>& rhs) const {
        bool val = true;
        for (int i = 0; i < D; i++) {
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
        for (int i = 0; i < D; i++) {
            data[i] = 0;
        }
        return *this;
    }
    
    /*! \fn template <class U, class V> bool InRange(U low, V high)
     *  \brief return true if all elements are on [a, b] */
    template <class U, class V>
    bool InRange(U low, V high) {
        bool val = true;
        for (int i = 0; i < D; i++) {
            val = val && (data[i] >= low && data[i] <= high);
        }
        return val;
    }
    
    /*! \fn template <class U, class V> bool InRange(SmallVec<U, D> low, SmallVec<V, D> high)
     *  \brief overloading InRange by taking SmallVec-type arguments instead of arithmetic type */
    template <class U, class V>
    bool InRange(SmallVec<U, D> low, SmallVec<V, D> high) {
        bool val = true;
        for (int i = 0; i < D; i++) {
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
    for (int i = 0; i < D; i++) {
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
    for (int i = 0; i < D; i++) {
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
    for (int i = 0; i < D; i++) {
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
    for (int i = 0; i < D; i++) {
        tmp.data[i] = lhs.data[i] / rhs;
    }
    return tmp;
}


/***********************************/
/********** MortonKey Part *********/
/***********************************/

/*
 * A MortonKey is a 128-bit (16 byte) number whose leftmost 32 bits are the particle index; the remaining 96 bits are the three integer coordinates interleaved into a Morton key.
 * N.B.: modern C++ offers bitset. If you are dealing with super large amout of particles (larger than 2^32 ~ 4.3 billion), then you might want to switch __uint128_t to bitset.
 */


/*! \struct template <int D> struct Orthant
 *  \brief orthant is the generalization of quadrant and octant and even hyperoctant in n>3 dimensions. This struct gives the directions in each dimension (+1 or -1). We will need them while building trees.
 *  \tparam D dimension of this vector */
template <int D>
struct Orthant {
    static const SmallVec<int, D> orthants[1<<D];
};

/*
 * Now define functions to output morton key
 */

/*! \fn template <typename T> void OutBinary(std::ostream &stream, T x)
 *  \brief output an integer number bit by bit */
template <typename T>
void OutBinary(std::ostream &stream, T x) {
    std::bitset<sizeof(T)*8> bits(x);
    stream << bits;
}

/*! \class BaseMortonKey
 *  \brief base class for MortonKey, use 128 bit key for all dimensions */
class BaseMortonKey {
private:
    /*
     * Below are useful constants worked for dilate3_32
     */
    
    /*! \var __uint128_t m1
     *  \brief binary: {0...63...0} 1 {0...63...0} 1 */
    __uint128_t m1;
    
    /*! \var __uint128_t m2
     *  \brief binary: {0...63...0} 1 {0...31...0} 1 {0...31...0} 1 */
    __uint128_t m2;
    
    /*! \var __uint128_t c1
     *  \brief binary: {1...32...1}{0...64...0}{1...32...1} */
    __uint128_t c1;
    
    /*! \var __uint128_t c2
     *  \brief binary: {0...16...0} {{1...16...1}{0...32...0}}x2 {1...16...1} */
    __uint128_t c2;
    
    /*! \var __uint128_t c3
     *  \brief binary: {{1...8...1}{0...16...0}}x5 {1...8...1} */
    __uint128_t c3;
    
    /*! \var __uint128_t c4
     *  \brief binary: {000011110000}x10 {00001111} */
    __uint128_t c4;
    
    /*! \var __uint128_t c5
     *  \brief binary: {110000}x21 {11} */
    __uint128_t c5;
    
    /*! \var __uint128_t c6
     *  \brief binary: {01} {001}x42 */
    __uint128_t c6;
    
    /*! \var __uint128_t upper32mask0
     *  \brief binary {0...32...0}{1...96...1} */
    __uint128_t upper32mask0;
    
    /*
     * Magic numbers for double2int: note that the double 0x1p+0(=1 in decimal) cannot be converted this way, so that the range of numbers is, strictly,  [0, 1). This two magic number make sure that the only 2^32 particles can be distinguished in one dimension. In other word, the minimum distance between two particles (or to be recognized as two particles) in one dimension case is 1.0/2^32 = 2.3283e-10.
     */
    
    /*! \var static constexpr double MAGIC = 6755399441055744.0
     *  \brief MAGIC = 2^52 + 2^51 = (0x0018000000000000)_16 */
    static constexpr double MAGIC = 6755399441055744.0;
    /*! \var static constexpr double MAXIMUMINTEGER = 4294967294.0
     *  \brief MAGIC = 2^32 - 2 = (0x00000000FFFFFFFE)_16 */
    static constexpr double MAXIMUMINTEGER = 4294967294.0;
    
public:
    // define name alias for intuitive definition
#ifdef OLDCPP
    typedef __uint128_t morton_key;
#else // OLDCPP
    using morton_key = __uint128_t;
#endif // OLDCPP
    
    /*! \fn BaseMortonKey()
     *  \brief constructor */
    BaseMortonKey();
    
    /*! \fn ~BaseMortonKey()
     *  \brief destructor */
    ~BaseMortonKey();
    
    /*! \fn __uint32_t Double2Int(double d)
     *  \brief convert a double on [0, 1) to an unsigned 32 bit integer */
    __uint32_t Double2Int(double d);
    
    /*! \fn void InitializeMortonConstants()
     *  \brief initialize constants used in future calculations */
    void InitializeMortonConstants();
    
    /*! \fn inline int Key8Level(morton_key &m_key, int &level)
     *  \brief extract info (three digits) of specific level from the 96-bit key */
    inline int Key8Level(morton_key m_key, int level);
    
    /*! \fn void OutKey(std::ostream &stream, morton_key m_key)
     *  \brief output the particle index and its key */
    void OutKey(std::ostream &stream, morton_key m_key);
    
    /*! \fn inline int ParIndex(morton_key m_key)
     *  \brief return the particle index from the Morton Key */
    inline int ParIndex(morton_key m_key);
    
    /*! \fn morton_key Dilate3_Int32(int pos)
     *  \brief spread the bits of pos 3 apart: i.e., {1011} becomes {001 000 001 001} */
    morton_key Dilate3_Int32(int pos);
    
};

/*
 * A functor, or a function object, is an object that can behave like a function. This is done by defining operator()() of the class. In this case, implement operator()() as a comparison function.
 */

/*! \struct AscendingMorton
 *  \brief define a functor similar to std::greater<T>() to compare the 96-bit morton key */
struct AscendingMorton {
    bool operator() (BaseMortonKey::morton_key x, BaseMortonKey::morton_key y) {
        return ( (x<<32) < (y<<32) );
    }
};

/*! \class template <int D> MortonKey
 *  \brief  */
template <int D>
class MortonKey : public BaseMortonKey {
private:
    
public:
    
#ifdef OLDCPP
    typedef SmallVec<int, D> ivec;
    typedef SmallVec<float, D> fvec;
    typedef SmallVec<double, D> dvec;
#else // OLDCPP
    using ivec = SmallVec<int, D>;
    using fvec = SmallVec<float, D>;
    using dvec = SmallVec<double, D>;
#endif // OLDCPP
    
    /*! \var dvec scale
     *  \brief scale the box length to 1 */
    dvec scale;
    
    /*! \var dvec boxmin, boxmax
     *  \brief the bounding box in user coordinates to be mapped to [0, 1)^D */
    dvec boxmin, boxmax;
    
    /*! \fn InitMortonKey(dvec __boxmin, dvec __boxmax);
     *  \brief initialize the space scale, set base for calculations of Morton Keys */
    void InitMortonKey(dvec __boxmin, dvec __boxmax) {
        boxmin = __boxmin;
        boxmax = __boxmax;
        for (int d = 0; d < D; d++) {
            scale[d] = 1.0 / (boxmax[d] - boxmin[d]);
        }
    }
    
    /*! \fn template <class U> Morton(SmallVec<U, D> pos, int index)
     *  \brief convert a position vector pos and particle index into a 128-bit Morton Key */
    template <class U>
    morton_key Morton(SmallVec<U, D> pos, int index) {
        dvec pos_scaled = pos - boxmin;
        for (int d = 0; d < D; d++) {
            pos_scaled[d] *= scale[d];
        }
        
        SmallVec<__uint32_t, D> int_pos;
        for (int d = 0; d < D; d++) {
            int_pos[d] = Double2Int(pos_scaled[d]);
        }
        return Morton(int_pos, index);
    }
    
    /*! \fn Morton(SmallVec<__uint32_t, D> pos, int index)
     *  \brief overloading Morton above for __uint32_t */
    morton_key Morton(SmallVec<__uint32_t, D> pos, int index) {
        morton_key result = (static_cast<__uint128_t>(index))<<96;
        for (int d = 0; d < D; d++) {
            result |= (Dilate3_Int32(pos[d])<<d);
        }
    }
    
};

/********************************/
/********** BHTree Part *********/
/********************************/





#endif /* tree_hpp */

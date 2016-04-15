//
//  tree.hpp
//  PLATO: PLAneTesimal locatOr
//
//  Created by Rixin Li on 3/11/16.
//  Copyright © 2016 Rixin Li. All rights reserved.
//
//  Some of the following code are based on programs written by Dr. Philip Pinto during his course (ASTR596) in fall 2015.
//  Descriptions and comments are written by Rixin Li.

#ifndef tree_hpp
#define tree_hpp

#include "global.hpp"

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
    
    /*! \fn template <class U> explicit SmallVec( const U& scalar)
     *  \brief overloading constructor to broadcast a scalar, e.g., unit_vec = SmallVec<double, 3>(1.0). For keyword "explicit", read explanations below */
    template <class U>
    explicit SmallVec(const U& scalar) {
        for (int i = 0; i != D; i++) {
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
        for (int i = 0; i != D; i++) {
            result[i] = -data[i];
        }
        return result;
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

/***********************************/
/********** Particle Part **********/
/***********************************/

/*! \class template <int D> Particle
 *  \brief data set for one particle
 *  \tparam D dimension of data */
template <int D>
class Particle {
private:
    
public:
    /*! \var SmallVec<double, D> r, v
     *  \brief position vector r and velocity vector v */
    SmallVec<double, D> pos, vec;
    
    /*! \var double m
     *  \brief particle mass in code unit */
    double mass;
    
    /*! \var double radius
     *  \brief particle radius */
    double radius;
    
    /*! \var int id
     *  \brief particle ID in total particle set */
    int id;
};

/*! \class template <int D> ParticleSet
 *  \brief data for the entire particle set
 *  \tparam D dimension of data */
template <int D>
class ParticleSet {
private:
    
public:
    /*! \var int num_particle
     *  \brief number of particles (must < 2^32-2) */
    __uint32_t num_particle;
    
    /*! \var int num_type
     *  \brief number of particle types */
    int num_type;
    
    /*! \var double coor_lim
     *  \brief coordinate limits for grid and domain.
     *  It is in the order of grid limits (x1l, x1u, x2l, x2u, x3l, x3u) and domain limits (x1dl, x1du, x2dl, x2du, x3dl, x3du), where l means lower limit, u means upper limit, d means domain */
    double coor_lim[12];
    
    /*! \var std::vector<double> type_info
     *  \brief info of different types in lis file */
    std::vector<double> type_info;
    
    /*! \var double time
     *  \brief current time in simulation */
    double time;
    
    /*! \var double dt
     *  \brief current time step */
    double dt;
    
    /*! \var Particle<D> particles
     *  \brief particle set */
    Particle<D> *particles;
    
    /*! \fn ParticleSet()
     *  \brief constructor */
    ParticleSet() : particles(nullptr) {}
    
    /*! \fn ~ParticleSet()
     *  \brief destructor */
    ~ParticleSet() {
        delete [] particles;
        particles = nullptr;
    }
    
    /*! \fn Particle<D> operator[] (const size_t i) const
     *  \brief define operator[] for easy access */
    Particle<D> operator[] (const size_t i) const {
        assert(i < num_particle);
        return *(particles+i);
    }
    
    /*! \fn Particle<D>& operator[] (const size_t i)
     *  \brief overload operator[] for element modification */
    Particle<D>& operator[] (const size_t i) {
        assert(i < num_particle);
        return *(particles+i);
    }
    
    /*! \fn void Reset()
     *  \brief delete particle set */
    void Reset() {
        delete [] particles;
        particles = nullptr;
    }
    
    /*! \fn void AllocateSpace(int N)
     *  \brief allcoate space for particles */
    void AllocateSpace(int N) {
        Reset();
        num_particle = N;
        particles = new Particle<D>[N];
    }
    
    /*! \fn void ReadLisFile(std::string filename)
     *  \brief read particle data from *.lis file */
    void ReadLisFile(std::string filename);
    
};

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
    
    /*! \fn inline int ParticleIndex(morton_key m_key)
     *  \brief return the particle index from the Morton Key */
    inline int ParticleIndex(morton_key m_key);
    
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
        for (int d = 0; d != D; d++) {
            scale[d] = 1.0 / (boxmax[d] - boxmin[d]);
        }
    }
    
    /*! \fn template <class U> Morton(SmallVec<U, D> pos, int index)
     *  \brief convert a position vector pos and particle index into a 128-bit Morton Key */
    template <class U>
    morton_key Morton(SmallVec<U, D> pos, int index) {
        dvec pos_scaled = pos - boxmin;
        for (int d = 0; d != D; d++) {
            pos_scaled[d] *= scale[d];
        }
        
        SmallVec<__uint32_t, D> int_pos;
        for (int d = 0; d != D; d++) {
            int_pos[d] = Double2Int(pos_scaled[d]);
        }
        return Morton(int_pos, index);
    }
    
    /*! \fn Morton(SmallVec<__uint32_t, D> pos, int index)
     *  \brief overloading Morton above for __uint32_t */
    morton_key Morton(SmallVec<__uint32_t, D> pos, int index) {
        morton_key result = (static_cast<__uint128_t>(index))<<96;
        for (int d = 0; d != D; d++) {
            result |= (Dilate3_Int32(pos[d])<<d);
        }
    }
    
};

/********************************/
/********** BHTree Part *********/
/********************************/

/*! \class template <int D> class BHtree : public MortonKey<D>
 *  \brief BHtree is the tree class for organizing particle data. In 3D, it's similar to octree.
 *  \tparam D dimension of this vector */
template <int D>
class BHtree : public MortonKey<D> {
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
    
    /*! \struct InternalParticle
     *  \brief necessary particle data for tree */
    struct InternalParticle {
        /*! \var dvec pos
         *  \brief the coordinates of particle position */
        dvec pos;
        
        /*! \var double mass
         *  \brief the particle mass */
        double mass;
        
        /*! \var __uint32_t id
         *  \brief particle index */
        __uint32_t id;
    };
    
    /*! \struct TreeNode
     *  \brief tree node structure */
    struct TreeNode {
        /*! \var dvec center
         *  \brief center coordinates of a node */
        dvec center;
        
        /*! \var double half_width
         *  \brief half the width of a node */
        double half_width;
        
        /*! \var __uint32_t begin
         *  \brief the beginning particle index in node */
        __uint32_t begin;
        
        /*! \var __uint32_t end
         *  \brief the ending particle index in node, notice that this "end" follows the C++ tradition and serves as the off-the-end iterator */
        __uint32_t end;
        
        /*! \var __uint32_t parent
         *  \brief the parent node's index */
        __uint32_t parent;
        
        /*! \var __uint32_t first_daughter;
         *  \brief the index of first daughter node */
        __uint32_t first_daughter;
        
        /*! \var __uint16_t orthant
         *  \brief orthant is the daughter direction from parent */
        __uint16_t orthant;
        
        /*! \var __uint8_t num_daughter
         *  \brief the number of daughters */
        __uint8_t num_daughter;
        
        /*! \var __uint8_t level
         *  \brief level in tree */
        __uint8_t level;
    };
    
    /*! \var static const int max_level = 32
     *  \brief the maximum levels of this tree is 32 */
    static const int max_level = 32;
    
    /*! \var typename MortonKey<D>::morton_key *morton
     *  \brief store all the morton key of particles */
    typename MortonKey<D>::morton_key *morton;
    
    /*! \var int num_particle
     *  \brief number of particles (must < 2^32-2) 
     *  This must < 2^32-1-1 = 0xffffffff-1, because (*TreeNode)->end means the off-the-end iterator, so we need one more number than the total number of particles. Anyway, if we are dealing with more particles than that, we should adapt our tools and use more advanced Morton Keys */
    __uint32_t num_particle;
    
    /*! \var InternalParticle *particle_list
     *  \brief a list of particle data */
    InternalParticle *particle_list;
    
    /*! \var TreeNode *tree;
     *  \brief the whole tree data */
    TreeNode *tree;
    
    /*! \var int num_leaf_set, *leaf_set
     *  \brief (the number of) leaf nodes */
    int num_leaf_nodes, *leaf_nodes;
    
    /*! \var int num_nodes, *node2leaf
     *  \brief TBD */
    int num_nodes, *node2leaf;
    
    /*! \var int max_leaf_size
     *  \brief max number of  */
    int max_leaf_size;
    
    /*! \var int max_daughters
     *  \brief max number of daughters 2^D */
    int max_daughters;
    
    /*! \var dvec root_center
     *  \brief the center coordinates of root node */
    dvec root_center;
    
    /*! \var int root
     *  \brief root node */
    int root;
    
    /*! \var int root_level
     *  \brief the level of root node */
    int root_level;
    
    /*! \var int node_ptr
     *  \brief use an integer as the pointer of particle (since the index of particle is int) */
    int node_ptr;
    
    /*! \var double half_width
     *  \brief half width of the whole tree structure */
    double half_width;
    
    /*! \var int level_count[max_level]
     *  \brief TBD */
    int level_count[max_level];
    
    /*! \var int level_ptr[max_level]
     *  \brief TBD */
    int level_ptr[max_level];
    
    /*! \var void *heaps;
     *  \brief TBD */
    void *heaps;
    
    /*! \var double const to_diagonal;
     *  \brief const used in SphereNodeIntersect */
    double const to_diagonal = sqrt(D);
    
    /*! \fn BHtree()
     *  \brief constructor, about the member initializer lists, refer to http://en.cppreference.com/w/cpp/language/initializer_list */
    BHtree() : tree(nullptr), morton(nullptr), leaf_nodes(nullptr), node2leaf(nullptr), particle_list(nullptr) {
        max_daughters = (1<<D);
        root = 0;
        root_level = 1;
        // heaps TBD
    }
    
    /*! \fn Reset()
     *  \brief release memory and reset */
    void Reset() {
        num_particle = 0;
        
        /*
         * Delete operator releases the memory of objects allocated by new operator; delete [] --> new []. Since delete operator will perform nullptr check first, we can safely use it without check. An expression with the delete[] operator, first calls the appropriate destructors for each element in the array (if these are of a class type), and then calls an array deallocation function (refer to http://www.cplusplus.com/reference/new/operator%20delete[]/).
         * With -std=c++11 and above, we should try to use shared_ptr for smart memory managements. Mark here and implement later.
         */
        
        delete [] particle_list;
        particle_list = nullptr;
        
        delete [] tree;
        tree = nullptr;
        
        delete [] morton;
        morton = nullptr;
        
        delete [] leaf_nodes;
        leaf_nodes = nullptr;
        
        delete [] node2leaf;
        node2leaf = nullptr;
    }
    
    /*! \fn ~BHtree()
     *  \brief destructor */
    ~BHtree() {
        Reset();
        // heaps TBD
    }
    
    /*! \fn void SortPoints()
     *  \brief sort points by morton key and then copy back to particle list */
    void SortPoints() {
        InternalParticle *tmp = new InternalParticle[num_particle];
        for (int i = 0; i != num_particle; i++) {
            tmp[i] = particle_list[ParticleIndex(morton[i])];
        }
        std::memcpy(particle_list, tmp, sizeof(InternalParticle)*num_particle);
        delete [] tmp;
    }
    
    /*! \fn void CountNodesLeaves(int const level, int __begin, int const __end)
     *  \brief traverse the tree and count nodes and leaves */
    void CountNodesLeaves(int const __level, int __begin, int const __end) { // double underscore here is to avoid confusion with TreeNode member or InternalParticle member
        int orthant = Key8Level(morton[__begin], __level);
        while ( (orthant < max_daughters) && (__begin < __end)) {
            int count = 0;
            while (__begin < __end) {
                if (Key8Level(morton[__begin], __level) == orthant ) {
                    __begin++;
                    count++;
                } else {
                    // already sorted, if false, then just break
                    break;
                }
            }
            assert(count > 0); // this is weird ???
            
            level_count[__level]++;
            num_nodes++;
            
            if (count <= max_leaf_size) {
                num_leaf_nodes++; // only nodes with leaves < max_leaf_size are leaves
            } else {
                CountNodesLeaves(__level+1, __begin-count, __begin-1);
            }
            
            if (__begin < __end) {
                orthant = Key8Level(morton[__begin], __level); // search next data
            }
        }
    }
    
    /*! \fn void FillTree(int const level, int __begin, int const __end, int const parent, dvec const center, double const __half_width)
     *  \brief fill the tree with data */
    void FillTree(int const __level, int __begin, int const __end, int const __parent, dvec const __center, double const __half_width) { // double underscore here is to avoid confusion with TreeNode member or InternalParticle member
        assert(__level < max_level);
        assert(__end > __begin); // note if this will cause bug
        assert(tree[__parent].first_daughter == 0);
        assert(tree[__parent].num_daughter == 0); // parent shouldn't have any daughters
        
        int orthant = Key8Level(morton[__begin], __level);
        while (__begin < __end) {
            assert( orthant < max_daughters);
            
            // count number of particles in this orthant
            int count = 0;
            while (__begin < __end) {
                if (Key8Level(morton[__begin], __level) == orthant) {
                    __begin++;
                    count++;
                } else {
                    break;
                }
            }
            assert(count > 0);
            
            // get daughter node number in tree
            int daughter = level_ptr[__level];
            level_ptr[__level]++;
            
            if (tree[__parent].first_daughter == 0) {
                // first daughter
                assert( tree[__parent].num_daughter = 0);
                tree[__parent].first_daughter = daughter;
                tree[__parent].num_daughter = 1;
            } else {
                // subsequent daughters
                tree[__parent].num_daughter++;
                assert(tree[__parent].num_daughter <= max_daughters);
            }
            
            TreeNode *p = &tree[daughter];
            p->level = __level + 1;
            p->parent = __parent;
            p->begin = __begin - count;
            p->end = __begin;
            p->half_width = __half_width;
            for (int d = 0; d != D; d++) {
                p->center[d] = __center[d] + __half_width * Orthant<D>::orthant[orthant][d];
            }
            p->orthant = orthant;
            p->first_daughter = 0;
            p->num_daughter = 0;
            node_ptr++;
            assert(node_ptr < num_nodes);
            
            if (count <= max_leaf_size) {
                // node with <= max_leaf_size particles is a leaf
                leaf_nodes[num_leaf_nodes] = daughter;
                node2leaf[daughter] = num_leaf_nodes;
                num_leaf_nodes++;
            } else {
                // node with > max_leaf_size particles is a branch
                FillTree(p->level, __begin-count, __begin-1, daughter, p->center, 0.5*__half_width);
            }
            
            // now next daughter of this node
            if (__begin < __end) {
                orthant = Key8Level(morton[__begin], __level);
            }
        }
    }
    
    /*! \fn void BuildTree(dvec const __center, double const __half_width, ParticleSet<D> const &particle_set, int const __max_leaf_size)
     *  \brief build tree from particle data */
    void BuildTree(dvec const __center, double const __half_width, ParticleSet<D> const &particle_set, int const __max_leaf_size) { // double underscore here is to avoid confusion with all sorts of members
        Reset();
        half_width = __half_width;
        root_center = __center;
        max_leaf_size = __max_leaf_size;
        
        assert(particle_set.num_particle < (__uint32_t)0xffffffff);
        num_particle = particle_set.num_particle;
        particle_list = new InternalParticle[num_particle];
        
        for (int i = 0; i != num_particle; i++) {
            particle_list[i].pos = particle_set[i].pos;
            particle_list[i].mass = particle_set[i].mass;
            particle_list[i].id = i; // original index of particles
        }
        
        // compute Morton Keys and sort particle_list by Morton order
        this->InitMortonKey(root_center-dvec(half_width), root_center+dvec(half_width));
        morton = new typename MortonKey<D>::morton_key[num_particle];
        
        for (int i = 0; i != num_particle; i++) {
            morton[i] = Morton(particle_list[i].pos, i);
        }
        std::sort(&(morton[0]), &(morton[num_particle]), AscendingMorton());
        for (int i = 0; i != num_particle-1; i++) {
            assert((morton[i]<<32) < (morton[i+1]<<32));
        }
        SortPoints();
        
        num_leaf_nodes = 0;
        level_count[0] = 1;
        for (int level = 1; level != max_level; level++) {
            level_count[level] = 0;
        }
        
        // run through the data once to determine space required for tree
        num_nodes = 1; // root contribute one
        CountNodesLeaves(root_level, 0, num_particle);
        
        assert(num_nodes == std::accumulate(level_count, level_count+max_level, 0));
        
        // allocate space for tree, leaf_nodes, and node2leaf index mapping
        node2leaf = new int[num_nodes];
        for (int i = 0; i != num_nodes; i++) {
            node2leaf[i] = -1;
        }
        leaf_nodes = new int[num_leaf_nodes];
        tree = new TreeNode[num_nodes];
        
        level_ptr[0] = 0;
        for (int level = 1; level != max_level; level++) {
            level_ptr[level] = level_ptr[level-1] + level_count[level-1];
        }
        node_ptr = 0;
        TreeNode *p = &tree[root];
        p->first_daughter = 0;
        p->orthant = 0;
        p->num_daughter = 0;
        p->level = root_level;
        p->center = root_center;
        p->half_width = half_width;
        p->begin = 0;
        p->end = num_particle;
        p->parent = -1;
        
        // run through the data again to build the tree
        num_leaf_nodes = 0;
        FillTree(root_level, 0, num_particle, root, root_center, 0.5*half_width);
        assert(node_ptr + 1 == num_nodes);
        delete [] morton;
    }
    
    /*! \fn bool Within(dvec const __pos, dvec const node_center, double const __half_width)
     *  \brief determine if a particle is within certain distance of a node center */
    bool Within(dvec const __pos, dvec const node_center, double const __half_width) {
        double epsilon = 1.0e-8;
        for (int d = 0; d != D; d++) {
            if ( !(__pos[d] >= node_center[d] - __half_width - epsilon && __pos[d] <= node_center[d] + __half_width + epsilon)) {
                return false;
            }
        }
        return true;
    }
    
    /*! \fn void CheckTree(int const node, int const __level, dvec const node_center, double const __half_width)
     *  \brief traverse the tree and check whether each point is within the node it is supposed to be. Also check levels/widths/centers */
    void CheckTree(int const node, int const __level, dvec const node_center, double const __half_width) {
        assert(tree[node].level == __level);
        assert(tree[node].half_width == __half_width);
        
        for (int p = tree[node].begin; p != tree[node].end; p++) {
            if (!Within(particle_list[p].pos, node_center, __half_width)) {
                io_ops->error_message << "Particle " << particle_list[p].pos << " outside node " << node_center << " with width " << 2*tree[node].half_width;
            }
        }
        
        for (int daughter = tree[node].first_daughter; daughter != tree[node].first_daughter + tree[node].num_daughter; daughter++) {
            dvec tmp_center = node_center;
            for (int d = 0; d != D; d++) {
                tmp_center[d] += 0.5 * __half_width * Orthant<D>::orthant[tree[daughter].orthant][d];
                dvec daughter_center = tree[daughter].center;
                assert(tmp_center == daughter_center);
                CheckTree(daughter, __level+1, daughter_center, 0.5*__half_width);
            }
        }
    }
    
    /*! \fn inline bool IsLeaf(int const node)
     *  \breif determine if a node is leaf */
    inline bool IsLeaf(int const node) {
        assert(node < num_nodes);
        return (tree[node].num_daughter == 0); // N.B., TreeNode has begin and end
    }
    
    /*! \fn inline int NodeSize(int const node)
     *  \brief return the number of particles in a node */
    inline int NodeSize(int const node) {
        assert(node < num_nodes);
        return tree[node].end - tree[node].begin; // note the end is off-the-end iterator
    }
    
    /*! \fn inline bool InNode(dvec const __pos, int const node)
     *  \brief determine if a particle is in a node */
    inline bool InNode(dvec const __pos, int const node) {
        return Within(__pos, tree[node].center, tree[node].half_width);
    }
    
    /*! \fn __uint32_t Key2Leaf(typename MortonKey<D>::morton_key const __morton, int const node, int const __level)
     *  \brief given a Morton Key, find the node number of the leaf cell where this key belongs */
    __uint32_t Key2Leaf(typename MortonKey<D>::morton_key const __morton, int const node, int const __level) {
        // if a leaf, just return answer
        if (IsLeaf(node)) {
            return node;
        }
        
        // else recurse into the correct direction
        int orthant = Key8Level(__morton, __level);
        int daughter = -1;
        for (int d = tree[node].first_daughter; d != tree[node].first_daughter + tree[node].num_daughter; d++) {
            if (tree[d].orthant == orthant) {
                daughter = Key2Leaf(__morton, d, __level+1);
                break;
            }
        }
        if (daughter == -1) {
            io_ops->error_message << "Key2Leaf: leaf cell doesn't exist in tree." << std::endl;
            assert(daughter >= 0);
        }
        return daughter;
    }
    
    /*! \fn __uint32_t Pos2Node(dvec const __pos)
     *  \brief given position, find the index of node containing it */
    __uint32_t Pos2Node(dvec const __pos) {
        assert( Within(__pos, root_center, half_width));
        typename MortonKey<D>::morton_key __morton = Morton(__pos, 0); // index doesn't matter here, just give 0
        return Key2Leaf(__morton, root, root_level);
    }
    
    /*! \fn bool SphereNodeIntersect(dvec const __center, double const r, int const node)
     *  \brief return true if any part of node is within the sphere (__center, r) */
    bool SphereNodeIntersect(dvec const __center, double const r, int const node) {
        assert(node < num_nodes);
        double c2c = (tree[node].center - __center).Norm2();
        double tmp_distance = tree[node].half_width * to_diagonal + r;
        
        // check if node is outside the sphere
        if (c2c > tmp_distance * tmp_distance) {
            return false;
        }
        
        // check if node center is inside the sphere
        if (c2c < r || c2c < tree[node].half_width) {
            return true;
        }
        
        // now do exact check for intersection
        // notice that when we extend each side, the space is divided into multiple sub-space, and the value for each sub-space is different
        dvec pos_min = tree[node].center - dvec(tree[node].half_width);
        dvec pos_max = tree[node].center + dvec(tree[node].half_width);
        
        double mindist2 = 0;
        for (int d = 0; d != D; d++) {
            if (__center[d] < pos_min[d]) {
                mindist2 += (__center[d] - pos_min[d]) * (__center[d] - pos_min[d]);
            } else if (__center[d] > pos_max[d]) {
                mindist2 += (__center[d] - pos_max[d]) * (__center[d] - pos_max[d]);
            }
        }
        
        return mindist2 <= r*r;
    }
    
    
    
    
    
    
    
};
































#endif /* tree_hpp */

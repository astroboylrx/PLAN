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
    
    
    /*
     * Warning: Keep extra mind and cautions here. I haven't fully understand the usage of enable_if.
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
        progIO->error_message << "Error: Use an array with wrong size to construct class SmallVec." << std::endl;
        progIO->Output(std::cerr, progIO->error_message, __even_more_output, __all_processors);
        exit(2);
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
        progIO->error_message << "Error: Use an array with wrong size to construct class SmallVec." << std::endl;
        progIO->Output(std::cerr, progIO->error_message, __even_more_output, __all_processors);
        exit(2);
    }
    
    template <class U, typename std::enable_if<std::is_pointer<U>::value, int>::type = 0>
    SmallVec(U arg) {
        progIO->error_message << "Error: Please don't provide a pointer to construct class SmallVec. Whether it is a pointer to scalar or a pointer to array remains unknown. " << std::endl;
        progIO->Output(std::cerr, progIO->error_message, __even_more_output, __all_processors);
        exit(2);
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
    
    template <class U>
    SmallVec<T, D>& operator = (const SmallVec<U, D>& rhs) {
        for (int i = 0; i != D; i++) {
            data[i] = rhs.data[i];
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
    SmallVec<float, D> pos, vec;
    
    /*! \var int property_index
     *  \brief index of particle properties */
    int property_index;
    
    /*! \var float density
     *  \brief local particle density */
    float density;
    
    /*! \var int id
     *  \brief particle ID in total particle set */
    int id;
};

/*! \class template <int D> Plantesimal
 *  \brief data for one plantesimal structure */
template <int D>
class Plantesimal {
private:
    
public:
    /*! \var std::vector<__uint32_t> indices
     *  \brief indices of particles that belongs to this plantesimals */
    std::vector<__uint32_t> indices;
    
    /*! \var SmallVec<double, D> center
     *  \brief geometric center */
    SmallVec<double, D> center;
    
    /*! \var double R_Hill
     *  \brief Hill Radius */
    double R_Hill;
    
    /*! \var double mass
     *  \brief Hill Radius */
    double mass;
    
    /*! \var double Jz
     *  \brief angular momentum in z-direction */
    double Jz;
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
    
    /*! \var float coor_lim
     *  \brief coordinate limits for grid and domain.
     *  It is in the order of grid limits (x1l, x1u, x2l, x2u, x3l, x3u) and domain limits (x1dl, x1du, x2dl, x2du, x3dl, x3du), where l means lower limit, u means upper limit, d means domain */
    float coor_lim[12];
    
    /*! \var std::vector<float> type_info
     *  \brief info of different types in lis file */
    std::vector<float> type_info;
    
    /*! \var float time
     *  \brief current time in simulation */
    float time;
    
    /*! \var float dt
     *  \brief current time step */
    float dt;
    
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
    
    /*! \fn void ReadMultipleLisFile(std::vector<std::string>::iterator begin, std::vector<std::string>::iterator end)
     *  \brief read particle data from a series of *.lis file created by each processor
     *  Assume that one cpu core read one snapshot once */
    void ReadMultipleLisFile(std::vector<std::string>::iterator begin, std::vector<std::string>::iterator end) {
        std::ifstream lis_file;
        long tmp_num_particle;
        
        // First step, obtain the box limit from RootMin and RootMax and count the total particle numbers
        lis_file.open(begin->c_str(), std::ios::binary);
        if (lis_file.is_open()) {
            lis_file.read(reinterpret_cast<char*>(coor_lim), 12*sizeof(float));
            progIO->log_info << *begin << ", x1l = " << coor_lim[0] << ", x1u = " << coor_lim[1]
            << ", x2l = " << coor_lim[2] << ", x2u = " << coor_lim[3]
            << ", x3l = " << coor_lim[4] << ", x3u = " << coor_lim[5]
            << ", x1dl = " << coor_lim[6] << ", x1du = " << coor_lim[7]
            << ", x2dl = " << coor_lim[8] << ", x2du = " << coor_lim[9]
            << ", x3dl = " << coor_lim[10] << ", x3du = " << coor_lim[11] << "\n";
            lis_file.read(reinterpret_cast<char*>(&num_type), sizeof(int));
            progIO->log_info << "num_type = " << num_type;
            type_info.resize(num_type);
            
            for (int i = 0; i != num_type; i++) {
                lis_file.read(reinterpret_cast<char*>(&type_info[i]), sizeof(float));
                progIO->log_info << ": type_info[" << i << "] = " << type_info[i];
            }
            progIO->log_info << "; || ";
            lis_file.read(reinterpret_cast<char*>(&time), sizeof(float));
            lis_file.read(reinterpret_cast<char*>(&dt), sizeof(float));
            progIO->log_info << "time = " << time << ", dt = " << dt;
            lis_file.read(reinterpret_cast<char*>(&tmp_num_particle), sizeof(long));
            num_particle = static_cast<__uint32_t>(tmp_num_particle);
            lis_file.close();
        } else { // if (lis_file.is_open())
            progIO->error_message << "Error: Failed to open file " << begin->c_str() << std::endl;
            progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
            exit(3); // cannot open file
        }
        
        for (std::vector<std::string>::iterator it = (begin+1); it != end; it++) {
            lis_file.open(it->c_str(), std::ios::binary);
            if (lis_file.is_open()) {
                lis_file.seekg((14+num_type)*sizeof(float)+sizeof(int), std::ios::beg);
                lis_file.read(reinterpret_cast<char*>(&tmp_num_particle), sizeof(long));
                num_particle += static_cast<__uint32_t>(tmp_num_particle);
                lis_file.close();
            } else {
                progIO->error_message << "Error: Failed to open file " << it->c_str() << std::endl;
                progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                exit(3); // cannot open file
            }
        }
        
        progIO->log_info << ", total num_particle = " << num_particle << std::endl;
        progIO->Output(std::clog, progIO->log_info, __even_more_output, __all_processors);
        
        AllocateSpace(num_particle);
        
        // Second step, extend box limit to include ghost zone
        
        
        // Thrid step, read particle data and duplicate ghost particles
        __uint32_t tmp_id = 0;
        Particle<D> *p;
        size_t triple_float = 3 * sizeof(float);
        size_t one_float = sizeof(float);
        size_t one_fil = sizeof(float) + sizeof(int) + sizeof(long);
        
        for (std::vector<std::string>::iterator it = begin; it != end; it++) {
            lis_file.open(it->c_str(), std::ios::binary);
            if (lis_file.is_open()) {
                lis_file.seekg((14+num_type)*sizeof(float)+sizeof(int), std::ios::beg);
                lis_file.read(reinterpret_cast<char*>(&tmp_num_particle), sizeof(long));
                std::stringstream content;
                content << lis_file.rdbuf();
                std::string tmp_str = content.str();
                const char *tmp_char = tmp_str.data();
                for (int i = 0; i != tmp_num_particle; i++) {
                    p = &particles[tmp_id];
                    std::memcpy((char*)p->pos.data, tmp_char, triple_float);
                    std::advance(tmp_char, triple_float);
                    std::memcpy((char*)p->vec.data, tmp_char, triple_float);
                    std::advance(tmp_char, triple_float);
                    std::memcpy((char*)&p->density, tmp_char, one_float);
                    std::advance(tmp_char, one_float);
                    std::memcpy((char*)&p->property_index, tmp_char, one_float);
                    std::advance(tmp_char, one_fil);
                    p->id = tmp_id++;
                }
                lis_file.close();
            } else { // if (lis_file.is_open())
                progIO->error_message << "Error: Failed to open file " << it->c_str() << std::endl;
                progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                exit(3); // cannot open file
            }
        }
        
        progIO->log_info << "id = " << particles[567].id << ", property_index = " << particles[567].property_index << ", rad = " << particles[567].density << ", pos = " << particles[567].pos << ", v = " << particles[567].vec << std::endl;
        progIO->Output(std::clog, progIO->log_info, __even_more_output, __all_processors);
        
    }
    
    /*! \fn void ReadSingleLisFile(std::vector<std::string>::iterator it)
     *  \brief read particle data from one combined lis file from all processors */
    void ReadSingleLisFile(std::vector<std::string>::iterator it) {
        std::ifstream lis_file;
        long tmp_num_particle;
        
        lis_file.open(it->c_str(), std::ios::binary);
        if (lis_file.is_open()) {
            lis_file.read(reinterpret_cast<char*>(coor_lim), 12*sizeof(float));
            progIO->log_info << *it << ", x1l = " << coor_lim[0] << ", x1u = " << coor_lim[1]
            << ", x2l = " << coor_lim[2] << ", x2u = " << coor_lim[3]
            << ", x3l = " << coor_lim[4] << ", x3u = " << coor_lim[5]
            << ", x1dl = " << coor_lim[6] << ", x1du = " << coor_lim[7]
            << ", x2dl = " << coor_lim[8] << ", x2du = " << coor_lim[9]
            << ", x3dl = " << coor_lim[10] << ", x3du = " << coor_lim[11] << "\n";
            lis_file.read(reinterpret_cast<char*>(&num_type), sizeof(int));
            progIO->log_info << "num_type = " << num_type;
            type_info.resize(num_type);
            for (int i = 0; i != num_type; i++) {
                lis_file.read(reinterpret_cast<char*>(&type_info[i]), sizeof(float));
                progIO->log_info << ": type_info[" << i << "] = " << type_info[i];
            }
            progIO->log_info << "; || ";
            lis_file.read(reinterpret_cast<char*>(&time), sizeof(float));
            lis_file.read(reinterpret_cast<char*>(&dt), sizeof(float));
            progIO->log_info << "time = " << time << ", dt = " << dt;
            lis_file.read(reinterpret_cast<char*>(&tmp_num_particle), sizeof(long));
            num_particle = static_cast<__uint32_t>(tmp_num_particle);
            progIO->log_info << ", total num_particle = " << num_particle << std::endl;
            progIO->Output(std::clog, progIO->log_info, __even_more_output, __all_processors);
            
            AllocateSpace(num_particle);
            
            // Thrid step, read particle data and duplicate ghost particles
            __uint32_t tmp_id = 0;
            Particle<D> *p;
            size_t triple_float = 3 * sizeof(float);
            size_t one_float = sizeof(float);
            size_t one_int = sizeof(int);
            size_t one_fil = one_float + one_int + sizeof(long);
            
            std::stringstream content;
            content << lis_file.rdbuf();
            std::string tmp_str = content.str();
            const char *tmp_char = tmp_str.data();
            for (int i = 0; i != tmp_num_particle; i++) {
                p = &particles[tmp_id];
                std::memcpy((char*)p->pos.data, tmp_char, triple_float);
                std::advance(tmp_char, triple_float);
                std::memcpy((char*)p->vec.data, tmp_char, triple_float);
                std::advance(tmp_char, triple_float);
                std::memcpy((char*)&p->density, tmp_char, one_float);
                std::advance(tmp_char, one_float);
                std::memcpy((char*)&p->property_index, tmp_char, one_int);
                std::advance(tmp_char, one_fil);
                p->id = tmp_id++;
            }
            lis_file.close();
        } else { // if (lis_file.is_open())
            progIO->error_message << "Error: Failed to open file " << it->c_str() << std::endl;
            progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
            exit(3); // cannot open file
        }
        
        progIO->log_info << "id = " << particles[567].id << ", property_index = " << particles[567].property_index << ", rad = " << particles[567].density << ", pos = " << particles[567].pos << ", v = " << particles[567].vec << std::endl;
        progIO->Output(std::clog, progIO->log_info, __even_more_output, __all_processors);
    }
    
    
    /*! \fn void ReadLisFile(int loop_count)
     *  \brief read particle data from *.lis file */
    void ReadLisFile(int loop_count) {
        timer[__tmp_used_timer].StartTimer();
        if (progIO->flags.combined_flag) {
            ReadSingleLisFile(progIO->file_name.lis_data_file_name.begin()+loop_count);
        } else {
            std::vector<std::string>::iterator file_head = progIO->file_name.lis_data_file_name.begin();
            ReadMultipleLisFile(file_head + loop_count * progIO->num_cpu,
                                file_head + loop_count * progIO->num_cpu + progIO->num_cpu);
        }
        timer[__tmp_used_timer].StopTimer();
        progIO->log_info << "Reading loop_count (" << loop_count << ") cost " << timer[__tmp_used_timer].GiveTime() << " seconds\n";
        progIO->Output(std::clog, progIO->log_info, __more_output, __all_processors);
        
        progIO->physical_quantities[loop_count].time = time;
        progIO->physical_quantities[loop_count].dt = dt;
        
        // initialize particle scale heights
        progIO->physical_quantities[loop_count].particle_scale_height.resize(num_type);
        for (auto &item : progIO->physical_quantities[loop_count].particle_scale_height) {
            item = 0;
        }
        // calculate particle scale heights and find out maximum particle density
        for (__uint32_t i = 0; i != num_particle; i++) {
            progIO->physical_quantities[loop_count].particle_scale_height[particles[i].property_index] += particles[i].pos[2] * particles[i].pos[2];
            progIO->physical_quantities[loop_count].max_particle_density = std::max(progIO->physical_quantities[loop_count].max_particle_density, particles[i].density);
        }
        progIO->log_info << "loop_count = " << loop_count << ", max(rho_p) = " << progIO->physical_quantities[loop_count].max_particle_density;
        int tmp_index = 0;
        for (auto &item : progIO->physical_quantities[loop_count].particle_scale_height) {
            item = std::sqrt(item/num_particle);
            progIO->log_info << ", H_p[" << tmp_index++ << "] = " << item;
        }
        progIO->log_info << std::endl;
        progIO->Output(std::clog, progIO->log_info, __more_output, __all_processors);
        
    }
    
    /*! \fn void MakeGhostParticles(const double ghost_zone_ratio)
     *  \brief make ghost particles for ghost zone */
    void MakeGhostParticles(const double ghost_zone_ratio) {
        
    }
    
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
    /*! \var using morton_key = __uint128_t
     *  \brief define a type equivalent to __uint128_t for morton key */
    using morton_key = __uint128_t;
    
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
    inline int Key8Level(morton_key m_key, int level) {
        int shr = 93 - 3 * (level - 1);
        return (m_key>>shr) & 7UL; // 7UL = {0...60...0}{0111}
    }
    
    /*! \fn void OutKey(std::ostream &stream, morton_key m_key)
     *  \brief output the particle index and its key */
    void OutKey(std::ostream &stream, morton_key m_key);
    
    /*! \fn inline int ParticleIndex(morton_key m_key)
     *  \brief return the particle index from the Morton Key */
    inline int ParticleIndex(morton_key m_key) {
        return (m_key>>96);
    }
    
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
    /*! \var using ivec = SmallVec<int, D>
     *  \brief define a vector of int type */
    using ivec = SmallVec<int, D>;
    
    /*! \var using fvec = SmallVec<float, D>
     *  \brief define a vector of float type */
    using fvec = SmallVec<float, D>;
    
    /*! \var using dvec = SmallVec<double, D>
     *  \brief define a vector of double type */
    using dvec = SmallVec<double, D>;
    
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
        return result;
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
    /*! \var using ivec = SmallVec<int, D>
     *  \brief define a vector of int type */
    using ivec = SmallVec<int, D>;
    
    /*! \var using fvec = SmallVec<float, D>
     *  \brief define a vector of float type */
    using fvec = SmallVec<float, D>;
    
    /*! \var using dvec = SmallVec<double, D>
     *  \brief define a vector of double type */
    using dvec = SmallVec<double, D>;
    
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
    
    /*! \var std::vector<std::pair<int, double>> heaps;
     *  \brief a heap stores k-nearest neighbours */
    std::vector<std::pair<int, double>> heaps;
    
    /*! \var const double to_diagonal;
     *  \brief const used in SphereNodeIntersect */
    const double to_diagonal = sqrt(D);
    
    /*! \fn BHtree()
     *  \brief constructor, about the member initializer lists, refer to http://en.cppreference.com/w/cpp/language/initializer_list */
    BHtree() : tree(nullptr), morton(nullptr), leaf_nodes(nullptr), node2leaf(nullptr), particle_list(nullptr) {
        max_daughters = (1<<D);
        root = 0;
        root_level = 1;
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
        ClearHeaps();
    }
    
    /*! \fn void SortPoints()
     *  \brief sort points by morton key and then copy back to particle list */
    void SortPoints() {
        InternalParticle *tmp = new InternalParticle[num_particle];
        for (int i = 0; i != num_particle; i++) {
            tmp[i] = particle_list[this->ParticleIndex(morton[i])];
        }
        std::memcpy(particle_list, tmp, sizeof(InternalParticle)*num_particle);
        delete [] tmp;
    }
    
    /*! \fn void CountNodesLeaves(const int level, int __begin, const int __end)
     *  \brief traverse the tree and count nodes and leaves */
    void CountNodesLeaves(const int __level, int __begin, const int __end) { // double underscore here is to avoid confusion with TreeNode member or InternalParticle member
        int orthant = this->Key8Level(morton[__begin], __level);
        while ( (orthant < max_daughters) && (__begin < __end)) {
            int count = 0;
            while (__begin < __end) {
                if (this->Key8Level(morton[__begin], __level) == orthant ) {
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
                orthant = this->Key8Level(morton[__begin], __level); // search next data
            }
        }
    }
    
    /*! \fn void FillTree(const int level, int __begin, const int __end, const int parent, const dvec center, const double __half_width)
     *  \brief fill the tree with data */
    void FillTree(const int __level, int __begin, const int __end, const int __parent, const dvec __center, const double __half_width) { // double underscore here is to avoid confusion with TreeNode member or InternalParticle member
        assert(__level < max_level);
        assert(__end > __begin); // note if this will cause bug
        assert(tree[__parent].first_daughter == 0);
        assert(tree[__parent].num_daughter == 0); // parent shouldn't have any daughters
        
        int orthant = this->Key8Level(morton[__begin], __level);
        while (__begin < __end) {
            assert( orthant < max_daughters);
            
            // count number of particles in this orthant
            int count = 0;
            while (__begin < __end) {
                if (this->Key8Level(morton[__begin], __level) == orthant) {
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
                assert(tree[__parent].num_daughter = 0);
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
                p->center[d] = __center[d] + __half_width * Orthant<D>::orthants[orthant][d];
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
                orthant = this->Key8Level(morton[__begin], __level);
            }
        }
    }
    
    /*! \fn void BuildTree(const dvec __center, const double __half_width, ParticleSet<D> &particle_set, const int __max_leaf_size)
     *  \brief build tree from particle data */
    void BuildTree(const dvec __center, const double __half_width, ParticleSet<D> &particle_set, const int __max_leaf_size) { // double underscore here is to avoid confusion with all sorts of members
        Reset();
        half_width = __half_width;
        root_center = __center;
        max_leaf_size = __max_leaf_size;
        
        assert(particle_set.num_particle < (__uint32_t)0xffffffff);
        num_particle = particle_set.num_particle;
        particle_list = new InternalParticle[num_particle];
        
        for (int i = 0; i != num_particle; i++) {
            particle_list[i].pos = particle_set[i].pos;
            particle_list[i].mass = particle_set[i].property_index;
            particle_list[i].id = i; // original index of particles
        }
        
        // compute Morton Keys and sort particle_list by Morton order
        this->InitMortonKey(root_center-dvec(half_width), root_center+dvec(half_width));
        morton = new typename MortonKey<D>::morton_key[num_particle];
        
        for (int i = 0; i != num_particle; i++) {
            morton[i] = this->Morton(particle_list[i].pos, i);
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
    
    /*! \fn bool Within(const dvec __pos, const dvec node_center, const double __half_width)
     *  \brief determine if a particle is within certain distance of a node center */
    bool Within(const dvec __pos, const dvec node_center, const double __half_width) {
        double epsilon = 1.0e-8;
        for (int d = 0; d != D; d++) {
            if ( !(__pos[d] >= node_center[d] - __half_width - epsilon && __pos[d] <= node_center[d] + __half_width + epsilon)) {
                return false;
            }
        }
        return true;
    }
    
    /*! \fn void CheckTree(const int node, const int __level, const dvec node_center, const double __half_width)
     *  \brief traverse the tree and check whether each point is within the node it is supposed to be. Also check levels/widths/centers */
    void CheckTree(const int node, const int __level, const dvec node_center, const double __half_width) {
        assert(tree[node].level == __level);
        assert(tree[node].half_width == __half_width);
        
        for (int p = tree[node].begin; p != tree[node].end; p++) {
            if (!Within(particle_list[p].pos, node_center, __half_width)) {
                progIO->error_message << "Particle " << particle_list[p].pos << " outside node " << node_center << " with width " << 2*tree[node].half_width;
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
    
    /*! \fn inline bool IsLeaf(const int node)
     *  \breif determine if a node is leaf */
    inline bool IsLeaf(const int node) {
        assert(node < num_nodes);
        return (tree[node].num_daughter == 0); // N.B., TreeNode has begin and end
    }
    
    /*! \fn inline int NodeSize(const int node)
     *  \brief return the number of particles in a node */
    inline int NodeSize(const int node) {
        assert(node < num_nodes);
        return tree[node].end - tree[node].begin; // note the end is off-the-end iterator
    }
    
    /*! \fn inline bool InNode(const dvec __pos, const int node)
     *  \brief determine if a particle is in a node */
    inline bool InNode(const dvec __pos, const int node) {
        return Within(__pos, tree[node].center, tree[node].half_width);
    }
    
    /*! \fn __uint32_t Key2Leaf(typename MortonKey<D>::morton_key const __morton, const int node, const int __level)
     *  \brief given a Morton Key, find the node number of the leaf cell where this key belongs */
    __uint32_t Key2Leaf(typename MortonKey<D>::morton_key const __morton, const int node, const int __level) {
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
            progIO->error_message << "Key2Leaf: leaf cell doesn't exist in tree." << std::endl;
            assert(daughter >= 0);
        }
        return daughter;
    }
    
    /*! \fn __uint32_t Pos2Node(const dvec __pos)
     *  \brief given position, find the index of node containing it */
    __uint32_t Pos2Node(const dvec __pos) {
        assert( Within(__pos, root_center, half_width));
        typename MortonKey<D>::morton_key __morton = Morton(__pos, 0); // index doesn't matter here, just give 0
        return Key2Leaf(__morton, root, root_level);
    }
    
    /*! \fn bool SphereNodeIntersect(const dvec __center, const double r, const int node)
     *  \brief return true if any part of node is within the sphere (__center, r) */
    bool SphereNodeIntersect(const dvec __center, const double r, const int node) {
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
    
    /*! \struct template <typename T1, typename T2> struct less_second
     *  \brief served as comparison method for heaps */
    template <typename T1, typename T2>
    struct less_second {
        typedef std::pair<T1, T2> type;
        bool operator ()(type const& a, type const& b) const {
            return a.second < b.second;
        }
    };
    
    /*! \fn void Add2Heaps(const int knn, const int i, const double dr2)
     *  \brief add element to heaps */
    void Add2Heaps(const int knn, const int i, const double dr2) {
        if (heaps.size() < knn) {
            heaps.push_back(std::pair<int, double>(i, dr2));
            std::push_heap(heaps.begin(), heaps.end(), less_second<int, double>());
        } else {
            if (dr2 < heaps.front().second) {
                std::pop_heap(heaps.begin(), heaps.end(), less_second<int, double>());
                heaps.pop_back();
                heaps.push_back(std::pair<int, double>(i, dr2));
                std::push_heap(heaps.begin(), heaps.end(), less_second<int, double>());
            }
        }
    }
    
    /*! \fn inline void ClearHeaps()
     *  \brief release memory of heaps */
    inline void ClearHeaps() {
        // force clear and reallocation
        std::vector<std::pair<int, double>>().swap(heaps);
    }
    
    /*! \fn void RecursiveKNN(const dvec __pos, const int node, const double dist, const int knn)
     *  \brief do recursive KNN search to traverse the tree */
    void RecursiveKNN(const dvec __pos, const int node, const double dist, const int knn) {
        if (SphereNodeIntersect(__pos, dist, node)) {
            if (IsLeaf(node)) {
                for (int p = tree[node].begin; p != tree[node].end; p++) {
                    Add2Heaps(knn, p, (__pos-particle_list[p].pos).Norm2());
                }
            } else {
                for (int d = tree[node].first_daughter; d != tree[node].first_daughter + tree[node].num_daughter; d++) {
                    RecursiveKNN(__pos, d, dist, knn);
                }
            }
        }
    }
    
    /*! \fn void KNN_Search(const dvec __pos, const int knn, double &radius_knn, int *indices)
     *  \brief given position, perform k-nearest neighbours search and return radius and particle indices */
    void KNN_Search(const dvec __pos, const int knn, double &radius_knn, int *indices) {
        assert(knn <= num_particle);
        
        if (heaps.size() != 0) {
            ClearHeaps(); // clear memory first
        }
        heaps.reserve(knn);
        
        double max_dr2 = 0;
        // obatin a estimated distance firstly
        if (Within(__pos, root_center, half_width)) {
            // if within the box
            int node = Pos2Node(__pos);
            if (NodeSize(node) == 1) {
                node = tree[node].parent;
            }
            for (int p = tree[node].begin; p != tree[node].end; p++) {
                max_dr2 = std::max(max_dr2, (__pos - particle_list[p].pos).Norm2());
            }
        } else {
            // if point is outside the entire box, use the min distance to box boundary
            for (int d = 0; d < D; d++) {
                double dx = MaxOf(root_center[d]-half_width - __pos[d], 0.0, __pos[d] - root_center[d] - half_width);
                max_dr2 += dx * dx;
            }
        }
        
        double max_dr = sqrt(max_dr2);
        // now use the distance guess to proceed
        do {
            ClearHeaps();
            heaps.reserve(knn);
            TraverseKNN(__pos, root, max_dr, knn);
            max_dr2 *= 2;
        } while (heaps.size() < knn);
        
        // Phil did TraverseKNN once more, but I don't see the necessity...
        
        // get original particle id
        for (int i = 0; i != heaps.size(); i++) {
            indices[i] = particle_list[heaps[i].first].id;
        }
        radius_knn = sqrt(heaps.front().second);
    }
    
    /*! \fn void RecursiveBallSearch(const dvec __pos, int node, const double radius, int *indices, int &count)
     *  \brief do recursive ball search to traverse the tree */
    void RecursiveBallSearch(const dvec __pos, int node, const double radius, int *indices, int &count) {
        if (SphereNodeIntersect(__pos, radius, node)) {
            if (IsLeaf(node)) {
                for (int p = tree[node].begin; p != tree[node].end; p++) {
                    if ((__pos - particle_list[p].pos).Norm2() <= radius * radius) {
                        indices[count++] = p;
                    }
                }
            } else {
                for (int d = tree[node].first_daughter; d != tree[node].first_daughter + tree[node].num_daughter; d++) {
                    RecursiveBallSearch(__pos, d, radius, indices, count);
                }
            }
        }
    }
    
    /*! \fn void BallSearch(const dvec __center, const double radius, int *indices, int &count)
     *  \brief perform a search within a sphere */
    void BallSearch (const dvec __center, const double radius, int *indices, int &count) {
        count = 0;
        RecursiveBallSearch(__center, root, radius, indices, count);
        // convert morton key to original particle id
        for (int i = 0; i < count; i++) {
            indices[i] = particle_list[indices[i]].id;
        }
    }
    
    /*! \fn void FindPlanetesimals()
     *  \brief find plantesimals inside data */
    void FindPlanetesimals() {
        // Get high density regions by sorting dpar
        
        
        // Find neighbors + removing from list
        
        
    }
    
    
    
    
};
































#endif /* tree_hpp */

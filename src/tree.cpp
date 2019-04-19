//
//  tree.cpp
//  PLAN: PLanetesimal ANalyzer
//
//  Created by Rixin Li on 3/11/16.
//  Copyright © 2016 Rixin Li. All rights reserved.
//

/*! \file tree.cpp
 *  \brief contains function definitions for tree-related class */

#include "tree.hpp"

/***********************************************/
/***** template specialization for Orthant *****/
/***********************************************/

/*! \var template <> const SmallVec<int, 3> Orthant<3>::orthants[1U<<3]
 *  \brief constant orthants used under 3D */
template <> // template specialization need template<>. (reference: http://en.cppreference.com/w/cpp/language/template_specialization )
const SmallVec<int, 3> Orthant<3>::orthants[1U<<3] = {
    {-1, -1, -1},
    { 1, -1, -1},
    {-1,  1, -1},
    { 1,  1, -1},
    {-1, -1,  1},
    { 1, -1,  1},
    {-1,  1,  1},
    { 1,  1,  1}
};

/*
 * template specialization shouldn't appear in the header files since it will cause "duplicate symbol" error.
 * But if you inline it, putting into header files is fine. 
 */

/*! \var template <> const SmallVec<int, 2> Orthant<2>::orthants[1U<<2]
 *  \brief constant orthants used under 2D */
template <>
const SmallVec<int, 2> Orthant<2>::orthants[1U<<2] = {
    {-1, -1},
    { 1, -1},
    {-1,  1},
    { 1,  1}
};

/*! \var template <> const SmallVec<int, 1> Orthant<1>::orthants[1U<<1]
 *  \brief constant orthants used under 1D */
template<>
const SmallVec<int, 1> Orthant<1>::orthants[1U<<1] = {
    -1, 1
};

/*************************************************/
/***** template specialization for OutBinary *****/
/*************************************************/

/*! \fn template<> void OutBinary<BaseMortonKey::uint128_t>(std::ostream &stream, BaseMortonKey::uint128_t x)
 *  \brief make a template specialization for uint128_t since the initilization of bitset only support 64-bit integer */
template<>
void OutBinary<BaseMortonKey::uint128_t>(std::ostream &stream, BaseMortonKey::uint128_t x)
{
    OutBinary<uint64_t>(stream, static_cast<uint64_t>(x>>64));
    OutBinary<uint64_t>(stream, static_cast<uint64_t>(x));
}


/*! \fn template<> void OutBinary<float>(std::ostream &stream, float x)
 *  \brief make a template specialization for float */
template<>
void OutBinary<float>(std::ostream &stream, float x)
{
    uint32_t *bits = reinterpret_cast<uint32_t *>(&x);
    OutBinary<uint32_t>(stream, *bits);
}

/*! \fn template<> void OutBinary<double>(std::ostream &stream, double x)
 *  \brief make a template specialization for double */
template<>
void OutBinary<double>(std::ostream &stream, double x)
{
    uint64_t *bits = reinterpret_cast<uint64_t *>(&x);
    OutBinary<uint64_t>(stream, *bits);
}

/*************************/
/***** BaseMortonKey *****/
/*************************/

/*! \fn BaseMortonKey()
 *  \brief constructor */
BaseMortonKey::BaseMortonKey()
{
    InitializeMortonConstants();
}

/*! \fn ~BaseMortonKey()
 *  \brief destructor */
BaseMortonKey::~BaseMortonKey()
{
    ;
}

/*! \fn uint32_t Double2Int(double d)
 *  \brief convert a double on [0, 1) to an unsigned 32 bit integer */
uint32_t BaseMortonKey::Double2Int(double d)
{
    static_assert(sizeof(double) == sizeof(uint64_t), "TypeError: sizeof(double) != 8\n");
    double tmp_double = d * MAXIMUMINTEGER + MAGIC;
    uint64_t tmp_int;
    std::memcpy(&tmp_int, &tmp_double, sizeof(double));
    return static_cast<uint32_t>(tmp_int); // obtained rightmost 32-bit without breaking strict-aliasing rules
}

/*! \fn void InitializeMortonConstants()
 *  \brief initialize constants used in future calculations */
void BaseMortonKey::InitializeMortonConstants()
{
    
    uint128_t one = 1;
    m1 = (one<<64) + 1; /*!< {0...63...0} 1 {0...63...0} 1 */
    m2 = (one<<64) + (one<<32) + 1; // {0...63...0} 1 {0...31...0} 1 {0...31...0} 1
    
    uint128_t x;
    
    x = 0xffffffffUL; // {0...96...0}{1...32...1}
    c1 = (x<<96)+x; // {1...32...1}{0...64...0}{1...32...1}
    upper32mask0 = ~(x<<96); // {0...32...0}{1...96...1}
    
    x = 0xffffUL; // {0...112...0}{1...16...1}
    c2 = (x<<96) + (x<<48) + x; // {0...16...0}{1...16...1}{0...32...0}{1...16...1}{0...32...0}{1...16...1}
    
    x = 0xffUL; // {0...120...0}{1...8...1}
    c3 = (x<<120) + (x<<96) + (x<<72) + (x<<48) + (x<<24) + x; // {{1...8...1}{0...16...0}}x5 {1...8...1}
    
    x = 0xf00fUL; // {0...112...0}{1111000000001111}
    c4 = (x<<120) + (x<<96) + (x<<72) + (x<<48) + (x<<24) + x; // {000011110000}x10 {00001111}
    
    x = 0xc30c3UL; // {11000011000011000011}
    c5  = (x<<120) + (x<<96) + (x<<72) + (x<<48) + (x<<24) + x; // {110000}x21 {11}
    
    x = 0x249249UL; // {001001001001001001001001}
    c6 =  (x<<120) + (x<<96) + (x<<72) + (x<<48) + (x<<24) + x; // {01} {001}x42
    
}

/*! \fn void OutKey(std::ostream stream, morton_key m_key)
 *  \brief output the particle index and its key as octal digits */
void BaseMortonKey::OutKey(std::ostream &stream, morton_key m_key)
{
    auto pindex = static_cast<unsigned int>(m_key>>96);
    stream << "Particle index: " << std::setw(10) << pindex;
    
    stream << "; Morton Key (octal): ";
    for (int i = 1; i != 33; i++) {
        stream << Key8Level(m_key, i);
    }
}

/*! \fn morton_key Dilate3_Int32(int pos)
 *  \brief spread the bits of pos 3 apart: i.e., {1011} becomes {001 000 001 001} */
BaseMortonKey::morton_key BaseMortonKey::Dilate3_Int32(int pos)
{
    morton_key r = pos;
    // you can OutBinary each step to see what is happening
    r = (r*m1) & c1;
    r = (r*m2) & c2;
    r = (r*0x100010001UL) & c3;
    r = (r*0x10101UL)     & c4;
    r = (r*0x00111UL)     & c5;
    r = (r*0x00015UL)     & c6;
    r = r & upper32mask0;
    return r;
}













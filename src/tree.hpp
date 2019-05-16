//
//  tree.hpp
//  PLAN: PLanetesimal ANalyzer
//
//  Created by Rixin Li on 3/11/16.
//  Copyright © 2016 Rixin Li. All rights reserved.
//
//  Some of the following code are based on programs written by Dr. Philip Pinto during his course (ASTR596) in fall 2015.
//  Descriptions and comments are written by Rixin Li.

/*! \file tree.hpp
 *  \brief provide objects related to SmallVec, particle, tree and planetesimal */

#ifndef tree_hpp
#define tree_hpp

#include "global.hpp"

// declare all important classes first

template <int D>
class Planetesimal;                 /*!< data for one planetesimal structure */

template <class T, int D>
class DataSet;              /*!< data set that encapsulates other data classes */


/******************************/
/********** Vtk Part **********/
/******************************/

/*
 * Given a type T and an lvalue expression x, the following two expressions for lvalue references have different syntax but the same semantics:
 *     reinterpret_cast<T&>(x)
 *     *reinterpret_cast<T*>(&(x))
 * Given a type T and an lvalue expression x, the following two expressions for rvalue references have different syntax but the same semantics:
 *     reinterpret_cast<T&&>(x)
 *     static_cast<T&&>(*reinterpret_cast<T*>(&(x)))
 */

/*! \fn template<typename T, typename std::enable_if<(sizeof(T)==2), int>::type = 0> T endian_reverse(T &x)
 *  \brief wrapping endian reverse functions for types containing 2 bytes */
template<typename T, typename std::enable_if<(sizeof(T)==sizeof(uint16_t)), int>::type = 0>
T endian_reverse(T &x) {
    uint16_t tmp;
    std::memcpy(&tmp, &x, sizeof(T));
#ifdef __GNUC__
    tmp = __builtin_bswap16(tmp);
#else /* __GNUC__ */
    tmp = boost::endian::endian_reverse(tmp);
#endif /* __GNUC__ */
    std::memcpy(&x, &tmp, sizeof(T));
    return x;
}
/*! \fn template<typename T, typename std::enable_if<(sizeof(T)==4), int>::type = 0> T endian_reverse(T &x)
 *  \brief wrapping endian reverse functions for types containing 4 bytes */
template<typename T, typename std::enable_if<(sizeof(T)==sizeof(uint32_t)), int>::type = 0>
T endian_reverse(T &x) {
    uint32_t tmp;
    std::memcpy(&tmp, &x, sizeof(T));
#ifdef __GNUC__
    // __bswap_32(x) may static_cast float into uint32_t internally, which we don't want
    tmp = __builtin_bswap32(tmp);
#else /* __GNUC__ */
    tmp = boost::endian::endian_reverse(tmp);
#endif /* __GNUC__ */
    std::memcpy(&x, &tmp, sizeof(T));
    return x;
}

/*! \fn template<typename T, typename std::enable_if<(sizeof(T)==8), int>::type = 0> T endian_reverse(T &x)
 *  \brief wrapping endian reverse functions for types containing 8 bytes */
template<typename T, typename std::enable_if<(sizeof(T)==sizeof(uint64_t)), int>::type = 0>
T endian_reverse(T &x) {
    uint64_t tmp;
    std::memcpy(&tmp, &x, sizeof(T));
#ifdef __GNUC__
    tmp = __builtin_bswap64(tmp);
#else /* __GNUC__ */
    tmp = boost::endian::endian_reverse(tmp);
#endif /* __GNUC__ */
    std::memcpy(&x, &tmp, sizeof(T));
    return x;
}

/*! \class template <class T, int D> VtkDataScalar
 *  \brief general scalar data set
 *  \tparam T type of data
 *  \tparam D dimension of data */
template <class T, int D>
class VtkDataScalar {
private:
    
public:
    /*! \var std::string data_name
     *  \brief data name */
    std::string data_name;
    
    /*! \var std::string data_type
     *  \brief data type */
    std::string data_type;
    
    /*! \var int num_components
     *  \brief number of components, default is 1 */
    int num_components;
    
    /*! \var std::string table_name
     *  \brief table name, default is "default" */
    std::string table_name;
    
    /*! \var std::streampos pos
     *  \brief the position of data in the file */
    std::streampos pos;
    
    /*! \var SmallVec<int, D> num_cells
     *  \brief the number of cells in each dimension
     *  the number of cells is in fact accessible via multi_array; RL: consider deprecate this in the future */
    SmallVec<int, D> num_cells;
    
    /*! \alias using array_type = boost::multi_array<T, D>
     *  \brief typedef an array type */
    using array_type = boost::multi_array<T, D>;
    
    /*! \alias using shape_type = boost::array<typename array_type::index, D>
     *  \brief typedef a shape type used to resize data */
    using shape_type = typename boost::array<typename array_type::index, D>;
    
    /*! \alias using view_type = typename array_type::template array_view<D>::type
     *  \brief typedef a view type to slice data in D dimensions */
    using view_type = typename array_type::template array_view<D>::type;
    
    /*! \alias template <typename std::enable_if<(D > 1), int>::type = 0> using view_r1d_type = typename array_type::template array_view<D-1>::type
     *  \brief typedef a view type to slice data in D-1 dimensions if possible
     *  When you use this type, remember the angle brackets, e.g., VtkDataScalar<float, dim>::view_r1d_type<> */
    template <typename std::enable_if<(D > 1), int>::type = 0>
    using view_r1d_type = typename array_type::template array_view<D-1>::type;
    
    /*! \alias template <typename std::enable_if<(D > 2), int>::type = 0> using view_r2d_type = typename array_type::template array_view<D-2>::type
     *  \brief typedef a view type to slice data in D-2 dimensions if possible 
     *  When you use this type, remember the angle brackets, e.g., VtkDataScalar<float, dim>::view_r2d_type<> */
    template <typename std::enable_if<(D > 2), int>::type = 0>
    using view_r2d_type = typename array_type::template array_view<D-2>::type;
    
    /*! \var array_type data
     *  \brief multi-dimensional array to hold data */
    array_type data;
    
    /*! \var shape_type shape
     *  \brief the shape of this data array */
    shape_type shape;
    
};

/*! \class template <class T, int D> VtkDataVector
 *  \brief general vector data set
 *  \tparam T type of data
 *  \tparam D dimension of data */
template <class T, int D>
class VtkDataVector {
private:
    
public:
    /*! \var std::string data_name
     *  \brief data name */
    std::string data_name;
    
    /*! \var std::string data_type
     *  \brief data type */
    std::string data_type;
    
    /*! \var std::streampos pos
     *  \brief the position of data in the file */
    std::streampos pos;
    
    /*! \var SmallVec<int, D> num_cells
     *  \brief the number of cells in each dimension
     *  the number of cells is in fact accessible via multi_array; RL: consider deprecate this in the future */
    SmallVec<int, D> num_cells;
    
    /*! \alias using array_type = boost::multi_array<T, D+1>
     *  \brief typedef an array type */
    using array_type = boost::multi_array<T, D+1>;
    
    /*! \alias using shape_type = boost::array<typename array_type::index, D+1>
     *  \brief typedef a shape type used to resize data */
    using shape_type = typename boost::array<typename array_type::index, D+1>;
    
    /*! \alias using view_type = array_type::array_view<D+1>::type
     *  \brief typedef a view type to slice data in D dimensions (+1 dimension due to vector) */
    using view_type = typename array_type::template array_view<D+1>::type;
    
    /*! \alias template <typename std::enable_if<(D > 1), int>::type = 0> using view_r1d_type = typename array_type::template array_view<D>::type
     *  \brief typedef a view type to slice data in D-1 dimensions (+1 dimension due to vector)
     *  When you use this type, remember the angle brackets, e.g., VtkDataVector<float, dim>::view_r1d_type<> */
    template <typename std::enable_if<(D > 1), int>::type = 0>
    using view_r1d_type = typename array_type::template array_view<D>::type;
    
    /*! \alias template <typename std::enable_if<(D > 2), int>::type = 0> using view_r2d_type = typename array_type::template array_view<D-1>::type
     *  \brief typedef a view type to slice data in D-2 dimensions (+1 dimension due to vector)
     *  When you use this type, remember the angle brackets, e.g., VtkDataVector<float, dim>::view_r2d_type<> */
    template <typename std::enable_if<(D > 2), int>::type = 0>
    using view_r2d_type = typename array_type::template array_view<D-1>::type;
    
    /*! \var array_type data
     *  \brief multi-dimensional array to hold data */
    array_type data;
    
    /*! \var shape_type shape
     *  \brief the shape of this data array */
    shape_type shape;
    
};


/*! \class template <class T, int D> VtkData
 *  \brief contains data from vtk files
 *  \tparam T type of data
 *  \tparam D dimension of data */
template <class T, int D>
class VtkData {
private:
    
public:
    /*! \var std::string version
     *  \brief version of vtk standard */
    std::string version;
    
    /*! \var std::string header
     *  \breif file header */
    std::string header;
    
    /*! \var std::string file_format
     *  \breif file format, ACSII OR BINARY */
    std::string file_format;
    
    /*! \var std::string dataset_structure
     *  \breif dataset structure, usually STRUCTURED_POINTS in our data
     *  An example:
     *  DATASET STRUCTURED_POINTS
     *  DIMENSIONS n_x n_y n_z
     *  ORIGIN x y z
     *  SPACING s_x s_y s_z
     */
    std::string dataset_structure;
    
    /*! \var double time
     *  \brief current time in simulation */
    double time;
    
    /*! \var SmallVec<int, D> num_cells
     *  \brief number of cells in each dimension */
    SmallVec<int, D> num_cells {SmallVec<int, D>(0)};
    
    /*! \var SmallVec<double, D> origin
     *  \brief coordinate of origin point */
    SmallVec<double, D> origin {SmallVec<double, D>(-0.1)};
    
    /*! \var SmallVec<double, D> spacing
     *  \brief cell spacing, dx/dy/dz... */
    SmallVec<double, D> spacing {SmallVec<double, D>(0.003125)};
    
    /*! \var long num_cell_data
     *  \brief number of CELL_DATA, should be equal to the product of dimensions */
    long num_cell_data;
    
    /*! \var std::map<std::string, VtkDataScalar<T, D>> scalar_data;
     *  \brief mapping data name to data */
    std::map<std::string, VtkDataScalar<T, D>> scalar_data;
    
    /*! \var std::map<std::string, VtkDataScalar<T, D>>::iterator sca_it
     *  \brief an iterator to the elements in scalar_data 
     *  sca_it seems not to be a good name */
    typename std::map<std::string, VtkDataScalar<T, D>>::iterator sca_it;
    
    /*! \var std::map<std::string, VtkDataVector<T, D>> vector_data;
     *  \brief mapping data name to data */
    std::map<std::string, VtkDataVector<T, D>> vector_data;
    
    /*! \var std::map<std::string, VtkDataVector<T, D>>::iterator vec_it
     *  \brief an iterator to the elements in vector_data */
    typename std::map<std::string, VtkDataVector<T, D>>::iterator vec_it;

    /*! \var int shape_changed_flag
     *  \brief set this flag to re-construct the cell_center coordinates */
    int shape_changed_flag;

    /*! \var typename VtkDataVector<T, D>::array_type cell_center;
     *  \brief store the coordinates of cell center points (following the order of [...][z][y][x])
     */
    typename VtkDataVector<T, D>::array_type cell_center;
    
    /*! \fn VtkData()
     *  \brief constructor */
    VtkData() {
        ;
    }
    
    /*! \fn ~VtkData()
     *  \brief destructor */
    ~VtkData() {
        ;
    }

    /*
     * The following four functions are used to iterate through an object satisfying the Boost MultiArray concept (e.g., multi_array or multi_array_view) and apply a function to each element
     */
    
    /*! \fn template<class U, class F> typename std::enable_if<(U::dimensionality==1), void>::type IterateArrayView(U& array, F f)
     *  \brief function specialization to iterate the final dimension (f is specialized to accept one argument) */
    template<class U, class F>
    typename std::enable_if<(U::dimensionality==1), void>::type IterateBoostMultiArrayConcept(U& array, F f) {
        for (auto& element : array) {
            f(element);
        }
    }
    
    /*! \fn template<class U, class F> typename std::enable_if<(U::dimensionality>1), void>::type IterateBoostMultiArrayConcept(U& array, F f)
     *  \brief function to iterate over an Boost MultiArray concept object (f is specialized to accept one argument) */
    template<class U, class F>
    typename std::enable_if<(U::dimensionality>1), void>::type IterateBoostMultiArrayConcept(U& array, F f) {
        for (auto element : array) {
            IterateBoostMultiArrayConcept<decltype(element), F>(element, f);
        }
    }
    
    /*! \fn template<class U, class F, class... Args> typename std::enable_if<(U::dimensionality==1), void>::type IterateBoostMultiArrayConcept(U& array, F f, Args& ...args)
     *  \brief function specialization to iterate the final dimension (f takes multiple arguments) */
    template<class U, class F, class... Args>
    typename std::enable_if<(U::dimensionality==1), void>::type IterateBoostMultiArrayConcept(U& array, F f, Args& ...args) {
        for (auto& element : array) {
            f(element, args...);
        }
    }
    
    /*! \fn template<class U, class F, class... Args> typename std::enable_if<(U::dimensionality>1), void>::type IterateBoostMultiArrayConcept(U& array, F f, Args& ...args)
     *  \brief function to iterate over an Boost MultiArray concept object (f takes multiple arguments) */
    template<class U, class F, class... Args>
    typename std::enable_if<(U::dimensionality>1), void>::type IterateBoostMultiArrayConcept(U& array, F f, Args& ...args) {
        for (auto element : array) {
            IterateBoostMultiArrayConcept<decltype(element), F, Args...>(element, f, args...);
        }
    }
    
    /*
     * The following three blocks are used to extract a multi_array_view from a multi_array without konwing the Dimension until compilation.
     */
    
    
    /*! \functor template<typename RangeArrayType, size_t Dimension> IndicesBuilder
     *  \brief a functor to build indices */
    template<typename RangeArrayType, size_t Dimension>
    struct IndicesBuilder {
        // Recursively invoke the functor for the next lowest dimension and add the next range.
        static auto Build(const RangeArrayType& range) -> decltype(IndicesBuilder<RangeArrayType, Dimension - 1>::Build(range)[range[Dimension - 1]]) {
            return IndicesBuilder<RangeArrayType, Dimension - 1>::Build(range)[range[Dimension - 1]];
        }
    };
    
    /*! \functor template<typename RangeArrayType> IndicesBuilder<RangeArrayType, 1>
     *  \brief a functor specialization to terminate recursion when Dimension=1 */
    template<typename RangeArrayType>
    struct IndicesBuilder<RangeArrayType, 1> {
        /*
         * In C++11, there are two syntax for function declaration:
         *     return-type identifier ( argument-declarations... )
         * and
         *     auto identifier ( argument-declarations... ) -> return_type
         * They are equivalent. But with the later one, you can specify the return_type using decltype(...), where ... is/are only declared in the argument-declarations. For example:
         *     template <typename T1, typename T2>
         *     auto compose(T1 a, T2 b) -> decltype(a + b);
         */
        static auto Build(const RangeArrayType& range) -> decltype(boost::indices[range[0]]) {
            return boost::indices[range[0]];
        }
    };

    /*! \fn template <typename U, size_t Dimension> typename boost::multi_array<U, Dimension>::template array_view<Dimension>::type ExtractSubArrayView(boost::multi_array<U, Dimension>& array, const boost::array<size_t, Dimension>& corner, const boost::array<size_t, Dimension>& subarray_size)
     *  \brief function to extract a view of subarray (corner, subarray_size) from the master array */
    template <typename U, size_t Dimension>
    typename boost::multi_array<U, Dimension>::template array_view<Dimension>::type ExtractSubArrayView(boost::multi_array<U, Dimension>& array, const boost::array<size_t, Dimension>& corner, const boost::array<size_t, Dimension>& subarray_size) {

        using array_type = boost::multi_array<U, Dimension>;
        using range_type = typename array_type::index_range;

        // Build a random-access container with the ranges.
        std::vector<range_type> range;
        for (size_t i = 0; i != Dimension; ++i) {
            range.push_back(range_type(corner[i], corner[i]+subarray_size[i]));
        }

        // Use the helper functor to build the index object.
        auto index = IndicesBuilder<decltype(range), Dimension>::Build(range);

        typename array_type::template array_view<Dimension>::type view = array[index];
        return view;
    }


    /*! \fn void ReadSingleVtkFile(std::vector<std::string>::iterator it)
     *  \brief read cell data from one combined vtk file from all processors */
    void ReadSingleVtkFile(std::vector<std::string>::iterator it) {
        std::ifstream vtk_file;
        std::string tmp_string;
        vtk_file.open(it->c_str(), std::ios::binary);
        
        const size_t D_T_size = D * sizeof(T);
        const size_t one_T_size = sizeof(T);
        
        if (vtk_file.is_open()) {
            // --> getline
            std::getline(vtk_file, version);
            if (version.compare("# vtk DataFile Version 3.0") != 0 && version.compare("# vtk DataFile Version 2.0") != 0) {
                progIO->log_info << "Warning: First line of " << *it << " is " << version << std::endl;
            }
            // --> getline
            std::getline(vtk_file, header);
            if (header.find("CONSERVED") != std::string::npos) {
                size_t time_pos = header.find("time= ");
                // stod() reads to the end of the number, so we can ignore those afterward
                time = std::stod(header.substr(time_pos+6));
                //time = strtod((header.substr(time_pos+6, 12)).c_str(), NULL); // c++0x
            } else {
                progIO->error_message << "Error: Expect CONSERVED, but read: " << header << std::endl;
                progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                exit(3); // cannot open file
            }
            // --> getline
            std::getline(vtk_file, file_format);
            if (file_format.compare("BINARY") != 0) {
                progIO->error_message << "Error: Unsupported file format: " << file_format << std::endl;
                progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                exit(3); // cannot open file
            }
            // --> getline
            std::getline(vtk_file, dataset_structure);
            if (dataset_structure.compare("DATASET STRUCTURED_POINTS") != 0) {
                progIO->error_message << "Error: Unsupported dataset structure: " << dataset_structure << std::endl;
                progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                exit(3); // cannot open file
            }
            // --> getline
            std::getline(vtk_file, tmp_string, ' ');
            if (tmp_string.compare("DIMENSIONS") == 0) {
                std::istringstream iss;
                std::getline(vtk_file, tmp_string);
                iss.str(tmp_string);
                int tmp_index = 0;
                decltype(num_cells) tmp_num_cells = num_cells;
                // iss >> num_cells[0] >> num_cells[1] >> num_cells[2];
                while (!iss.eof()) {
                    iss >> num_cells[tmp_index++]; // making this independent of dimensions
                }

                // We want to store the number of grid cells, not the number of grid cell corners
                while (tmp_index > 0) {
                    num_cells[--tmp_index]--; // note tmp_index starts from D
                }
                if (num_cells != tmp_num_cells) {
                    shape_changed_flag = 1;
                }
            } else {
                progIO->error_message << "Error: No dimensions info" << std::endl;
                progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                exit(3); // cannot open file
            }
            // --> getline
            std::getline(vtk_file, tmp_string, ' ');
            if (tmp_string.compare("ORIGIN") == 0) {
                std::istringstream iss;
                std::getline(vtk_file, tmp_string);
                iss.str(tmp_string);
                int tmp_index = 0;
                while (!iss.eof()) {
                    iss >> origin[tmp_index++] >> std::ws; // eliminate trailing whitespace
                }
            } else {
                progIO->error_message << "Error: No origin info" << std::endl;
                progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                exit(3); // cannot open file
            }
            // --> getline
            std::getline(vtk_file, tmp_string, ' ');
            if (tmp_string.compare("SPACING") == 0) {
                std::istringstream iss;
                std::getline(vtk_file, tmp_string);
                iss.str(tmp_string);
                int tmp_index = 0;
                while (!iss.eof()) {
                    iss >> spacing[tmp_index++] >> std::ws; // eliminate trailing whitespace
                }
            } else {
                progIO->error_message << "Error: No spacing info" << std::endl;
                progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                exit(3); // cannot open file
            }
            // --> getline
            std::getline(vtk_file, tmp_string, ' ');
            if (tmp_string.compare("CELL_DATA") == 0) {
                std::istringstream iss;
                std::getline(vtk_file, tmp_string);
                iss.str(tmp_string);
                iss >> num_cell_data;
                // RL: note that using accumulate() can take advantage of template int D
                SmallVec<long, D> tmp_num_cells = num_cells; // to accommodate the type of num_cell_data
                long tmp_num_cell_data = std::accumulate(tmp_num_cells.data, tmp_num_cells.data+D, 1, std::multiplies<long>());
                if (num_cell_data != tmp_num_cell_data) {
                    progIO->error_message << "Nx*Ny*Nz = " << tmp_num_cell_data << "!= Cell_Data = " << num_cell_data << std::endl;
                    progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                    exit(3); // cannot open file
                }
            } else {
                progIO->error_message << "Error: No info about the number of CELL_DATA" << std::endl;
                progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                exit(3); // cannot open file
            }
            
            progIO->log_info << *it << ", time = " << time << ", " << file_format << ", " << dataset_structure << ", num_cells = " << num_cells
                << ", origin = " << origin << ", spacing = " << spacing << ", CELL_DATA = " << num_cell_data << std::endl;
            progIO->log_info << "data names:";
            
            while (!vtk_file.eof()) {
                // --> getline for SCALARS or VECTORS
                std::getline(vtk_file, tmp_string, ' ');
                if (tmp_string[0] == '\n') {
                    tmp_string.erase(0, 1); // remove possible leading '\n'
                }
                if (tmp_string.compare("SCALARS") == 0) {
                    // --> getline for dataname
                    std::getline(vtk_file, tmp_string, ' ');
                    
                    sca_it = scalar_data.find(tmp_string);
                    if (sca_it == scalar_data.end()) {
                        scalar_data[tmp_string] = VtkDataScalar<T, D>();
                        sca_it = scalar_data.find(tmp_string);
                        sca_it->second.data_name = tmp_string;
                        sca_it->second.num_cells = num_cells;
                        for (int i = 0; i != D; i++) {
                            sca_it->second.shape[i] = num_cells[D-1-i];
                        }
                        //std::copy(num_cells.data, num_cells.data+D, sca_it->second.shape.begin());
                        sca_it->second.data.resize(sca_it->second.shape);
                    } else {
                        if (num_cells != sca_it->second.num_cells) {
                            for (int i = 0; i != D; i++) {
                                sca_it->second.shape[i] = num_cells[D-1-i];
                            }
                            sca_it->second.num_cells = num_cells;
                            sca_it->second.data.resize(sca_it->second.shape);
                        }
                    }
                    
                    progIO->log_info << " | " << sca_it->second.data_name;
                    
                    // --> getline
                    std::getline(vtk_file, tmp_string);
                    size_t ws_pos = tmp_string.find_first_of(' '); // ws = whitespace
                    if (ws_pos != std::string::npos) { // in case of the existence of numComp
                        std::istringstream iss;
                        iss.str(tmp_string);
                        // std::ws is used to extracts as many whitespace chars as possible from the current position
                        // but notice that std::basic_istream objects have the std::skipws flag set as default: this applies a similar effect before the formatted extraction operations
                        //iss >> sca_it->second.data_type >> std::ws >> sca_it->second.num_components;
                        iss >> sca_it->second.data_type >> sca_it->second.num_components;
                    } else {
                        sca_it->second.data_type = tmp_string;
                        sca_it->second.num_components = 1;
                    }
                    if (sca_it->second.data_type.compare("float") != 0) {
                        progIO->error_message << "Error: Expected float format, found " << tmp_string << std::endl;
                        progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                        exit(3); // cannot open file
                    }
                    
                    // --> getline
                    std::getline(vtk_file, tmp_string);
                    if (tmp_string.compare("LOOKUP_TABLE default") != 0) {
                        progIO->error_message << "Error: Expected \"LOOKUP_TABLE default\", unsupportted file" << std::endl;
                        progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                        exit(3); // cannot open file
                    }
                    sca_it->second.table_name = "default";
                    sca_it->second.pos = vtk_file.tellg();
                    
                    std::streampos length = one_T_size * num_cell_data;
                    vtk_file.read(reinterpret_cast<char*>(sca_it->second.data.data()), length);
                    
                    for (auto it = sca_it->second.data.data(); it != sca_it->second.data.data()+sca_it->second.data.num_elements(); it++) {
                        *it = endian_reverse<T>(*it);
                    }
                    
                } else if (tmp_string.compare("VECTORS") == 0) { // if (tmp_string.compare("SCALARS") == 0)
                    // --> getline for dataname
                    std::getline(vtk_file, tmp_string, ' ');
                    
                    vec_it = vector_data.find(tmp_string);
                    if (vec_it == vector_data.end()) {
                        vector_data[tmp_string] = VtkDataVector<T, D>();
                        vec_it = vector_data.find(tmp_string);
                        vec_it->second.data_name = tmp_string;
                        vec_it->second.num_cells = num_cells;
                        for (int i = 0; i != D; i++) {
                            vec_it->second.shape[i] = num_cells[D-1-i];
                        }
                        vec_it->second.shape[D] = D; // vector!
                        vec_it->second.data.resize(vec_it->second.shape);
                    } else {
                        if (num_cells != vec_it->second.num_cells) {
                            for (int i = 0; i != D; i++) {
                                sca_it->second.shape[i] = num_cells[D-1-i];
                            }
                            vec_it->second.shape[D] = D; // vector!
                            vec_it->second.num_cells = num_cells;
                            vec_it->second.data.resize(vec_it->second.shape);
                        }
                    }
                    progIO->log_info << " | " << vec_it->second.data_name;
                    
                    // --> getline
                    std::getline(vtk_file, tmp_string);
                    vec_it->second.data_type = tmp_string;
                    if (vec_it->second.data_type.compare("float") != 0) {
                        progIO->error_message << "Error: Expected float format, found " << tmp_string << std::endl;
                        progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                        exit(3); // cannot open file
                    }
                    
                    vec_it->second.pos = vtk_file.tellg();
                    
                    std::streampos length = D_T_size * num_cell_data;
                    vtk_file.read(reinterpret_cast<char*>(vec_it->second.data.data()), length);
                    for (auto it = vec_it->second.data.data(); it != vec_it->second.data.data()+vec_it->second.data.num_elements(); it++) {
                        *it = endian_reverse<T>(*it);
                    }
                    
                } else { // if (tmp_string.compare("SCALARS") == 0)
                    // it seems vtk file has a new empty line in the end
                    if (tmp_string.length() != 0) {
                        progIO->error_message << "Error: Expected SCALARS or VECTORS, found " << tmp_string << std::endl;
                        progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                        exit(3); // cannot open file
                    }
                }
            }
            progIO->log_info << " |" << std::endl;
            progIO->Output(std::clog, progIO->log_info, __even_more_output, __all_processors);
            vtk_file.close();
        } else { // if (vtk_file.is_open())
            progIO->error_message << "Error: Failed to open file " << it->c_str() << std::endl;
            progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
            exit(3); // cannot open file
        }

    }

    /*! \fn void ReadMultipleVtkFile(std::vector<std::string>::iterator begin, std::vector<std::string>::iterator end)
     *  \brief read cell data from a series of *.vtk file created by each processor
     *  Assume that one cpu core read one snapshot once */
    void ReadMultipleVtkFile(std::vector<std::string>::iterator begin, std::vector<std::string>::iterator end) {
        /*
         * read in all files' data
         * identify NGrid[D] by origin
         * create multi_array views to accept data
         */
        std::ifstream vtk_file;
        std::string tmp_string;
        std::vector<SmallVec<double, D>> origins;   /*!< origin coordinates of each files */
        std::vector<SmallVec<int, D>> dimensions;   /*!< num_cells of each file */
        std::vector<SmallVec<double, D>> endings;   /*!< diagonal corner opposite to the grid origin */
        std::vector<long> grid_cells;               /*!< num_cell_data of each file */
        SmallVec<double, D> ending;                 /*!< diagonal corner opposite to the total mesh origin */
        num_cell_data = 0;

        const size_t D_T_size = D * sizeof(T);
        const size_t one_T_size = sizeof(T);

        vtk_file.open(begin->c_str(), std::ios::binary);
        if (vtk_file.is_open()) {
            // --> getline
            std::getline(vtk_file, version);
            if (version.compare("# vtk DataFile Version 3.0") != 0 &&
                version.compare("# vtk DataFile Version 2.0") != 0) {
                progIO->log_info << "Warning: First line of " << *begin << " is " << version << std::endl;
            }
            // --> getline
            std::getline(vtk_file, header);
            if (header.find("CONSERVED") != std::string::npos) {
                size_t time_pos = header.find("time= ");
                time = std::stod(header.substr(time_pos + 6));
            } else {
                progIO->error_message << "Error in " << *begin << ": Expect CONSERVED, but read: " << header << std::endl;
                progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                exit(3); // cannot open file
            }
            // --> getline
            std::getline(vtk_file, file_format);
            if (file_format.compare("BINARY") != 0) {
                progIO->error_message << "Error in " << *begin << ": Unsupported file format: " << file_format << std::endl;
                progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                exit(3); // cannot open file
            }
            // --> getline
            std::getline(vtk_file, dataset_structure);
            if (dataset_structure.compare("DATASET STRUCTURED_POINTS") != 0) {
                progIO->error_message << "Error in " << *begin << ": Unsupported dataset structure: " << dataset_structure << std::endl;
                progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                exit(3); // cannot open file
            }
            // --> getline
            std::getline(vtk_file, tmp_string, ' ');
            if (tmp_string.compare("DIMENSIONS") == 0) {
                std::istringstream iss;
                std::getline(vtk_file, tmp_string);
                iss.str(tmp_string);
                int tmp_index = 0;
                SmallVec<int, D> tmp_dimensions;
                while (!iss.eof()) {
                    iss >> tmp_dimensions[tmp_index++]; // making this independent of dimensionality
                }
                while (tmp_index > 0) {
                    tmp_dimensions[--tmp_index]--; // note tmp_index starts from D
                }
                dimensions.push_back(tmp_dimensions);
            } else {
                progIO->error_message << "Error in " << *begin << ": No dimensions info" << std::endl;
                progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                exit(3); // cannot open file
            }
            // --> getline
            std::getline(vtk_file, tmp_string, ' ');
            if (tmp_string.compare("ORIGIN") == 0) {
                std::istringstream iss;
                std::getline(vtk_file, tmp_string);
                iss.str(tmp_string);
                int tmp_index = 0;
                SmallVec<double, D> tmp_origin;
                while (!iss.eof()) {
                    iss >> tmp_origin[tmp_index++] >> std::ws; // eliminate trailing whitespace
                }
                origins.push_back(tmp_origin);
            } else {
                progIO->error_message << "Error in " << *begin << ": No origin info" << std::endl;
                progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                exit(3); // cannot open file
            }
            // --> getline
            std::getline(vtk_file, tmp_string, ' ');
            if (tmp_string.compare("SPACING") == 0) {
                std::istringstream iss;
                std::getline(vtk_file, tmp_string);
                iss.str(tmp_string);
                int tmp_index = 0;
                while (!iss.eof()) {
                    iss >> spacing[tmp_index++] >> std::ws; // eliminate trailing whitespace
                }
            } else {
                progIO->error_message << "Error in " << *begin << ": No spacing info" << std::endl;
                progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                exit(3); // cannot open file
            }
            endings.push_back(origins[0] + dimensions[0] * spacing);
            // --> getline
            std::getline(vtk_file, tmp_string, ' ');
            if (tmp_string.compare("CELL_DATA") == 0) {
                std::istringstream iss;
                std::getline(vtk_file, tmp_string);
                iss.str(tmp_string);
                long tmp_num_cell_data;
                iss >> tmp_num_cell_data;
                grid_cells.push_back(tmp_num_cell_data);
                num_cell_data += tmp_num_cell_data;
            } else {
                progIO->error_message << "Error in " << *begin << ": No info about the number of CELL_DATA" << std::endl;
                progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                exit(3); // cannot open file
            }

            progIO->log_info << *begin << ", time = " << time << ", " << file_format << ", " << dataset_structure << ", num_cells = " << dimensions[0]
                             << ", origin = " << origins[0] << ", spacing = " << spacing << ", CELL_DATA = " << num_cell_data << std::endl;
            progIO->log_info << "data names:";

            while (!vtk_file.eof()) {
                // --> getline for SCALARS or VECTORS
                std::getline(vtk_file, tmp_string, ' ');
                if (tmp_string[0] == '\n') {
                    tmp_string.erase(0, 1); // remove possible leading '\n'
                }
                if (tmp_string.compare("SCALARS") == 0) {
                    // --> getline for dataname
                    std::getline(vtk_file, tmp_string, ' ');

                    sca_it = scalar_data.find(tmp_string);
                    if (sca_it == scalar_data.end()) {
                        scalar_data[tmp_string] = VtkDataScalar<T, D>();
                        sca_it = scalar_data.find(tmp_string);
                        sca_it->second.data_name = tmp_string;
                    }
                    progIO->log_info << " | " << sca_it->second.data_name;

                    // --> getline
                    std::getline(vtk_file, tmp_string);
                    size_t ws_pos = tmp_string.find_first_of(' '); // ws = whitespace
                    if (ws_pos != std::string::npos) { // in case of the existence of numComp
                        std::istringstream iss;
                        iss.str(tmp_string);
                        iss >> sca_it->second.data_type >> sca_it->second.num_components;
                    } else {
                        sca_it->second.data_type = tmp_string;
                        sca_it->second.num_components = 1;
                    }
                    if (sca_it->second.data_type.compare("float") != 0) {
                        progIO->error_message << "Error in " << *begin << ": Expected float format, found " << tmp_string << std::endl;
                        progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                        exit(3); // cannot open file
                    }

                    // --> getline
                    std::getline(vtk_file, tmp_string);
                    if (tmp_string.compare("LOOKUP_TABLE default") != 0) {
                        progIO->error_message << "Error in " << *begin << ": Expected \"LOOKUP_TABLE default\", unsupportted file" << std::endl;
                        progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                        exit(3); // cannot open file
                    }
                    sca_it->second.table_name = "default";
                    sca_it->second.pos = vtk_file.tellg();
                    std::streampos length = one_T_size * num_cell_data;
                    vtk_file.seekg(length, vtk_file.cur);

                } else if (tmp_string.compare("VECTORS") == 0) { // if (tmp_string.compare("SCALARS") == 0)
                    // --> getline for dataname
                    std::getline(vtk_file, tmp_string, ' ');

                    vec_it = vector_data.find(tmp_string);
                    if (vec_it == vector_data.end()) {
                        vector_data[tmp_string] = VtkDataVector<T, D>();
                        vec_it = vector_data.find(tmp_string);
                        vec_it->second.data_name = tmp_string;
                    }
                    progIO->log_info << " | " << vec_it->second.data_name;

                    // --> getline
                    std::getline(vtk_file, tmp_string);
                    vec_it->second.data_type = tmp_string;
                    if (vec_it->second.data_type.compare("float") != 0) {
                        progIO->error_message << "Error in " << *begin << ": Expected float format, found " << tmp_string << std::endl;
                        progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                        exit(3); // cannot open file
                    }

                    vec_it->second.pos = vtk_file.tellg();
                    std::streampos length = D_T_size * num_cell_data;
                    vtk_file.seekg(length, vtk_file.cur);

                } else { // if (tmp_string.compare("SCALARS") == 0)
                    // it seems vtk file has a new empty line in the end
                    if (tmp_string.length() != 0) {
                        progIO->error_message << "Error: Expected SCALARS or VECTORS, found " << tmp_string << std::endl;
                        progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                        exit(3); // cannot open file
                    }
                }
            }
            progIO->log_info << " |" << std::endl;
            progIO->Output(std::clog, progIO->log_info, __even_more_output, __all_processors);
            vtk_file.close();
        } else { // if (vtk_file.is_open())
            progIO->error_message << "Error: Failed to open file " << begin->c_str() << std::endl;
            progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
            exit(3); // cannot open file
        }

        for (auto it = (begin+1); it != end; it++) {
            vtk_file.open(it->c_str(), std::ios::binary);
            if (vtk_file.is_open()) {
                std::getline(vtk_file, tmp_string); // version
                std::getline(vtk_file, tmp_string); // header
                std::getline(vtk_file, tmp_string); // file_format
                std::getline(vtk_file, tmp_string); // data_structure
                // --> getline
                std::getline(vtk_file, tmp_string, ' ');
                if (tmp_string.compare("DIMENSIONS") == 0) {
                    std::istringstream iss;
                    std::getline(vtk_file, tmp_string);
                    iss.str(tmp_string);
                    int tmp_index = 0;
                    SmallVec<int, D> tmp_dimensions;
                    while (!iss.eof()) {
                        iss >> tmp_dimensions[tmp_index++];
                    }
                    while (tmp_index > 0) {
                        tmp_dimensions[--tmp_index]--; // note tmp_index starts from D
                    }
                    dimensions.push_back(tmp_dimensions);
                } else {
                    progIO->error_message << "Error in " << *it << ": No dimensions info" << std::endl;
                    progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                    exit(3); // cannot open file
                }
                // --> getline
                std::getline(vtk_file, tmp_string, ' ');
                if (tmp_string.compare("ORIGIN") == 0) {
                    std::istringstream iss;
                    std::getline(vtk_file, tmp_string);
                    iss.str(tmp_string);
                    int tmp_index = 0;
                    SmallVec<double, D> tmp_origin;
                    while (!iss.eof()) {
                        iss >> tmp_origin[tmp_index++] >> std::ws;
                    }
                    origins.push_back(tmp_origin);
                } else {
                    progIO->error_message << "Error in " << *it << ": No origin info" << std::endl;
                    progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                    exit(3); // cannot open file
                }
                endings.push_back(*(--origins.end()) + *(--dimensions.end()) * spacing);
                // RL: assuming spacing is all the same
                std::getline(vtk_file, tmp_string);

                // --> getline
                std::getline(vtk_file, tmp_string, ' ');
                if (tmp_string.compare("CELL_DATA") == 0) {
                    std::istringstream iss;
                    std::getline(vtk_file, tmp_string);
                    iss.str(tmp_string);
                    long tmp_num_cell_data;
                    iss >> tmp_num_cell_data;
                    grid_cells.push_back(tmp_num_cell_data);
                    num_cell_data += tmp_num_cell_data;
                }
                vtk_file.close();
            } else { // if (vtk_file.is_open())
                progIO->error_message << "Error: Failed to open file " << it->c_str() << std::endl;
                progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                exit(3); // cannot open file
            }
        }

        origin = *std::min_element(origins.begin(), origins.end(), SmallVecLessEq<double, double, D>);
        ending = *std::max_element(endings.begin(), endings.end(), SmallVecLessEq<double, double, D>);
        SmallVec<double, D> tmp_double_num_cells = (ending - origin) / spacing;

        decltype(num_cells) tmp_num_cells = num_cells;
        for (int i = 0; i != D; i++) {
            num_cells[i] = static_cast<int>(std::lrint(tmp_double_num_cells[i]));
        }
        if (num_cells != tmp_num_cells) {
            shape_changed_flag = 1;
        }

        // RL: note that using accumulate() can take advantage of template int D
        SmallVec<long, D> tmp_long_num_cells = num_cells; // to accommodate the type of num_cell_data
        long tmp_num_cell_data = std::accumulate(tmp_long_num_cells.data, tmp_long_num_cells.data+D, 1, std::multiplies<long>());
        if (num_cell_data != tmp_num_cell_data) {
            progIO->error_message << "Nx*Ny*Nz = " << tmp_num_cell_data << "!= Cell_Data = " << num_cell_data << std::endl;
            progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
            exit(3); // cannot open file
        }

        for (auto it = scalar_data.begin(); it != scalar_data.end(); it++) {
            if (it->second.num_cells != num_cells) {
                it->second.num_cells = num_cells;
                for (int i = 0; i != D; i++) {
                    it->second.shape[i] = num_cells[D-1-i];
                }
                it->second.data.resize(it->second.shape);
            }
        }

        for (auto it = vector_data.begin(); it != vector_data.end(); it++) {
            if (it->second.num_cells != num_cells) {
                it->second.num_cells = num_cells;
                for (int i = 0; i != D; i++) {
                    it->second.shape[i] = num_cells[D-1-i];
                }
                it->second.shape[D] = D; // vector!
                it->second.data.resize(it->second.shape);
            }
        }
        int file_count = 0;
        for (auto it = begin; it != end; it++) {
            vtk_file.open(it->c_str(), std::ios::binary);
            if (vtk_file.is_open()) {
                std::getline(vtk_file, tmp_string); // version
                std::getline(vtk_file, tmp_string); // header
                std::getline(vtk_file, tmp_string); // file_format
                std::getline(vtk_file, tmp_string); // data_structure
                std::getline(vtk_file, tmp_string); // dimensions
                std::getline(vtk_file, tmp_string); // origin
                std::getline(vtk_file, tmp_string); // spacing
                std::getline(vtk_file, tmp_string); // cell_data
                
                // shared info by scalar and vector
                SmallVec<double, D> tmp_double_corner = (origins[file_count] - origin) / spacing;
                
                boost::array<size_t, D> scalar_corner, scalar_subarray_size;
                boost::array<size_t, D+1> vector_corner, vector_subarray_size;
                boost::multi_array<T, 1> tmp_scalar_data (boost::extents[grid_cells[file_count]]);
                boost::multi_array<T, 1> tmp_vector_data (boost::extents[grid_cells[file_count]*D]);
                for (int i = 0; i != D; i++) {
                    scalar_corner[i] = static_cast<size_t>(std::lrint(tmp_double_corner[D-1-i]));
                    scalar_subarray_size[i] = dimensions[file_count][D-1-i];
                    vector_corner[i] = static_cast<size_t>(std::lrint(tmp_double_corner[D-1-i]));
                    vector_subarray_size[i] = dimensions[file_count][D-1-i];
                }
                vector_corner[D] = 0;
                vector_subarray_size[D] = D;
                
                while (!vtk_file.eof()) {
                    // --> getline for SCALARS or VECTORS
                    std::getline(vtk_file, tmp_string, ' ');
                    if (tmp_string[0] == '\n') {
                        tmp_string.erase(0, 1); // remove possible leading '\n'
                    }
                    if (tmp_string.compare("SCALARS") == 0) {
                        // --> getline for dataname
                        std::getline(vtk_file, tmp_string, ' ');
                        
                        sca_it = scalar_data.find(tmp_string);
                        if (sca_it == scalar_data.end()) {
                            progIO->error_message << "Error in " << *it << ": find unknown data named " << tmp_string << std::endl;
                            progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                            exit(3); // cannot open file
                        }
                        // --> getline for data_type and LOOKUP_TABLE default
                        std::getline(vtk_file, tmp_string);
                        std::getline(vtk_file, tmp_string);
                        
                        std::streamsize length = one_T_size * grid_cells[file_count];
                        vtk_file.read(reinterpret_cast<char*>(tmp_scalar_data.data()), length);
                        
                        auto read_item =  tmp_scalar_data.data();
                        typename VtkDataScalar<T, D>::view_type tmp_view = ExtractSubArrayView(sca_it->second.data, scalar_corner, scalar_subarray_size);
                        
                        IterateBoostMultiArrayConcept(tmp_view, [](T& element, decltype(read_item)& tmp_data)->void {
                            element = endian_reverse<T>(*tmp_data);
                            tmp_data++;
                        }, read_item);
                    } else if (tmp_string.compare("VECTORS") == 0) { // if (tmp_string.compare("SCALARS") == 0)
                        // --> getline for dataname
                        std::getline(vtk_file, tmp_string, ' ');
                        
                        vec_it = vector_data.find(tmp_string);
                        if (vec_it == vector_data.end()) {
                            progIO->error_message << "Error in " << *it << ": find unknown data named " << tmp_string << std::endl;
                            progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                            exit(3); // cannot open file
                        }
                        // --> getline for data_type
                        std::getline(vtk_file, tmp_string);
                        
                        std::streamsize length = D_T_size * grid_cells[file_count];
                        vtk_file.read(reinterpret_cast<char*>(tmp_vector_data.data()), length);
                        
                        auto read_item = tmp_vector_data.data();
                        typename VtkDataVector<T, D>::view_type tmp_view = ExtractSubArrayView(vec_it->second.data, vector_corner, vector_subarray_size);
                        
                        IterateBoostMultiArrayConcept(tmp_view, [](T& element, decltype(read_item)& tmp_data)->void {
                            element = endian_reverse<T>(*tmp_data);
                            tmp_data++;
                        }, read_item);
                    } else { // if (tmp_string.compare("SCALARS") == 0)
                        // it seems vtk file has a new empty line in the end
                        if (tmp_string.length() != 0) {
                            progIO->error_message << "Error: Expected SCALARS or VECTORS, found " << tmp_string << std::endl;
                            progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                            exit(3); // cannot open file
                        }
                    }
                }
                vtk_file.close();
                
            } else { // if (vtk_file.is_open())
                progIO->error_message << "Error: Failed to open file " << it->c_str() << std::endl;
                progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);

            }
            file_count++;
        }

    }

    /*! \fn void ReadVtkFile(int loop_count)
     *  \brief read cell data from *.vtk files */
    void ReadVtkFile(int loop_count) {

        if (progIO->flags.combined_flag) {
            ReadSingleVtkFile(progIO->file_name.vtk_data_file_name.begin()+loop_count);
        } else {
            std::vector<std::string>::iterator file_head = progIO->file_name.vtk_data_file_name.begin();
            ReadMultipleVtkFile(file_head + loop_count * progIO->num_cpus,
                                file_head + loop_count * progIO->num_cpus + progIO->num_cpus);
        }


        // Now construct the cell center
        if (shape_changed_flag) {
            typename VtkDataVector<T, D>::shape_type tmp_shape;
            for (int i = 0; i != D; i++) {
                tmp_shape[i] = num_cells[D-1-i];
            }
            tmp_shape[D] = D; // vector!
            cell_center.resize(tmp_shape);
            
            long element_count = 0;
            // use long type to avoid possible overflow
            SmallVec<long, D> index_limits = num_cells;
            SmallVec<long, D> limit_products;
            SmallVec<long, D> tmp_indices;
            SmallVec<T, D> tmp_cell_center;
            SmallVec<T, D> tmp_origin = origin + spacing * SmallVec<T, D>(0.5);
            size_t tmp_size = D * sizeof(T);
            std::partial_sum(index_limits.data, index_limits.data+D, limit_products.data, std::multiplies<long>());
            for (auto item = cell_center.data(); item != cell_center.data()+cell_center.num_elements(); ) {
                tmp_indices[0] = element_count % limit_products[0];
                for (int dim = D-1; dim != 0; dim--) {
                    long tmp_element_count = element_count; // need to consider all higher dimensions
                    for (int d = D-1; d != dim ; d--) {
                        tmp_element_count -= tmp_indices[d] * limit_products[d-1]; // notice d-1
                    }
                    tmp_indices[dim] = tmp_element_count / limit_products[dim-1]; // notice dim-1
                }
                tmp_cell_center = tmp_origin + spacing * tmp_indices;
                std::memcpy(item, tmp_cell_center.data, tmp_size);
                std::advance(item, D);
                element_count++;
            }
            progIO->log_info << "cell_center constructed, [0][0][0] = (";
            std::copy(cell_center.data(), cell_center.data()+D, std::ostream_iterator<T>(progIO->log_info, ", "));
            progIO->log_info << ")" << std::endl;
            progIO->Output(std::clog, progIO->log_info, __even_more_output, __master_only);

            // reset this flag
            shape_changed_flag = 0;
        } // if (shape_changed_flag)

        progIO->numerical_parameters.cell_length = spacing;
        progIO->numerical_parameters.cell_volume = std::accumulate(spacing.data, spacing.data+D, 1.0, std::multiplies<double>());

        // RL: tried to check time, but vtk and lis have different precisions on time info
        // RL: consider check the numerical_parameters.box_min/box_max with info from vtk files
    }

    
};


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
    // RL: sorting members by alignment reduce the struct size (here 80->72)
    // ref: https://stackoverflow.com/a/119128/4009531

    /*! \var SmallVec<double, D> r, v
     *  \brief position vector r and velocity vector v */
    SmallVec<double, D> pos, vel;

    /*! \var double density
     *  \brief local particle density */
    double density;

    /*! \var uint32_t id_in_run
     *  \brief particle ID in the simulation */
    uint32_t id_in_run;

    /*! \var uint32_t id
     *  \brief particle ID in total particle set (sometimes id_in_run is not contiguous) */
    uint32_t id;

    /*! \var int property_index
     *  \brief index of particle properties */
    int property_index;

    /*! \var uint16_t cpu_id
     *  \brief the processor ID this particle belongs to */
    uint16_t cpu_id;
    
    /*! \fn Particle<D>& operator = (const Particle &rhs)
     *  \brief assignment operator =
    Particle<D>& operator = (const Particle<D> &rhs) {
        pos = rhs.pos;
        vel = rhs.vel;
        property_index = rhs.property_index;
        density  = rhs.density;
        id_in_run = rhs.id_in_run;
        cpu_id = rhs.cpu_id;
        id = rhs.id;
        return *this;
    } */
};

/*! \class template <int D> ParticleSet
 *  \brief data for the entire particle set
 *  \tparam D dimension of data */
template <int D>
class ParticleSet {
private:
    
public:
    /*! \alias using ivec = SmallVec<int, D>
     *  \brief define a vector of int type */
    using ivec = SmallVec<int, D>;
    
    /*! \alias using fvec = SmallVec<float, D>
     *  \brief define a vector of float type */
    using fvec = SmallVec<float, D>;
    
    /*! \alias using dvec = SmallVec<double, D>
     *  \brief define a vector of double type */
    using dvec = SmallVec<double, D>;
    
    /*! \var int num_particles
     *  \brief number of particles */
    uint32_t num_particles;
    
    /*! \var int num_ghost_particles
     *  \brief number of ghost particles */
    uint32_t num_ghost_particles;
    
    /*! \var int num_ghost_particles
     *  \brief number of total particles (must < 2^32-2) */
    uint32_t num_total_particles;
    
    /*! \var int num_types
     *  \brief number of particle types */
    unsigned int num_types;
    
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

    /*! \var double *new_densities
     *  \brief new densities calculated from KNN search */
    double *new_densities;
    
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
        assert(i < num_total_particles);
        return *(particles+i);
    }
    
    /*! \fn Particle<D>& operator[] (const size_t i)
     *  \brief overload operator[] for element modification */
    Particle<D>& operator[] (const size_t i) {
        assert(i < num_total_particles);
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
        particles = new Particle<D>[N];
    }
    
    /*! \fn void ReadMultipleLisFile(std::vector<std::string>::iterator begin, std::vector<std::string>::iterator end)
     *  \brief read particle data from a series of *.lis file created by each processor
     *  Assume that one cpu core read one snapshot once */
    void ReadMultipleLisFile(std::vector<std::string>::iterator begin, std::vector<std::string>::iterator end) {
        std::ifstream lis_file;
        long tmp_num_particles;
        float tmp_coor_lim[12], tmp_float_value, tmp_float_vector[D];
        // First step, obtain the box limit from RootMin and RootMax and count the total particle numbers
        lis_file.open(begin->c_str(), std::ios::binary);
        if (lis_file.is_open()) {
            lis_file.read(reinterpret_cast<char*>(tmp_coor_lim), 12*sizeof(float));
            for (int i = 0; i != 12; i++) {
                coor_lim[i] = static_cast<double>(tmp_coor_lim[i]);
            }
            progIO->log_info << *begin << ", x1l = " << coor_lim[0] << ", x1u = " << coor_lim[1]
            << ", x2l = " << coor_lim[2] << ", x2u = " << coor_lim[3]
            << ", x3l = " << coor_lim[4] << ", x3u = " << coor_lim[5]
            << ", x1dl = " << coor_lim[6] << ", x1du = " << coor_lim[7]
            << ", x2dl = " << coor_lim[8] << ", x2du = " << coor_lim[9]
            << ", x3dl = " << coor_lim[10] << ", x3du = " << coor_lim[11] << "\n";
            lis_file.read(reinterpret_cast<char*>(&num_types), sizeof(int));
            progIO->log_info << "num_types = " << num_types;
            type_info.resize(num_types);
            
            for (unsigned int i = 0; i != num_types; i++) {
                lis_file.read(reinterpret_cast<char*>(&tmp_float_value), sizeof(float));
                type_info[i] = static_cast<double>(tmp_float_value);
                progIO->log_info << ": type_info[" << i << "] = " << type_info[i];
            }
            progIO->log_info << "; || ";
            lis_file.read(reinterpret_cast<char*>(&tmp_float_value), sizeof(float));
            time = static_cast<double>(tmp_float_value);
            lis_file.read(reinterpret_cast<char*>(&tmp_float_value), sizeof(float));
            dt = static_cast<double>(tmp_float_value);
            progIO->log_info << "time = " << time << ", dt = " << dt;
            lis_file.read(reinterpret_cast<char*>(&tmp_num_particles), sizeof(long));
            num_particles = static_cast<uint32_t>(tmp_num_particles);
            lis_file.close();
        } else { // if (lis_file.is_open())
            progIO->error_message << "Error: Failed to open file " << begin->c_str() << std::endl;
            progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
            exit(3); // cannot open file
        }
        
        for (auto it = (begin+1); it != end; it++) {
            lis_file.open(it->c_str(), std::ios::binary);
            if (lis_file.is_open()) {
                lis_file.seekg((14+num_types)*sizeof(float)+sizeof(int), std::ios::beg);
                lis_file.read(reinterpret_cast<char*>(&tmp_num_particles), sizeof(long));
                num_particles += static_cast<uint32_t>(tmp_num_particles);
                lis_file.close();
            } else {
                progIO->error_message << "Error: Failed to open file " << it->c_str() << std::endl;
                progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                exit(3); // cannot open file
            }
        }
        
        num_total_particles = num_particles;
        //uint32_t tmp_num_particles_in_each_processor = num_total_particles / progIO->num_cpus;
        progIO->log_info << ", num_particles = " << num_particles << "; || ";
        
        AllocateSpace(num_particles);
        
        // Thrid step, read particle data
        uint32_t tmp_id = 0;
        unsigned long tmp_long;
        unsigned int tmp_int;
        Particle<D> *p;
        size_t D_float = D * sizeof(float);
        size_t one_float = sizeof(float);
        size_t one_int = sizeof(int);
        size_t one_long = sizeof(long);
        //size_t one_ili = one_int + sizeof(long) + one_int;
        
        for (auto it = begin; it != end; it++) {
            lis_file.open(it->c_str(), std::ios::binary);
            if (lis_file.is_open()) {
                lis_file.seekg((14+num_types)*sizeof(float)+sizeof(int), std::ios::beg);
                lis_file.read(reinterpret_cast<char*>(&tmp_num_particles), sizeof(long));
                std::stringstream content;
                content << lis_file.rdbuf();
                std::string tmp_str = content.str();
                const char *tmp_char = tmp_str.data();
                for (uint32_t i = 0; i != tmp_num_particles; i++) {
                    p = &particles[tmp_id];
                    std::memcpy((char*)tmp_float_vector, tmp_char, D_float);
                    for (int d = 0; d != D; d++) {
                        p->pos[d] = static_cast<double>(tmp_float_vector[d]);
                    }
                    std::advance(tmp_char, D_float);
                    std::memcpy((char*)tmp_float_vector, tmp_char, D_float);
                    for (int d = 0; d != D; d++) {
                        p->vel[d] = static_cast<double>(tmp_float_vector[d]);
                    }
                    std::advance(tmp_char, D_float);
                    std::memcpy((char*)&tmp_float_value, tmp_char, one_float);
                    p->density = static_cast<double>(tmp_float_value);
                    std::advance(tmp_char, one_float);
                    std::memcpy((char*)&p->property_index, tmp_char, one_int);

                    std::advance(tmp_char, one_int);
                    std::memcpy((char*)&tmp_long, tmp_char, one_long);
                    p->id_in_run = static_cast<uint32_t>(tmp_long);
                    std::advance(tmp_char, one_long);
                    std::memcpy((char*)&tmp_int, tmp_char, one_int);
                    p->cpu_id = static_cast<uint16_t>(tmp_int);
                    std::advance(tmp_char, one_int);

                    // RL: we do not use this id anymore, see reasons in ReadSingleLisFile
                    //p->id = p->cpu_id * tmp_num_particles_in_each_processor + p->id_in_run;

                    // RL: previously, we didn't consider the particle ID in the simulations
                    // Fortunately, sampling in Athena only output the first XXX particles in each processor
                    //std::advance(tmp_char, one_ili);
                    //p->id = tmp_id;
                    tmp_id++;
                }
                lis_file.close();
            } else { // if (lis_file.is_open())
                progIO->error_message << "Error: Failed to open file " << it->c_str() << std::endl;
                progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                exit(3); // cannot open file
            }
        }
        
        uint32_t tmp_index = num_particles - 1;
        progIO->log_info << "Last particle's info: id = " << particles[tmp_index].id << ", property_index = " << particles[tmp_index].property_index << ", rad = " << particles[tmp_index].density << ", pos = " << particles[tmp_index].pos << ", v = " << particles[tmp_index].vel << std::endl;
        progIO->Output(std::clog, progIO->log_info, __even_more_output, __all_processors);
        
    }
    
    /*! \fn void ReadSingleLisFile(std::vector<std::string>::iterator it)
     *  \brief read particle data from one combined lis file from all processors */
    void ReadSingleLisFile(std::vector<std::string>::iterator it) {
        std::ifstream lis_file;
        long tmp_num_particles;
        float tmp_coor_lim[12], tmp_float_value, tmp_float_vector[D];
        
        lis_file.open(it->c_str(), std::ios::binary);
        if (lis_file.is_open()) {
            lis_file.read(reinterpret_cast<char*>(tmp_coor_lim), 12*sizeof(float));
            for (int i = 0; i != 12; i++) {
                coor_lim[i] = static_cast<double>(tmp_coor_lim[i]);
            }
            progIO->log_info << *it << ", x1l = " << coor_lim[0] << ", x1u = " << coor_lim[1]
            << ", x2l = " << coor_lim[2] << ", x2u = " << coor_lim[3]
            << ", x3l = " << coor_lim[4] << ", x3u = " << coor_lim[5]
            << ", x1dl = " << coor_lim[6] << ", x1du = " << coor_lim[7]
            << ", x2dl = " << coor_lim[8] << ", x2du = " << coor_lim[9]
            << ", x3dl = " << coor_lim[10] << ", x3du = " << coor_lim[11] << "\n";
            lis_file.read(reinterpret_cast<char*>(&num_types), sizeof(int));
            progIO->log_info << "num_types = " << num_types;
            type_info.resize(num_types);
            for (unsigned int i = 0; i != num_types; i++) {
                lis_file.read(reinterpret_cast<char*>(&tmp_float_value), sizeof(float));
                type_info[i] = static_cast<double>(tmp_float_value);
                progIO->log_info << ": type_info[" << i << "] = " << type_info[i];
            }
            progIO->log_info << "; || ";
            lis_file.read(reinterpret_cast<char*>(&tmp_float_value), sizeof(float));
            time = static_cast<double>(tmp_float_value);
            lis_file.read(reinterpret_cast<char*>(&tmp_float_value), sizeof(float));
            dt = static_cast<double>(tmp_float_value);
            progIO->log_info << "time = " << time << ", dt = " << dt;
            lis_file.read(reinterpret_cast<char*>(&tmp_num_particles), sizeof(long));
            num_particles = static_cast<uint32_t>(tmp_num_particles);
            num_total_particles = num_particles;
            progIO->log_info << ", num_particles = " << num_particles << "; || ";
            
            AllocateSpace(num_particles);
            //uint32_t max_particle_id_in_run = 0;
            // RL: the old way has a poor support for simulations where not all cpus have particles, especially when some cpus have less particles than others (e.g., b/c their decomposed grid intersect with the designated region where Athena initializes particles)
            //uint32_t tmp_num_particles_in_each_processor = num_particles / progIO->num_cpus;
            
            // Third step, read particle data
            uint32_t tmp_id = 0; unsigned long tmp_long; unsigned int tmp_int;
            Particle<D> *p;
            size_t D_float = D * sizeof(float);
            size_t one_float = sizeof(float);
            size_t one_int = sizeof(int);
            size_t one_long = sizeof(long);
            //size_t one_ili = one_int + one_long + one_int;
            
            std::stringstream content;
            content << lis_file.rdbuf();
            std::string tmp_str = content.str();
            const char *tmp_char = tmp_str.data();
            for (uint32_t i = 0; i != tmp_num_particles; i++) {
                p = &particles[tmp_id];
                std::memcpy((char*)tmp_float_vector, tmp_char, D_float);
                for (int d = 0; d != D; d++) {
                    p->pos[d] = static_cast<double>(tmp_float_vector[d]);
                }
                std::advance(tmp_char, D_float);
                std::memcpy((char*)tmp_float_vector, tmp_char, D_float);
                for (int d = 0; d != D; d++) {
                    p->vel[d] = static_cast<double>(tmp_float_vector[d]);
                }
                std::advance(tmp_char, D_float);
                std::memcpy((char*)&tmp_float_value, tmp_char, one_float);
                p->density = static_cast<double>(tmp_float_value);
                std::advance(tmp_char, one_float);
                std::memcpy((char*)&p->property_index, tmp_char, one_int);

                std::advance(tmp_char, one_int);
                std::memcpy((char*)&tmp_long, tmp_char, one_long);
                p->id_in_run = static_cast<uint32_t>(tmp_long);
                std::advance(tmp_char, one_long);
                std::memcpy((char*)&tmp_int, tmp_char, one_int);
                p->cpu_id = static_cast<uint16_t>(tmp_int);
                std::advance(tmp_char, one_int);

                // RL: in this way, we ensure all particles have their own unique id and this method works for tracking particles even if some particles are missing (e.g., out of the box). However, if only cpus with a large cpu_id have particles, we may run into a final particle id larger than 2^32. Thus, this method requires the unique id stored in a long variable. A low-risky workaround may be: counting the number of cpu_id's by setting 1 or 0 to an int array[num_cpus], and then giving an index number with an increasing order to those cpus that have particles initially (e.g., {0, 0, 1, 2, 0, 0, 3, 4, 0, 5, 0, 6, 7, 8}). The final unique id would be array[cpu_id] * max_particle_id_in_run. Still, this workaround does not guarantee id < 2^32.
                //max_particle_id_in_run = std::max(max_particle_id_in_run, p->id_in_run);

                // RL: using 'num_particles / num_cpus' cannot reflect the real and unique id under some special circumstances (see above, before the definition of 'tmp_num_particles_in_each_processor'
                // p->id = p->cpu_id * tmp_num_particles_in_each_processor + p->id_in_run;

                // RL: previously, we didn't consider the particle ID in the simulations
                // Fortunately, sampling in Athena only output the first XXX particles in each processor
                //std::advance(tmp_char, one_ili);
                //p->id = tmp_id;
                tmp_id++;
            }
            lis_file.close();
        } else { // if (lis_file.is_open())
            progIO->error_message << "Error: Failed to open file " << it->c_str() << std::endl;
            progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
            exit(3); // cannot open file
        }
        
        uint32_t tmp_index = num_particles - 1;
        progIO->log_info << "Last particle's info: id = " << particles[tmp_index].id << ", property_index = " << particles[tmp_index].property_index << ", rad = " << particles[tmp_index].density << ", pos = " << particles[tmp_index].pos << ", v = " << particles[tmp_index].vel << std::endl;
        progIO->Output(std::clog, progIO->log_info, __even_more_output, __all_processors);
    }
    
    
    /*! \fn void ReadLisFile(int loop_count)
     *  \brief read particle data from *.lis files */
    void ReadLisFile(int loop_count) {
        auto &paras = progIO->numerical_parameters;

        if (progIO->flags.combined_flag) {
            ReadSingleLisFile(progIO->file_name.lis_data_file_name.begin()+loop_count);
        } else {
            auto file_head = progIO->file_name.lis_data_file_name.begin();
            ReadMultipleLisFile(file_head + loop_count * progIO->num_cpus,
                                file_head + loop_count * progIO->num_cpus + progIO->num_cpus);
        }

        // RL: given that we do not lose particles very often, we sort particles first by cpu_id and then by id_in_run. The resulting index is the unique id assigned to particles for the purpose of tracking particles along time.
        std::sort(particles, particles+num_particles, [](const Particle<D> &a, const Particle<D> &b) {
            if (a.cpu_id == b.cpu_id) {
                return a.id_in_run < b.id_in_run;
            }
            return a.cpu_id < b.cpu_id;
        });
        for (uint32_t i = 0; i != num_particles; i++) {
            particles[i].id = i;
        }
        progIO->physical_quantities[loop_count].time = time;
        progIO->physical_quantities[loop_count].dt = dt;

        if (loop_count == mpi->loop_begin) {
            paras.box_min = SmallVec<double, D>(coor_lim[6], coor_lim[8], coor_lim[10]);
            paras.box_max = SmallVec<double, D>(coor_lim[7], coor_lim[9], coor_lim[11]);
            paras.box_center = (paras.box_min + paras.box_max) / 2.0;
            paras.box_length = paras.box_max - paras.box_min;
            paras.CalculateNewParameters();

            paras.mass_per_particle.resize(num_types);
            paras.mass_fraction_per_species.resize(num_types);
            double tmp_sum = 0;
            for (unsigned int i = 0; i != num_types; i++) {
                paras.mass_fraction_per_species[i] = type_info[i];
                tmp_sum += type_info[i];
            }
            for (unsigned int i = 0; i != num_types; i++) {
                paras.mass_fraction_per_species[i] /= tmp_sum;
                paras.mass_total_code_units = paras.solid_to_gas_ratio * paras.box_length[0] * paras.box_length[1] * std::sqrt(2.*paras.PI);
                paras.mass_per_particle[i] = paras.mass_fraction_per_species[i] * paras.mass_total_code_units / num_particles;
            }
        }

        if (progIO->flags.user_defined_box_flag) {
            for (int i = 0; i != dim; i++) {
                if (progIO->user_box_min[i] == 0 && progIO->user_box_max[i] == 0) {
                    progIO->user_box_min[i] = paras.box_min[i];
                    progIO->user_box_max[i] = paras.box_max[i];
                }
            }
            if (paras.box_min.InRange(progIO->user_box_min, progIO->user_box_max) && paras.box_max.InRange(progIO->user_box_min, progIO->user_box_max)) {
                progIO->log_info << "User-defined coordinate limits are beyond the original box. Nothing to do." << std::endl;
            } else {
                progIO->log_info << "User-defined coordinate limits are in effect: min = " << progIO->user_box_min << "; max = " << progIO->user_box_max << ". Turning on No_Ghost flag is recommended. ";

                Particle<D> *user_selected_particles;
                user_selected_particles = new Particle<D>[num_particles];
                uint32_t num_user_selected_particles = 0;
                for (uint32_t i = 0; i != num_particles; i++) {
                    if (particles[i].pos.InRange(progIO->user_box_min, progIO->user_box_max)) {
                        user_selected_particles[num_user_selected_particles] = particles[i];
                        num_user_selected_particles++;
                    }
                }

                progIO->log_info << num_user_selected_particles << " particles are picked out. ";

                Reset();
                AllocateSpace(num_user_selected_particles);
                std::memcpy(particles, user_selected_particles, sizeof(Particle<D>)*num_user_selected_particles);

                num_particles = num_user_selected_particles;
                num_total_particles = num_particles;

                delete [] user_selected_particles;
                user_selected_particles = nullptr;
            }
            progIO->Output(std::clog, progIO->log_info, __more_output, __all_processors);
        }
    }

    /*! \fn void MakeGhostParticles(const NumericalParameters &paras)
     *  \brief make ghost particles for ghost zone based on ghost zone size
     *  Note this function assumes we are dealing with 3D data. We can implement more if there are other situations. */
    void MakeGhostParticles(NumericalParameters &paras) {
        if (progIO->flags.no_ghost_particle_flag) {
            return;
        }

        // to save the memory usage, we need to check how many ghost particles
        Particle<D> *ghost_particles_x, *ghost_particles_y;
        uint32_t tmp_id = num_particles, ghost_id = 0;
        double new_y = 0;

        dvec non_ghost_width = paras.box_half_width - paras.ghost_zone_width;
        for (int i = 0; i != D; i++) {
            if (non_ghost_width[i] < 0) {
                non_ghost_width[i] = 0;
            }
        }

        dvec non_ghost_min = paras.box_center - non_ghost_width;
        dvec non_ghost_max = paras.box_center + non_ghost_width;

        // First, we make ghost particles for radial direction which need shear mapping
        ghost_id = 0;
        for (uint32_t i = 0; i != num_particles; i++) {
            if (particles[i].pos[0] > non_ghost_max[0] || particles[i].pos[0] < non_ghost_min[0]) {
                ghost_id++;
            }
        }
        ghost_particles_x = new Particle<D>[ghost_id];
        ghost_id = 0;
        // N.B.: f(x, y, z) = f(x + Lx, y - (q Omega Lx) * t, z)
        //       f(x - Lx, y + (q Omega Lx) * t, z) = f(x, y, z)
        for (uint32_t i = 0; i != num_particles; i++) {
            if (particles[i].pos[0] > non_ghost_max[0]) {
                ghost_particles_x[ghost_id] = particles[i];
                ghost_particles_x[ghost_id].id = tmp_id++;
                ghost_particles_x[ghost_id].pos[0] -= paras.box_length[0];
                new_y = ghost_particles_x[ghost_id].pos[1] + paras.shear_speed * time;
                // new_y = new_y [- ymin] - int( (new_y - ymin) / L_Y ) * L_Y [+ ymin]
                ghost_particles_x[ghost_id].pos[1] = new_y - static_cast<int>((new_y - paras.box_min[1]) / paras.box_length[1]) * paras.box_length[1];
                ghost_id++;
            }
            if (particles[i].pos[0] < non_ghost_min[0]) {
                ghost_particles_x[ghost_id] = particles[i];
                ghost_particles_x[ghost_id].id = tmp_id++;
                ghost_particles_x[ghost_id].pos[0] += paras.box_length[0];
                new_y = ghost_particles_x[ghost_id].pos[1] - paras.shear_speed * time;
                // new_y = [ymax -] ( ([ymax -] new_y) + int( (ymax - new_y) / L_Y ) * L_Y )
                ghost_particles_x[ghost_id].pos[1] = new_y + static_cast<int>((paras.box_max[1] - new_y) / paras.box_length[1]) * paras.box_length[1];
                ghost_id++;
            }
        }

        // Second, we make ghost particles for other direction
        // note that ghost particles may also produce ghost particles
        uint32_t tmp_num_ghost_particles = ghost_id;
        ghost_id = 0;
        for (uint32_t i = 0; i != tmp_num_ghost_particles; i++) {
            if (ghost_particles_x[i].pos[1] < non_ghost_min[1] || ghost_particles_x[i].pos[1] > non_ghost_max[1]) {
                ghost_id++;
            }
        }
        for (uint32_t i = 0; i != num_particles; i++) {
            if (particles[i].pos[1] > non_ghost_max[1] || particles[i].pos[1] < non_ghost_min[1]) {
                ghost_id++;
            }
        }
        ghost_particles_y = new Particle<D>[ghost_id];
        ghost_id = 0;
        for (uint32_t i = 0; i != tmp_num_ghost_particles; i++) {
            if (ghost_particles_x[i].pos[1] < non_ghost_min[1]) {
                ghost_particles_y[ghost_id] = ghost_particles_x[i];
                ghost_particles_y[ghost_id].id = tmp_id++;
                ghost_particles_y[ghost_id].pos[1] += paras.box_length[1];
                ghost_id++;
            }
            if (ghost_particles_x[i].pos[1] > non_ghost_max[1]) {
                ghost_particles_y[ghost_id] = ghost_particles_x[i];
                ghost_particles_y[ghost_id].id = tmp_id++;
                ghost_particles_y[ghost_id].pos[1] -= paras.box_length[1];
                ghost_id++;
            }
        }
        for (uint32_t i = 0; i != num_particles; i++) {
            if (particles[i].pos[1] < non_ghost_min[1]) {
                ghost_particles_y[ghost_id] = particles[i];
                ghost_particles_y[ghost_id].id = tmp_id++;
                ghost_particles_y[ghost_id].pos[1] += paras.box_length[1];
                ghost_id++;
            }
            if (particles[i].pos[1] > non_ghost_max[1]) {
                ghost_particles_y[ghost_id] = particles[i];
                ghost_particles_y[ghost_id].id = tmp_id++;
                ghost_particles_y[ghost_id].pos[1] -= paras.box_length[1];
                ghost_id++;
            }
        }

        num_total_particles = tmp_id;
        num_ghost_particles = ghost_id + tmp_num_ghost_particles;
        assert(num_total_particles == num_particles + num_ghost_particles);

        progIO->log_info << "Finish making ghost particles: num_ghost_particles = " << num_ghost_particles << ", and now num_total_particles = " << num_total_particles << std::endl;
        progIO->Output(std::clog, progIO->log_info, __more_output, __all_processors);

        // Third, combine particles and ghost_particles
        Particle<D> *tmp_particles = new Particle<D>[num_particles];
        std::memcpy(tmp_particles, particles, sizeof(Particle<D>)*num_particles);
        AllocateSpace(num_total_particles);
        std::memcpy(particles, tmp_particles, sizeof(Particle<D>)*num_particles);
        // release memory immediately for better memory usage
        delete [] tmp_particles;
        tmp_particles = nullptr;

        std::memcpy(particles+num_particles, ghost_particles_x, sizeof(Particle<D>)*tmp_num_ghost_particles);
        std::memcpy(particles+num_particles+tmp_num_ghost_particles, ghost_particles_y, sizeof(Particle<D>)*ghost_id);

        /* this is a small check for ghost particles
        dvec box_limit = paras.box_max + paras.ghost_zone_width;
        for (uint32_t i = num_ghost_particles; i != num_total_particles; i++) {
            assert (particles[i].pos <= box_limit);
        }
        //*/

        // Four, release memory
        delete [] ghost_particles_x;
        ghost_particles_x = nullptr;
        delete [] ghost_particles_y;
        ghost_particles_y = nullptr;
    }

    /*! \fn void MakeFinerSurfaceDensityMap(const int Nx, const int Ny)
     *  \brief output the solid surface density with finer resolution */
    double** MakeFinerSurfaceDensityMap(const unsigned int Nx, const unsigned int Ny) {
        double **Sigma_ghost = nullptr;
        // clang gives error: variable-sized object may not be initialized
        double tmp_Sigma[Ny+4][Nx+4]; // = {0.0}; // 1 more cell each side as ghost zones
        double ccx[Nx+4], ccy[Ny+4], tmp, idx_origin, idy_origin; // cell-center-x/y
        double dx, dy, inv_dx, inv_dy, dx2, dy2, half_dx, half_dy, three_half_dx, three_half_dy;
        std::vector<double> sigma_per_particle;

        if (progIO->flags.user_defined_box_flag) {
            dx = (progIO->user_box_max[0] - progIO->user_box_min[0]) / Nx;
            dy = (progIO->user_box_max[1] - progIO->user_box_min[1]) / Ny;
        } else {
            dx = progIO->numerical_parameters.box_length[0] / Nx;
            dy = progIO->numerical_parameters.box_length[1] / Ny;
        }

        inv_dx = 1./dx; dx2 = dx*dx; half_dx = dx/2.; three_half_dx = 1.5*dx;
        inv_dy = 1./dy; dy2 = dy*dy; half_dy = dy/2.; three_half_dy = 1.5*dy;

        // usually, dx = dy
        progIO->numerical_parameters.ghost_zone_width = sn::dvec(dx, dy, 0);
        progIO->numerical_parameters.max_ghost_zone_width = std::max(dx, dy);
        MakeGhostParticles(progIO->numerical_parameters);

        if (progIO->flags.user_defined_box_flag) {
            tmp = progIO->user_box_min[0] - 2.5 * dx;
            idx_origin = progIO->user_box_min[0] - dx; // start from cell 1
        } else {
            tmp = progIO->numerical_parameters.box_min[0] - 2.5 * dx; // start from outside
            idx_origin = progIO->numerical_parameters.box_min[0] - dx; // start from cell 1
        }
        std::generate(ccx, ccx+Nx+4, [&tmp, &dx]() {
            tmp += dx;
            return tmp;
        });

        if (progIO->flags.user_defined_box_flag) {
            tmp = progIO->user_box_min[1] - 2.5 * dy;
            idy_origin = progIO->user_box_min[1] - dy; // start from cell 1
        } else {
            tmp = progIO->numerical_parameters.box_min[1] - 2.5 * dy; // start from outside
            idy_origin = progIO->numerical_parameters.box_min[1] - dy; // start from cell 1
        }
        std::generate(ccy, ccy+Ny+4, [&tmp, &dy]() {
            tmp += dy;
            return tmp;
        });

        sigma_per_particle.resize(num_types);
        for (unsigned int i = 0; i != num_types; i++) {
            sigma_per_particle[i] = progIO->numerical_parameters.mass_per_particle[i] / dx / dy;
        }

#ifndef OpenMP_ON
        boost::multi_array<double, 2> Sigma_ghost_private;
        Sigma_ghost_private.resize(boost::extents[Ny+4][Nx+4]);
        Sigma_ghost = new double *[Ny+4];
        Sigma_ghost[0] = Sigma_ghost_private.data();
        for (int i = 1; i != Ny+4; i++) {
            Sigma_ghost[i] = Sigma_ghost[i-1] + Nx+4;
        }
        Particle<D> *p;
        int idx, idy;
        double dist, weightx[3], weighty[3];
#else
        boost::multi_array<double, 3> Sigma_ghost_private;
        Sigma_ghost_private.resize(boost::extents[progIO->numerical_parameters.num_avail_threads][Ny+4][Nx+4]);

        omp_set_num_threads(progIO->numerical_parameters.num_avail_threads);
#pragma omp parallel private(Sigma_ghost)
        {
            int omp_myid = omp_get_thread_num();
            Sigma_ghost = new double *[Ny+4];
            Sigma_ghost[0] = Sigma_ghost_private.data() + omp_myid * (Ny+4) * (Nx+4);
            for (int i = 1; i != Ny + 4; i++) {
                Sigma_ghost[i] = Sigma_ghost[i - 1] + Nx + 4;
            }
            Particle<D> *p;
            int idx, idy;
            double dist, weightx[3], weighty[3];

#pragma omp for
#endif
            for (uint32_t i = 0; i < num_total_particles; i++) {
                p = &particles[i];
                idx = static_cast<int>(std::floor((p->pos[0] - idx_origin) * inv_dx));
                idy = static_cast<int>(std::floor((p->pos[1] - idy_origin) * inv_dy));

                if (progIO->flags.user_defined_box_flag) {
                    if (idx > Nx+1 || idx < 0) {
                        continue;
                    }
                    if (idy > Ny+1 || idy < 0) {
                        continue;
                    }
                } else {
                    if (idx == Nx+2) {
                        idx -= 1; // for exactly surface_max[0]
                    }
                    if (idy == Ny+2) {
                        idy -= 1; // for exactly surface_max[1]
                    }
                    if (idx == -1) {
                        idx = 0; // for exactly surface_min[0]
                    }
                    if (idy == -1) {
                        idy = 0; // for exactly surface_min[1]
                    }
                }

                for (int j = 0; j != 3; j++) {
                    dist = std::abs(p->pos[0] - ccx[idx + j]);
                    if (dist <= half_dx) {
                        weightx[j] = 0.75 - dist * dist / dx2;
                    } else if (dist < three_half_dx) {
                        weightx[j] = 0.5 * std::pow(1.5 - dist / dx, 2.);
                    } else {
                        weightx[j] = 0.;
                    }
                    dist = std::abs(p->pos[1] - ccy[idy + j]);
                    if (dist <= half_dy) {
                        weighty[j] = 0.75 - dist * dist / dy2;
                    } else if (dist < three_half_dy) {
                        weighty[j] = 0.5 * std::pow(1.5 - dist / dy, 2.);
                    } else {
                        weighty[j] = 0.;
                    }
                }

                // RL: check use
                //assert(std::abs(weightx[0]+weightx[1]+weightx[2]-1) > 3e-16);
                //assert(std::abs(weighty[0]+weighty[1]+weighty[2]-1) > 3e-16);

                for (int j = 0; j != 3; j++) {
                    for (int k = 0; k != 3; k++) {
                        Sigma_ghost[idy + j][idx + k] +=
                                sigma_per_particle[p->property_index] * weighty[j] * weightx[k];
                    }
                }
            }
#ifdef OpenMP_ON
        }
#endif
      std::memcpy(tmp_Sigma[0], Sigma_ghost_private.data(), sizeof(double)*(Ny+4)*(Nx+4));
#ifdef OpenMP_ON
        for (unsigned int i = 1; i != progIO->numerical_parameters.num_avail_threads; i++) {
            std::transform(tmp_Sigma[0], &tmp_Sigma[Ny+3][Nx+4], Sigma_ghost_private.data()+i*(Ny+4)*(Nx+4), tmp_Sigma[0], std::plus<double>());
        }
        Sigma_ghost_private.resize(boost::extents[0][0][0]);
#else
        Sigma_ghost_private.resize(boost::extents[0][0]);
#endif
        double **Sigma_p = new double *[Ny];
        Sigma_p[0] = new double[Ny * Nx];
        std::fill(Sigma_p[0], Sigma_p[0] + Ny * Nx, 0.0);
        std::memcpy(Sigma_p[0], &tmp_Sigma[2][2], sizeof(double)*Nx);
        for (int i = 1; i != Ny; i++) {
            Sigma_p[i] = Sigma_p[i - 1] + Nx;
            std::memcpy(Sigma_p[i], &tmp_Sigma[i+2][2], sizeof(double)*Nx);
        }

        /*
        for (int i = 0; i != Ny; i++) {
            Sigma_p[i][0] += tmp_Sigma[i][Nx+1];
            Sigma_p[i][Nx-1] += tmp_Sigma[i][0];
        }
        for (int i = 0; i != Nx; i++) {
            Sigma_p[0][i] += tmp_Sigma[Ny+1][i];
            Sigma_p[Ny-1][i] += tmp_Sigma[0][i];
        }
        Sigma_p[0][0] += tmp_Sigma[Ny+1][Nx+1];
        Sigma_p[Ny-1][Nx-1] += tmp_Sigma[0][0];
        Sigma_p[0][Nx-1] += tmp_Sigma[Ny+1][0];
        Sigma_p[Ny-1][0] += tmp_Sigma[0][Nx+1];
         */

        if (Sigma_ghost != nullptr){
            delete [] Sigma_ghost;
            Sigma_ghost = nullptr;
        }

        return Sigma_p;
    }

    /*! \fn void RebuildVtk(const int Nx, const int Ny, const int Nz, std::string filename)
     *  \brief output the vtk file with particle density and momentum */
    template <class T>
    void RebuildVtk(const unsigned int &Nx, const unsigned int &Ny, const unsigned int &Nz, std::string &filename) {
        auto &paras = progIO->numerical_parameters;
        VtkDataScalar<T, D> rhop;
        VtkDataVector<T, D> w;
        rhop.data.resize(boost::extents[Nz][Ny][Nx]);
        w.data.resize(boost::extents[Nz][Ny][Nx][3]);

        boost::multi_array<T, 3> rhop_ghost;
        boost::multi_array<T, 4> w_ghost;
        rhop_ghost.resize(boost::extents[Nz+4][Ny+4][Nx+4]);
        w_ghost.resize(boost::extents[Nz+4][Ny+4][Nx+4][3]);

        double ccx[Nx+4], ccy[Ny+4], ccz[Nz+4], tmp, idx_origin, idy_origin, idz_origin; // cell-center-x/y
        double dx, dy, dz, inv_dx, inv_dy, inv_dz, dx2, dy2, dz2;
        double half_dx, half_dy, half_dz, three_half_dx, three_half_dy, three_half_dz;
        std::vector<double> rhop_per_particle;

        dx = paras.box_length[0] / Nx;
        dy = paras.box_length[1] / Ny;
        dz = paras.box_length[2] / Nz;

        inv_dx = 1./dx; dx2 = dx*dx; half_dx = dx/2.; three_half_dx = 1.5*dx;
        inv_dy = 1./dy; dy2 = dy*dy; half_dy = dy/2.; three_half_dy = 1.5*dy;
        inv_dz = 1./dz; dz2 = dz*dz; half_dz = dz/2.; three_half_dz = 1.5*dz;

        // usually, dx = dy = dz
        paras.ghost_zone_width = sn::dvec(dx, dy, dz);
        paras.max_ghost_zone_width = MaxOf(dx, dy, dz);
        MakeGhostParticles(paras);

        tmp = paras.box_min[0] - 2.5 * dx; // start from outside
        idx_origin = paras.box_min[0] - dx; // start from cell 1
        std::generate(ccx, ccx+Nx+4, [&tmp, &dx]() {
            tmp += dx;
            return tmp;
        });

        tmp = paras.box_min[1] - 2.5 * dy; // start from outside
        idy_origin = paras.box_min[1] - dy; // start from cell 1
        std::generate(ccy, ccy+Ny+4, [&tmp, &dy]() {
            tmp += dy;
            return tmp;
        });

        tmp = paras.box_min[2] - 2.5 * dz; // start from outside
        idz_origin = paras.box_min[2] - dz; // start from cell 1
        std::generate(ccz, ccz+Nz+4, [&tmp, &dz]() {
            tmp += dz;
            return tmp;
        });

        rhop_per_particle.resize(num_types);
        for (unsigned int i = 0; i != num_types; i++) {
            rhop_per_particle[i] = paras.mass_per_particle[i] / dx / dy / dz;
        }

        Particle<D> *p;
        int idx, idy, idz;
        double dist, weightx[3], weighty[3], weightz[3];
        for (uint32_t i = 0; i < num_total_particles; i++) {
            p = &particles[i];
            idx = static_cast<int>(std::floor((p->pos[0] - idx_origin) * inv_dx));
            idy = static_cast<int>(std::floor((p->pos[1] - idy_origin) * inv_dy));
            idz = static_cast<int>(std::floor((p->pos[2] - idz_origin) * inv_dz));

            if (idx == Nx+2) {
                idx -= 1; // for exactly max[0]
            }
            if (idy == Ny+2) {
                idy -= 1; // for exactly max[1]
            }
            if (idz == Nz+2) {
                idz -= 1; // for exactly max[2]
            }
            if (idx == -1) {
                idx = 0; // for exactly min[0]
            }
            if (idy == -1) {
                idy = 0; // for exactly min[1]
            }
            if (idz == -1) {
                idz = 0; // for exactly min[2]
            }

            for (int j = 0; j != 3; j++) {
                dist = std::abs(p->pos[0] - ccx[idx + j]);
                if (dist <= half_dx) {
                    weightx[j] = 0.75 - dist * dist / dx2;
                } else if (dist < three_half_dx) {
                    weightx[j] = 0.5 * std::pow(1.5 - dist / dx, 2.);
                } else {
                    weightx[j] = 0.;
                }
                dist = std::abs(p->pos[1] - ccy[idy + j]);
                if (dist <= half_dy) {
                    weighty[j] = 0.75 - dist * dist / dy2;
                } else if (dist < three_half_dy) {
                    weighty[j] = 0.5 * std::pow(1.5 - dist / dy, 2.);
                } else {
                    weighty[j] = 0.;
                }
                dist = std::abs(p->pos[2] - ccz[idz + j]);
                if (dist <= half_dz) {
                    weightz[j] = 0.75 - dist * dist / dz2;
                } else if (dist < three_half_dz) {
                    weightz[j] = 0.5 * std::pow(1.5 - dist / dz, 2.);
                } else {
                    weightz[j] = 0.;
                }
            }

            for (int j = 0; j != 3; j++) {
                for (int k = 0; k != 3; k++) {
                    for (int l = 0; l != 3; l++) {
                        double tmp_weight = weightz[j] * weighty[k] * weightx[l];
                        rhop_ghost[idz + j][idy + k][idx + l] += rhop_per_particle[p->property_index] * tmp_weight;
                        for (int d = 0; d != 3; d++) {
                            w_ghost[idz + j][idy + k][idx + l][d] += rhop_per_particle[p->property_index] * tmp_weight * p->vel[d];
                        }
                    }
                }
            }
        }

        for (int iz = 0; iz != Nz; iz++) {
            for (int iy = 0; iy != Ny; iy++) {
                std::memcpy(&(rhop.data[iz][iy][0]), &(rhop_ghost[iz+2][iy+2][2]), sizeof(T)*Nx);
                std::memcpy(&(w.data[iz][iy][0][0]), &(w_ghost[iz+2][iy+2][2][0]), sizeof(T)*Nx*3);
            }
        }

        std::ofstream file_vtk;
        file_vtk.open(filename, std::ios::binary);
        if (file_vtk.is_open()) {
            progIO->out_content << "Writing to " << filename << std::endl;
            progIO->Output(std::cout, progIO->out_content, __normal_output, __all_processors);
            file_vtk << "# vtk DataFile Version 3.0" << std::endl;
            file_vtk << "CONSERVED vars at time= " << std::scientific << std::setprecision(6) << time
                     << ", level= 0, domain= 0" << std::endl;
            file_vtk << "BINARY\nDATASET STRUCTURED_POINTS\nDIMENSIONS " << std::fixed << Nx + 1
                     << " " << Ny + 1 << " " << Nz + 1 << std::endl;
            file_vtk << "ORIGIN" << std::scientific << std::setprecision(6) << std::setw(14) << paras.box_center[0]
                     << std::setw(14) << paras.box_center[1] << std::setw(14) << paras.box_center[2] << std::endl;
            file_vtk << "SPACING" << std::scientific << std::setprecision(6) << std::setw(13) << paras.cell_length[0]
                     << std::setw(13) << paras.cell_length[1] << std::setw(13) << paras.cell_length[2] << std::endl;
            file_vtk << "CELL_DATA " << std::fixed << Nx * Ny * Nz << std::endl;
            file_vtk << "SCALARS particle_density float\nLOOKUP_TABLE default" << std::endl;
            for (auto it = rhop.data.data(); it != rhop.data.data() + rhop.data.num_elements(); it++) {
                *it = endian_reverse<T>(*it);
            }
            file_vtk.write(reinterpret_cast<char *>(rhop.data.data()), sizeof(T) * Nx * Ny * Nz);
            file_vtk << std::endl;
            file_vtk << "VECTORS particle_momentum float" << std::endl;
            for (auto it = w.data.data(); it != w.data.data() + w.data.num_elements(); it++) {
                *it = endian_reverse<T>(*it);
            }
            file_vtk.write(reinterpret_cast<char *>(w.data.data()), D * sizeof(T) * Nx * Ny * Nz);
            file_vtk.close();
        } else {
            progIO->error_message << "Failed to open " << filename << std::endl;
            progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
        }

        rhop.data.resize(boost::extents[0][0][0]);
        w.data.resize(boost::extents[0][0][0][0]);
        rhop_ghost.resize(boost::extents[0][0][0]);
        w_ghost.resize(boost::extents[0][0][0][0]);

    }

};


/***********************************/
/********** MortonKey Part *********/
/***********************************/

/*
 * A MortonKey is a 128-bit (16 byte) number whose leftmost 32 bits are the particle index; the remaining 96 bits are the three integer coordinates interleaved into a Morton key.
 * N.B.: modern C++ offers bitset. If you are dealing with super large amout of particles (larger than 2^32 ~ 4.3 billion), then you might want to switch uint128_t to bitset.
 */


/*! \class template <int D> struct Orthant
 *  \brief orthant is the generalization of quadrant and octant and even hyperoctant in n>3 dimensions. This struct gives the directions in each dimension (+1 or -1). We will need them while building trees.
 *  \tparam D dimension of this vector */
template <int D>
struct Orthant {
    static const SmallVec<int, D> orthants[1U<<D];
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
public:
    /*! \alias using uint128_t = __uint128_t
     *  \brief stick to the C++ standard's integer format uint[N]_t */
#ifdef __GNUC__
    using uint128_t = __uint128_t;
#endif /* __GNUC __ */

private:
    /*
     * Below are useful constants worked for dilate3_32
     */

    /*! \var uint128_t m1
     *  \brief binary: {0...63...0} 1 {0...63...0} 1 */
    uint128_t m1;

    /*! \var uint128_t m2
     *  \brief binary: {0...63...0} 1 {0...31...0} 1 {0...31...0} 1 */
    uint128_t m2;

    /*! \var uint128_t c1
     *  \brief binary: {1...32...1}{0...64...0}{1...32...1} */
    uint128_t c1;

    /*! \var uint128_t c2
     *  \brief binary: {0...16...0} {{1...16...1}{0...32...0}}x2 {1...16...1} */
    uint128_t c2;

    /*! \var uint128_t c3
     *  \brief binary: {{1...8...1}{0...16...0}}x5 {1...8...1} */
    uint128_t c3;

    /*! \var uint128_t c4
     *  \brief binary: {000011110000}x10 {00001111} */
    uint128_t c4;

    /*! \var uint128_t c5
     *  \brief binary: {110000}x21 {11} */
    uint128_t c5;

    /*! \var uint128_t c6
     *  \brief binary: {01} {001}x42 */
    uint128_t c6;

    /*! \var uint128_t upper32mask0
     *  \brief binary {0...32...0}{1...96...1} */
    uint128_t upper32mask0;

    /*
     * Magic numbers for double2int: note that the double 0x1p+0(=1 in decimal) cannot be converted this way, so that the range of numbers is, strictly,  [0, 1). This two magic number make sure that the only 2^32 particles can be distinguished in one dimension. In other word, the minimum distance between two particles (or to be recognized as two particles) in one dimension case is 1.0/2^32 = 2.3283e-10.
     * N.B.: if the box is quite large, keep in mind that the resolution before rescaling is max(box_length)/2^32. Thus, if particles only concentrate at a certain region, it is better to focus on that region instead of the entire box.
     */

    /*! \var static constexpr double MAGIC = 6755399441055744.0
     *  \brief MAGIC = 2^52 + 2^51 = (0x0018000000000000)_16 */
    static constexpr double MAGIC = 6755399441055744.0;
    /*! \var static constexpr double MAXIMUMINTEGER = 4294967294.0
     *  \brief MAGIC = 2^32 - 2 = (0x00000000FFFFFFFE)_16 */
    static constexpr double MAXIMUMINTEGER = 4294967294.0;

public:
    // define name alias for intuitive definition
    /*! \alias using morton_key = uint128_t
     *  \brief define a type equivalent to uint128_t for morton key */
    using morton_key = uint128_t;

    /*! \fn BaseMortonKey()
     *  \brief constructor */
    BaseMortonKey();

    /*! \fn ~BaseMortonKey()
     *  \brief destructor */
    ~BaseMortonKey();

    /*! \fn uint32_t Double2Int(double d)
     *  \brief convert a double on [0, 1) to an unsigned 32 bit integer */
    uint32_t Double2Int(double d);

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
 *  \brief interface of morton key */
template <int D>
class MortonKey : public BaseMortonKey {
private:

public:
    /*! \alias using ivec = SmallVec<int, D>
     *  \brief define a vector of int type */
    using ivec = SmallVec<int, D>;

    /*! \alias using fvec = SmallVec<float, D>
     *  \brief define a vector of float type */
    using fvec = SmallVec<float, D>;

    /*! \alias using dvec = SmallVec<double, D>
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

    /*! \fn template <class U> Morton(const SmallVec<U, D> &pos, int index)
     *  \brief convert a position vector pos and particle index into a 128-bit Morton Key */
    template <class U>
    morton_key Morton(const SmallVec<U, D> &pos, uint32_t index) {
        dvec pos_scaled = pos - boxmin;
        for (int d = 0; d != D; d++) {
            pos_scaled[d] *= scale[d];
        }

        SmallVec<uint32_t, D> int_pos;
        for (int d = 0; d != D; d++) {
            int_pos[d] = Double2Int(pos_scaled[d]);
        }
        return Morton(int_pos, index);
    }

    /*! \fn Morton(const SmallVec<uint32_t, D> &pos, int index)
     *  \brief overloading Morton above for uint32_t */
    morton_key Morton(const SmallVec<uint32_t, D> &pos, uint32_t index) {
        morton_key result = (static_cast<uint128_t>(index))<<96;
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
    /*! \alias using ivec = SmallVec<int, D>
     *  \brief define a vector of int type */
    using ivec = SmallVec<int, D>;
    
    /*! \alias using fvec = SmallVec<float, D>
     *  \brief define a vector of float type */
    using fvec = SmallVec<float, D>;
    
    /*! \alias using dvec = SmallVec<double, D>
     *  \brief define a vector of double type */
    using dvec = SmallVec<double, D>;
    
    /*! \struct InternalParticle
     *  \brief necessary particle data for tree */
    struct InternalParticle {
        /*! \var dvec pos
         *  \brief the coordinates of particle position */
        dvec pos;

        /*! \var dvec vel
         *  \brief particle velocity */
        dvec vel;
        
        /*! \var double mass
         *  \brief the particle mass */
        double mass;

        /*! \var double ath_density
         *  \brief solid density of the local cell in ATHENA
         *  we may not need this since we can store them elsewhere and then ask for it when needed */
        double ath_density {0.0};

        /*! \var double new_density
         *  \brief new density calculated by KNN search */
        double new_density {0.0};
        
        /*! \var uint32_t id
         *  \brief the original particle index */
        uint32_t original_id;

        /*! \var uint32_t densest_neighbor_id
         *  \brief the original index of the neighbor with the highest density
         *  remove this variable for better memory usage
        uint32_t densest_neighbor_id {0};  */

        /*! \var uint32_t peak_id
         *  \brief the index of the particle at the peak of the planetesimal that this particle belongs to
         *  remove this variable for better memory usage
        uint32_t peak_id {0}; */

        /*! \var bool sink_flag
         *  \brief if this is a sink particle, then "True"; else "False" */
        bool sink_flag {false};

        /*! \var bool in_clump_flag
         *  \brief if this particle already belongs to a clump, then "True"; else "False" */
        bool in_clump_flag {false};
    };

    /*! \var std::unordered_map<uint32_t, std::vector<InternalParticle>> sink_particle_indices
     *  \brief record sink particles' indices
     *  In order to use the entire info (pos, vel, etc.) in the future, we map the leading particle's index in tree.particle_list with all the particles' indices in particle_set (including the leading one itself). We may need them for more accurate calculation. */
    std::unordered_map<uint32_t, std::vector<InternalParticle>> sink_particle_indices;

    /*! \struct TreeNode
     *  \brief tree node structure */
    struct TreeNode {
        /*! \var dvec center
         *  \brief center coordinates of a node */
        dvec center;
        
        /*! \var double half_width
         *  \brief half the width of a node */
        double half_width;
        
        /*! \var uint32_t begin
         *  \brief the beginning particle index in node */
        uint32_t begin;
        
        /*! \var uint32_t end
         *  \brief the ending particle index in node, notice that this "end" follows the C++ tradition and serves as the off-the-end iterator */
        uint32_t end;
        
        /*! \var uint32_t parent
         *  \brief the parent node's index */
        uint32_t parent;
        
        /*! \var uint32_t first_daughter;
         *  \brief the index of first daughter node */
        uint32_t first_daughter;
        
        /*! \var uint16_t orthant
         *  \brief orthant is the daughter direction from parent */
        uint16_t orthant;
        
        /*! \var uint8_t num_daughters
         *  \brief the number of daughters */
        uint8_t num_daughters;
        
        /*! \var uint8_t level
         *  \brief level in tree */
        uint8_t level;
    };
    
    /*! \var static const int max_level = 32
     *  \brief the maximum levels of this tree is 32 */
    static const int max_level = 32;
    
    /*! \var typename MortonKey<D>::morton_key *morton
     *  \brief store all the morton key of particles */
    typename MortonKey<D>::morton_key *morton;
    
    /*! \var int num_particles
     *  \brief number of particles (must < 2^32-2)
     *  This must < 2^32-1-1 = 0xffffffff-1, because (*TreeNode)->end means the off-the-end iterator, so we need one more number than the total number of particles. Anyway, if we are dealing with more particles than that, we should adapt our tools and use more advanced Morton Keys or hire more processors for one data set */
    uint32_t num_particles;

    /*! \var int enough_particle_resolution_flag
     *  \brief indicating if the number of particles is large enough in terms of grid horizontal resolution */
    int enough_particle_resolution_flag {0};
    
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
    int max_leaf_size {1U<<D};
    
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
    double half_width {0};

    /*! \var double epsilon
     *  \brief tolerance that a particle be outside a node
     *  The tolerance should be root_half_width / 2^32 */
    double epsilon {1. / 4294967296.};
    
    /*! \var int level_count[max_level]
     *  \brief count how many nodes in each level */
    int level_count[max_level];
    
    /*! \var int level_ptr[max_level]
     *  \brief used for building tree, point to currently empty node at each level
     *  The number of nodes in each level has been fixed after CountNodesLeaves() but before FillTree() */
    int level_ptr[max_level];
    
    /*! \var std::vector<std::pair<int, double>> heaps;
     *  \brief a heap stores k-nearest neighbours */
    std::vector<std::pair<int, double>> heaps;
    
    /*! \var const double to_diagonal;
     *  \brief const used in SphereNodeIntersect */
    const double to_diagonal = std::sqrt(D);

    dvec cell_length;
    dvec half_cell_length;
    dvec cell_length_squared;
    dvec three_half_cell_length;
    
    /*! \fn BHtree()
     *  \brief constructor, about the member initializer lists, refer to http://en.cppreference.com/w/cpp/language/initializer_list */
    BHtree() : morton(nullptr), particle_list(nullptr), tree(nullptr), leaf_nodes(nullptr), node2leaf(nullptr) {
        max_daughters = (1U<<D);
        root = 0;
        root_level = 1;
    }
    
    /*! \fn Reset()
     *  \brief release memory and reset */
    void Reset() {
        num_particles = 0;
        
        /*
         * Delete operator releases the memory of objects allocated by new operator; delete [] --> new []. Since delete operator will perform nullptr check first, we can safely use it without check. An expression with the delete[] operator, first calls the appropriate destructors for each element in the array (if these are of a class type), and then calls an array deallocation function (refer to http://www.cplusplus.com/reference/new/operator%20delete[]/).
         * With -std=c++11 and above, we should try to use shared_ptr and unique_ptr for smart memory managements. Mark here and implement later.
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

        sink_particle_indices.clear();
        assert(sink_particle_indices.size() == 0);
    }
    
    /*! \fn ~BHtree()
     *  \brief destructor */
    ~BHtree() {
        Reset();
        ClearHeaps();
    }

    /*! \fn void MakeSinkParticle(ParticleSet<D> &particle_set)
     *  \brief make a sink particle based on all the original particles' info (combine mass, vel, etc.) */
    void  MakeSinkParticle(ParticleSet<D> &particle_set) {
        for (auto it : sink_particle_indices) {
            double mass = 0.0;
            dvec vel(0);
            for (auto it_par : it.second) {
                //double tmp_mass = progIO->numerical_parameters.mass_per_particle[particle_set[it_par].property_index];
                vel += it_par.mass * it_par.vel;
                mass += it_par.mass;
            }
            vel /= mass;
            particle_list[it.first].mass = mass;
            particle_list[it.first].vel = vel;
        }
    }

    /*! \fn void SortPoints()
     *  \brief sort points by morton key and then copy back to particle list */
    void SortPoints() {
        InternalParticle *tmp = new InternalParticle[num_particles];
        typename MortonKey<D>::morton_key *tmp_morton = new typename MortonKey<D>::morton_key[num_particles];

        uint32_t new_index = 0; // used for tmp[...]
        for (uint32_t i = 0; i < num_particles; ) {
            // check if it belongs to a sink particle (at the same position by machine precision)
            tmp_morton[new_index] = morton[i];
            if (i < num_particles && (morton[i]<<32) == (morton[i+1]<<32)) {
                // if so, record its siblings
                tmp[new_index] = particle_list[this->ParticleIndex(morton[i])];
                tmp[new_index].sink_flag = true;
                auto it = sink_particle_indices.emplace(new_index, std::vector<InternalParticle>());
                assert(it.second);
                it.first->second.push_back(particle_list[this->ParticleIndex(morton[i])]);
                it.first->second.push_back(particle_list[this->ParticleIndex(morton[i+1])]);
                uint32_t tmp_i = i + 2;
                while (tmp_i < num_particles && (morton[i]<<32) == (morton[tmp_i]<<32)) {
                    it.first->second.push_back(particle_list[this->ParticleIndex(morton[tmp_i])]);
                    tmp_i++;
                }
                i = tmp_i;
            } else {
                // if not, move forward as normal
                tmp[new_index] = particle_list[this->ParticleIndex(morton[i])];
                i++;
            }
            new_index++;
        }
        delete [] particle_list;
        particle_list = nullptr;
        delete [] morton;
        morton = nullptr;

        num_particles = new_index;
        particle_list = new InternalParticle[num_particles]; // new_index should be the new amount
        morton = new typename MortonKey<D>::morton_key[num_particles];
        std::memcpy(particle_list, tmp, sizeof(InternalParticle)*num_particles);
        std::memcpy(morton, tmp_morton, sizeof(typename MortonKey<D>::morton_key)*num_particles);
        delete [] tmp;
        tmp = nullptr;
        delete [] tmp_morton;
        tmp_morton = nullptr;
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
            
            //if (__level >= max_level) std::cerr << "Hit max_level" << std::endl;
            level_count[__level]++;
            num_nodes++;
            
            if (count <= max_leaf_size || __level == max_level - 1) { // prevent using further levels
                num_leaf_nodes++; // only nodes with leaves < max_leaf_size are leaves
            } else {
                CountNodesLeaves(__level+1, __begin-count, __begin);
            }
            
            if (__begin < __end) {
                orthant = this->Key8Level(morton[__begin], __level); // search next data
            }
        }
    }
    
    /*! \fn void FillTree(const int level, int __begin, const int __end, const int parent, const dvec &__center, const double __half_width)
     *  \brief fill the tree with data */
    void FillTree(const int __level, int __begin, const int __end, const int __parent, const dvec &__center, const double __half_width) { // double underscore here is to avoid confusion with TreeNode member or InternalParticle member
        assert(__level < max_level);
        assert(__end > __begin); // note if this will cause bug
        assert(tree[__parent].first_daughter == 0);
        assert(tree[__parent].num_daughters == 0); // parent shouldn't have any daughters
        
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
            
            // get daughter node number from currently-empty node budget in tree
            int daughter = level_ptr[__level];
            level_ptr[__level]++;
            
            if (tree[__parent].first_daughter == 0) {
                // first daughter
                assert(tree[__parent].num_daughters == 0);
                tree[__parent].first_daughter = daughter;
                tree[__parent].num_daughters = 1;
            } else {
                // subsequent daughters
                tree[__parent].num_daughters++;
                assert(tree[__parent].num_daughters <= max_daughters);
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
            p->num_daughters = 0;
            node_ptr++;
            assert(node_ptr < num_nodes);
            
            if (count <= max_leaf_size || __level == max_level - 1) { // prevent using further levels
                // node with <= max_leaf_size particles is a leaf
                leaf_nodes[num_leaf_nodes] = daughter;
                node2leaf[daughter] = num_leaf_nodes;
                num_leaf_nodes++;
            } else {
                // node with > max_leaf_size particles is a branch
                FillTree(p->level, __begin-count, __begin, daughter, p->center, 0.5*__half_width);
            }
            
            // now next daughter of this node
            if (__begin < __end) {
                orthant = this->Key8Level(morton[__begin], __level);
            }
        }
    }
    
    /*! \fn void BuildTree(NumericalParameters &paras, ParticleSet<D> &particle_set, bool quiet=false, bool check_pos=true)
     *  \brief build tree from particle data */
    void BuildTree(NumericalParameters &paras, ParticleSet<D> &particle_set, bool quiet=false, bool check_pos=true) { // double underscore here is to avoid confusion with all sorts of members
        Reset();

        // RL: if particles only occupy a small region in the box, then this tree does not need to cover the entire box (ghost particles also accounted for).
        if (check_pos) {
            paras.particle_max = particle_set[0].pos;
            paras.particle_min = particle_set[0].pos;
            for (uint32_t i = 1; i != particle_set.num_total_particles; i++) {
                if (!SmallVecLessEq(paras.particle_min, particle_set[i].pos)) {
                    for (int d = 0; d != D; d++) {
                        paras.particle_min[d] = std::min(particle_set[i].pos[d], paras.particle_min[d]);
                    }
                }
                if (!SmallVecGreatEq(paras.particle_max, particle_set[i].pos)) {
                    for (int d = 0; d != D; d++) {
                        paras.particle_max[d] = std::max(particle_set[i].pos[d], paras.particle_max[d]);
                    }
                }
            }
            paras.max_particle_extent = 0.0;
            for (int d = 0; d != D; d++) {
                paras.max_particle_extent = std::max(paras.max_particle_extent, paras.particle_max[d] - paras.particle_min[d]);
            }

            half_width = paras.max_half_width + paras.max_ghost_zone_width;
            root_center = paras.box_center;
            if (paras.max_particle_extent < half_width) {
                // Round it for better analyses
                half_width = (std::ceil(paras.max_particle_extent * 100.) + 1) / 100.;
                for (int d = 0; d != D; d++) {
                    if (paras.particle_max[d] - paras.particle_min[d] < 0.5 * paras.box_length[d]) {
                        root_center[d] = (paras.particle_max[d] + paras.particle_min[d]) / 2.0;
                        root_center[d] = std::round(root_center[d] * 1000.) / 1000.;
                    } else {
                        root_center[d] = paras.box_center[d];
                    }
                }
                progIO->log_info << "According to the spatial distribution of particles, the root node of BH tree is now centered at " << root_center << ", and with a half_width of " << half_width << std::endl;
            }
        } else {
            // if root_center and half_width are not set, we use the default box info
            if (half_width < 1e-16) {
                half_width = paras.max_half_width + paras.max_ghost_zone_width;
                root_center = paras.box_center;
            }
        }

        // Since one processor deals with one entire snapshot, so we need to make sure the total number of particles is smaller than MAX(uint32_t). If needed in the future, we may find a way to deploy multiple processors to deal with the data in one single snapshot.
        assert(particle_set.num_total_particles < (uint32_t)0xffffffff);
        num_particles = particle_set.num_total_particles;
        particle_list = new InternalParticle[num_particles];
        
        for (uint32_t i = 0; i != num_particles; i++) {
            particle_list[i].pos = particle_set[i].pos;
            particle_list[i].vel = particle_set[i].vel;
            particle_list[i].mass = progIO->numerical_parameters.mass_per_particle[particle_set[i].property_index];
            particle_list[i].ath_density = particle_set[i].density;
            // now original id means the particle's permanent id
            particle_list[i].original_id = particle_set[i].id;
            //particle_list[i].original_id = i; // != particle_set[i].id, it is the original index of particles, after sorting by morton key, the order changes
        }
        particle_set.Reset(); // release memory
        // RL: Keep in mind that SortPoints() below and MakeGhostParticles() in ParticleSet all need twice the memory (when copying particle data to a new array); unless they are fixed, we don't need to migrate all the functionality of ParticleSet to BHtree

        // compute Morton Keys and sort particle_list by Morton order
        this->InitMortonKey(root_center-dvec(half_width), root_center+dvec(half_width));
        epsilon = 2 * half_width / std::pow(2, 31); // the index is 31 since the width is only half the box
        morton = new typename MortonKey<D>::morton_key[num_particles];
        
        for (uint32_t i = 0; i != num_particles; i++) {
            morton[i] = this->Morton(particle_list[i].pos, i);
        }
        std::sort(&(morton[0]), &(morton[num_particles]), AscendingMorton());
        
        // RL: I already found particles that are too close (in fact, even at the same position with the machine precision) in my data. There are two ways to deal with:
        // 1), tolerate it for now since max_leaf_size > 1 and disable the uniqueness check of morton keys below
        // 2), combine them into one sink particle and record any useful info
        // now we choose the second approach

        //for (uint32_t i = 0; i != num_particles-1; i++) {
        //    assert((morton[i]<<32) <= (morton[i+1]<<32));
        //}
        SortPoints();
        MakeSinkParticle(particle_set);

        for (uint32_t i = 0; i != num_particles-1; i++) {
            assert((morton[i]<<32) < (morton[i+1]<<32));
        }

        num_leaf_nodes = 0;
        level_count[0] = 1;
        for (int level = 1; level != max_level; level++) {
            level_count[level] = 0;
        }
        
        // run through the data once to determine space required for tree
        num_nodes = 1; // root contribute one
        CountNodesLeaves(root_level, 0, num_particles);
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
        p->num_daughters = 0;
        p->level = root_level;
        p->center = root_center;
        p->half_width = half_width;
        p->begin = 0;
        p->end = num_particles;
        p->parent = std::numeric_limits<uint32_t>::max(); // 4294967295
        
        // run through the data again to build the tree
        num_leaf_nodes = 0;
        FillTree(root_level, 0, num_particles, root, root_center, 0.5*half_width);
        assert(node_ptr + 1 == num_nodes);
        delete [] morton;
        morton = nullptr;

        CheckTree(root, root_level, root_center, half_width);

        // calculate members used to speed up density-kernel-calculations
        cell_length = progIO->numerical_parameters.cell_length;
        half_cell_length = cell_length / 2.;
        cell_length_squared = cell_length * cell_length;
        three_half_cell_length = half_cell_length * 3.;

        if (!quiet) {
            progIO->log_info << "Finish building a tree: num_nodes = " << num_nodes << ", num_sink_particles = " << sink_particle_indices.size();
            if (sink_particle_indices.size() > 0) {
                progIO->log_info << ", the largest sink particle contains " << std::max_element(sink_particle_indices.begin(), sink_particle_indices.end(), [](const typename decltype(sink_particle_indices)::value_type &a, const typename decltype(sink_particle_indices)::value_type &b) {
                    return a.second.size() < b.second.size();
                })->second.size() << " super particles. ";
            }
            progIO->log_info << std::endl;
            progIO->Output(std::clog, progIO->log_info, __more_output, __all_processors);
        }
    }
    
    /*! \fn bool Within(const dvec &__pos, const dvec &node_center, const double __half_width)
     *  \brief determine if a particle is within certain distance of a node center */
    bool Within(const dvec &__pos, const dvec &node_center, const double __half_width) {
        for (int d = 0; d != D; d++) {
            if ( !(__pos[d] >= node_center[d] - __half_width - epsilon && __pos[d] <= node_center[d] + __half_width + epsilon)) {
                return false;
            }
        }
        return true;
    }
    
    /*! \fn void CheckTree(const int node, const int __level, const dvec &node_center, const double __half_width)
     *  \brief traverse the tree and check whether each point is within the node it is supposed to be. Also check levels/widths/centers
     *  \todo if some particles are outside their nodes due to the lack of Morton Key spatial resolution, we need to relocate these particles into their nodes, as long as the number is not big. */
    void CheckTree(const int node, const int __level, const dvec &node_center, const double __half_width) {
        assert(tree[node].level == __level);
        assert(tree[node].half_width == __half_width);
        
        for (uint32_t p = tree[node].begin; p != tree[node].end; p++) {
            if (!Within(particle_list[p].pos, node_center, __half_width)) {
                progIO->error_message << "Particle " << particle_list[p].pos << " outside node " << node_center << " with width " << 2*tree[node].half_width << std::endl;
                progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
            }
        }
        
        for (uint32_t daughter = tree[node].first_daughter; daughter != tree[node].first_daughter + tree[node].num_daughters; daughter++) {
            dvec tmp_center = node_center;
            for (int d = 0; d != D; d++) {
                tmp_center[d] += 0.5 * __half_width * Orthant<D>::orthants[tree[daughter].orthant][d];
            }
            dvec daughter_center = tree[daughter].center;
            assert(tmp_center == daughter_center);
            CheckTree(daughter, __level+1, daughter_center, 0.5*__half_width);
        }
    }
    
    /*! \fn inline bool IsLeaf(const int node)
     *  \breif determine if a node is leaf */
    inline bool IsLeaf(const int node) {
        assert(node < num_nodes);
        return (tree[node].num_daughters == 0); // N.B., TreeNode has begin and end
    }
    
    /*! \fn inline int NodeSize(const int node)
     *  \brief return the number of particles in a node */
    inline unsigned int NodeSize(const int node) {
        assert(node < num_nodes);
        return tree[node].end - tree[node].begin; // note the end is off-the-end iterator
    }
    
    /*! \fn inline bool InNode(const dvec &__pos, const int node)
     *  \brief determine if a particle is in a node */
    inline bool InNode(const dvec &__pos, const int node) {
        return Within(__pos, tree[node].center, tree[node].half_width);
    }
    
    /*! \fn uint32_t Key2Leaf(typename MortonKey<D>::morton_key const __morton, const int node, const int __level)
     *  \brief given a Morton Key, find the node number of the leaf cell where this key belongs */
    uint32_t Key2Leaf(typename MortonKey<D>::morton_key const __morton, const int node, const int __level) {
        // if a leaf, just return answer
        if (IsLeaf(node)) {
            return node;
        }
        
        // else recurse into the correct direction
        int orthant = this->Key8Level(__morton, __level);
        int daughter = -1;
        for (uint32_t d = tree[node].first_daughter; d != tree[node].first_daughter + tree[node].num_daughters; d++) {
            if (tree[d].orthant == orthant) {
                daughter = Key2Leaf(__morton, d, __level+1);
                break;
            }
        }
        if (daughter == -1) {
            progIO->error_message << "Key2Leaf: leaf cell doesn't exist in tree." << std::endl;
            progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
            assert(daughter >= 0);
        }
        return daughter;
    }
    
    /*! \fn uint32_t Pos2Node(const dvec &__pos)
     *  \brief given position, find the index of node containing it */
    uint32_t Pos2Node(const dvec &__pos) {
        assert( Within(__pos, root_center, half_width));
        typename MortonKey<D>::morton_key __morton = this->Morton(__pos, 0); // index doesn't matter here, just give 0
        return Key2Leaf(__morton, root, root_level);
    }
    
    /*! \fn bool SphereContainNode(const dvec &__center, const double r, const int node)
     *  \brief return true if the entire node is within the sphere (__center, r) */
    bool SphereContainNode(const dvec &__center, const double r, const int node) {
        assert(node < num_nodes);
        SmallVec<double, D> max_distance = (tree[node].center - __center);
        double tmp_distance = tree[node].half_width * to_diagonal; // assume node is a cube
        
        for (int i = 0; i != D; i++) {
            max_distance[i] = std::fabs(max_distance[i]) + tmp_distance;
        }
        
        if (max_distance.Norm2() < r * r) {
            return true;
        } else {
            return false;
        }
    }
    
    /*! \fn bool SphereNodeIntersect(const dvec &__center, const double r, const int node)
     *  \brief return true if any part of node is within the sphere (__center, r) */
    bool SphereNodeIntersect(const dvec &__center, const double r, const int node) {
        assert(node < num_nodes);
        double c2c = (tree[node].center - __center).Norm2();
        double tmp_distance = tree[node].half_width * to_diagonal + r;
        
        // check if node is outside the sphere
        if (c2c > tmp_distance * tmp_distance) {
            return false;
        }
        
        // check if node center is inside the sphere or the spheric center is inside the node
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
    
    /*! \functor template <typename T1, typename T2> struct less_second
     *  \brief served as comparison method for heaps */
    template <typename T1, typename T2>
    struct less_second {
        typedef std::pair<T1, T2> type;
        bool operator ()(type const& a, type const& b) const {
            return a.second < b.second;
        }
    };
    
    /*! \fn void Add2Heaps(const unsigned int knn, const int i, const double dr2)
     *  \brief add element to heaps */
    void Add2Heaps(const unsigned int knn, const int i, const double dr2) {
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
    
    /*! \fn void RecursiveKNN(const dvec &__pos, const int node, const double dist, const int knn)
     *  \brief do recursive KNN search to traverse the tree */
    void RecursiveKNN(const dvec &__pos, const int node, const double dist, const int knn) {
        if (SphereContainNode(__pos, dist, node)) {
            for (uint32_t p = tree[node].begin; p != tree[node].end; p++) {
                Add2Heaps(knn, p, (__pos-particle_list[p].pos).Norm2());
            }
        } else if (SphereNodeIntersect(__pos, dist, node)) {
            if (IsLeaf(node)) {
                for (uint32_t p = tree[node].begin; p != tree[node].end; p++) {
                    Add2Heaps(knn, p, (__pos-particle_list[p].pos).Norm2());
                }
            } else {
                for (uint32_t d = tree[node].first_daughter; d != tree[node].first_daughter + tree[node].num_daughters; d++) {
                    RecursiveKNN(__pos, d, dist, knn);
                }
            }
        }
    }
    
    /*! \fn void KNN_Search(const dvec &__pos, const unsigned int knn, double &radius_knn, uint32_t *indices, bool in_order = false)
     *  \brief given position, perform k-nearest neighbours search and return radius and particle indices
     *  Notice that a search around one particle will include that particle in the results */
    void KNN_Search(const dvec &__pos, const unsigned int knn, double &radius_knn, uint32_t *indices, bool in_order = false) {
        assert(knn <= num_particles);
        
        if (heaps.size() != 0) {
            ClearHeaps(); // clear memory first
        }
        heaps.reserve(knn);
        
        double max_dr2 = 0;
        // obtain a estimated distance firstly
        if (Within(__pos, root_center, half_width)) {
            // if within the box
            int node = Pos2Node(__pos);
            while (NodeSize(node) < knn/4) {
                node = tree[node].parent;
            }
            //for (uint32_t p = tree[node].begin; p != tree[node].end; p++) {
            //    max_dr2 = std::max(max_dr2, (__pos - particle_list[p].pos).Norm2());
            //}
            max_dr2 = tree[node].half_width*tree[node].half_width;
        } else {
            // if point is outside the entire box, use the min distance to box boundary
            for (int d = 0; d < D; d++) {
                double dx = MaxOf(root_center[d]-half_width - __pos[d], 0.0, __pos[d] - root_center[d] - half_width);
                max_dr2 += dx * dx;
            }
        }
        
        double max_dr = std::sqrt(max_dr2);
        // now use the distance guess to proceed
        do {
            ClearHeaps();
            heaps.reserve(knn);
            RecursiveKNN(__pos, root, max_dr, knn);
            max_dr *= 2;
        } while (heaps.size() < knn);
        
        // Phil did TraverseKNN once more, but I don't see the necessity...

        radius_knn = std::sqrt(heaps.front().second);
        // if needed, we can sort them by distances
        if (in_order) {
            std::sort_heap(heaps.begin(), heaps.end(), less_second<int, double>());
        }

        // get original particle id, which turned out to be more complicated
        for (unsigned int i = 0; i != heaps.size(); i++) {
            indices[i] = heaps[i].first; // now get new particle id in tree data structure
        }
    }


    /*! \fn void Add2Heaps_OpenMP(const unsigned int knn, const int i, const double dr2, std::vector<std::pair<int, double>> &local_heaps)
     *  \brief add element to heaps */
    void Add2Heaps_OpenMP(const unsigned int knn, const int i, const double dr2, std::vector<std::pair<int, double>> &local_heaps) {
        if (local_heaps.size() < knn) {
            local_heaps.push_back(std::pair<int, double>(i, dr2));
            std::push_heap(local_heaps.begin(), local_heaps.end(), less_second<int, double>());
        } else {
            if (dr2 < local_heaps.front().second) {
                std::pop_heap(local_heaps.begin(), local_heaps.end(), less_second<int, double>());
                local_heaps.pop_back();
                local_heaps.push_back(std::pair<int, double>(i, dr2));
                std::push_heap(local_heaps.begin(), local_heaps.end(), less_second<int, double>());
            }
        }
    }

    /*! \fn void RecursiveKNN_OpenMP(const dvec &__pos, const int node, const double dist, const int knn,  std::vector<std::pair<int, double>> &local_heaps)
     *  \brief do recursive KNN search to traverse the tree */
    void RecursiveKNN_OpenMP(const dvec &__pos, const int node, const double dist, const int knn, std::vector<std::pair<int, double>> &local_heaps) {
        if (SphereContainNode(__pos, dist, node)) {
            for (uint32_t p = tree[node].begin; p != tree[node].end; p++) {
                Add2Heaps_OpenMP(knn, p, (__pos-particle_list[p].pos).Norm2(), local_heaps);
            }
        } else if (SphereNodeIntersect(__pos, dist, node)) {
            if (IsLeaf(node)) {
                for (uint32_t p = tree[node].begin; p != tree[node].end; p++) {
                    Add2Heaps_OpenMP(knn, p, (__pos-particle_list[p].pos).Norm2(), local_heaps);
                }
            } else {
                for (uint32_t d = tree[node].first_daughter; d != tree[node].first_daughter + tree[node].num_daughters; d++) {
                    RecursiveKNN_OpenMP(__pos, d, dist, knn, local_heaps);
                }
            }
        }
    }

    /*! \fn void KNN_SearchWithOpenMP(const dvec &__pos, const unsigned int knn, double &radius_knn, uint32_t *indices,  std::vector<std::pair<int, double>> &local_heaps)
     *  \brief given position, perform k-nearest neighbours search and return radius and particle indices
     *  Notice that a search around one particle will include that particle in the results */
    void KNN_Search_OpenMP(const dvec &__pos, const unsigned int knn, double &radius_knn, uint32_t *indices,  std::vector<std::pair<int, double>> &local_heaps, double estimated_distance=0) {
        assert(knn <= num_particles);

        if (local_heaps.size() != 0) {
            std::vector<std::pair<int, double>>().swap(local_heaps);
        }
        local_heaps.reserve(knn);

        double max_dr;
        if (estimated_distance == 0) {
            double max_dr2 = 0;
            // obtain a estimated distance firstly
            if (Within(__pos, root_center, half_width)) {
                // if within the box
                int node = Pos2Node(__pos);
                while (NodeSize(node) < knn/4) {
                    node = tree[node].parent;
                }
                //for (uint32_t p = tree[node].begin; p != tree[node].end; p++) {
                //    max_dr2 = std::max(max_dr2, (__pos - particle_list[p].pos).Norm2());
                //}
                max_dr2 = tree[node].half_width*tree[node].half_width;
            } else {
                // if point is outside the entire box, use the min distance to box boundary
                for (int d = 0; d < D; d++) {
                    double dx = MaxOf(root_center[d]-half_width - __pos[d], 0.0, __pos[d] - root_center[d] - half_width);
                    max_dr2 += dx * dx;
                }
            }

            max_dr = std::sqrt(max_dr2);
        } else {
            max_dr = estimated_distance;
        }

        // now use the distance guess to proceed
        do {
            std::vector<std::pair<int, double>>().swap(local_heaps);
            local_heaps.reserve(knn);
            RecursiveKNN_OpenMP(__pos, root, max_dr, knn, local_heaps);
            max_dr *= 2;
        } while (local_heaps.size() < knn);

        // Phil did TraverseKNN once more, but I don't see the necessity...

        radius_knn = std::sqrt(local_heaps.front().second);
        // sort them by distances
        std::sort_heap(local_heaps.begin(), local_heaps.end(), less_second<int, double>());

        // get original particle id, which turned out to be more complicated
        for (unsigned int i = 0; i != local_heaps.size(); i++) {
            indices[i] = local_heaps[i].first; // now get new particle id in tree data structure
        }
    }
    
    /*! \fn void RecursiveBallSearchCount(const dvec __pos, int node, const double radius, uint32_t &count)
     *  \brief do recursive ball search (only count numbers) to traverse the tree */
    void RecursiveBallSearchCount(const dvec &__pos, int node, const double radius, uint32_t &count) {
        if (SphereContainNode(__pos, radius, node)) {
            count += (tree[node].end - tree[node].begin);
        } else if (SphereNodeIntersect(__pos, radius, node)) {
            if (IsLeaf(node)) {
                for (uint32_t p = tree[node].begin; p != tree[node].end; p++) {
                    if ((__pos - particle_list[p].pos).Norm2() <= radius * radius) {
                        count++;
                    }
                }
            } else {
                for (uint32_t d = tree[node].first_daughter; d != tree[node].first_daughter + tree[node].num_daughters; d++) {
                    RecursiveBallSearchCount(__pos, d, radius, count);
                }
            }
        }
    }
    
    /*! \fn void RecursiveBallSearch(const dvec &__pos, int node, const double radius, uint32_t *indices, uint32_t &count)
     *  \brief do recursive ball search to traverse the tree */
    void RecursiveBallSearch(const dvec &__pos, int node, const double radius, uint32_t *indices, uint32_t &count) {
        if (SphereContainNode(__pos, radius, node)) {
            for (uint32_t p = tree[node].begin; p != tree[node].end; p++) {
                indices[count++] = p;
            }
        } else if (SphereNodeIntersect(__pos, radius, node)) {
            if (IsLeaf(node)) {
                for (uint32_t p = tree[node].begin; p != tree[node].end; p++) {
                    if ((__pos - particle_list[p].pos).Norm2() <= radius * radius) {
                        indices[count++] = p;
                    }
                }
            } else {
                for (uint32_t d = tree[node].first_daughter; d != tree[node].first_daughter + tree[node].num_daughters; d++) {
                    RecursiveBallSearch(__pos, d, radius, indices, count);
                }
            }
        }
    }
    
    /*! \fn void BallSearch(const dvec &__center, const double radius, uint32_t *indices, uint32_t &count)
     *  \brief perform a search within a sphere */
    void BallSearch (const dvec &__center, const double radius, uint32_t *indices, uint32_t &count) {
        count = 0;
        RecursiveBallSearch(__center, root, radius, indices, count);
        // convert morton key to original particle id (optional step, only useful while interacting with outsiders)
        //for (int i = 0; i < count; i++) {
        //    indices[i] = particle_list[indices[i]].original_id; // make sure what you want
        //}
    }
    
    
    /*! \fn dvec ShearedCenter2Center(const dvec &__center, const int node, dvec &max_distance, double shear_distance)
     *  \brief return the distance between spherical center and sheared node center */
    dvec ShearedCenter2Center(const dvec &__center, const int node, dvec &max_distance, double shear_distance) {
        assert(node < num_nodes);
        dvec tmp_node_center = tree[node].center;
        dvec c2c;
        for (int i = 0; i != D; i++) {
            c2c[i] = __center[i] - tmp_node_center[i];
            if (std::fabs(c2c[i]) > max_distance[i]) {
                if (c2c[i] > 0) {
                    c2c[i] -= progIO->numerical_parameters.box_length[i];
                } else {
                    c2c[i] += progIO->numerical_parameters.box_length[i];
                }
                if (i == 0 && D > 1) {
                    if (__center[0] < 0) { // RL: assume the center of simulation is dvec(0)
                        tmp_node_center[1] -= shear_distance;
                        tmp_node_center[1] += static_cast<int>((progIO->numerical_parameters.box_max[1]-tmp_node_center[1]) / progIO->numerical_parameters.box_length[1]) * progIO->numerical_parameters.box_length[1];
                    } else {
                        tmp_node_center[1] += shear_distance;
                        tmp_node_center[1] -= static_cast<int>((tmp_node_center[1]-progIO->numerical_parameters.box_min[1]) / progIO->numerical_parameters.box_length[1]) * progIO->numerical_parameters.box_length[1];
                    }
                }
            }
        }
        return c2c;
    }
    
    /*! \fn bool SphereContainNodeWithShear(const dvec &__center, const double r, const int node, const dvec &c2c)
     *  \brief return true if the entire node is within the sphere (__center, r) */
    bool SphereContainNodeWithShear(const dvec &__center, const double r, const int node, const dvec &c2c) {
        double tmp_distance = tree[node].half_width * to_diagonal; // assume node is a cube
        dvec tmp_dvec;
        for (int i = 0; i != D; i++) {
            tmp_dvec[i] = std::fabs(c2c[i])+tmp_distance;
        }
        if (tmp_dvec.Norm2() < r * r) {
            return true;
        } else {
            return false;
        }
    }
    
    /*! \fn bool SphereNodeIntersectWithShear(const dvec __center, const double r, const int node, const dvec &c2c)
     *  \brief return true if any part of node is within the sphere (__center, r) */
    bool SphereNodeIntersectWithShear(const dvec __center, const double r, const int node, const dvec &c2c) {
        double c2c_distance = c2c.Norm2();
        double tmp_distance = tree[node].half_width * to_diagonal + r;
        
        // check if node is outside the sphere
        if (c2c_distance > tmp_distance * tmp_distance) {
            return false;
        }
        
        // check if node center is inside the sphere
        if (c2c_distance < r || c2c_distance < tree[node].half_width) {
            return true;
        }
        
        // now do exact check for intersection
        // notice that when we extend each side, the space is divided into multiple sub-space, and the value for each sub-space is different
        dvec pos_min = __center - c2c - dvec(tree[node].half_width);
        dvec pos_max = __center - c2c + dvec(tree[node].half_width);
        
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
    
    /*! \fn void RecursiveBallSearchCountWithShear(const dvec __pos, int node, const double radius, uint32_t &count, dvec &max_distance, double shear_distance)
     *  \brief do recursive ball search (only count numbers) to traverse the tree */
    void RecursiveBallSearchCountWithShear(const dvec __pos, int node, const double radius, uint32_t &count, dvec &max_distance, double shear_distance) {
        dvec c2c = ShearedCenter2Center(__pos, node, max_distance, shear_distance);
        if (SphereContainNodeWithShear(__pos, radius, node, c2c)) {
            count += (tree[node].end - tree[node].begin);
        } else if (SphereNodeIntersectWithShear(__pos, radius, node, c2c)) {
            if (IsLeaf(node)) {
                for (uint32_t p = tree[node].begin; p != tree[node].end; p++) {
                    if ((__pos - particle_list[p].pos).Norm2() <= radius * radius) {
                        count++;
                    }
                }
            } else {
                for (uint32_t d = tree[node].first_daughter; d != tree[node].first_daughter + tree[node].num_daughters; d++) {
                    RecursiveBallSearchCountWithShear(__pos, d, radius, count, max_distance, shear_distance);
                }
            }
        }
    }

    /*! \fn double QuadraticSpline(dvec dist) const
     *  \brief calculate spatial weight using (1D Quadratic Spline kernel)^D */
    double QuadraticSpline(dvec dist) const {
        dist.AbsSelf();
        double weight = 1;
        for (int i = 0; i != D; i++) {
            if (dist[i] <= half_cell_length[i]) {
                weight *= (0.75 - dist[i]*dist[i]/cell_length_squared[i]);
            } else if (dist[i] < three_half_cell_length[i]) {
                weight *= 0.5 * std::pow(1.5 - dist[i]/cell_length[i], 2.);
            } else {
                return 0.;
            }
        }
        return weight;
    }

    /*! \fn template <class T> static double QseudoQuadraticSplinesKernel(const DataSet<T, D> &ds, const double radius, uint32_t self_id, uint32_t *indices, std::vector<std::pair<int, double>> &local_heaps)
     *  \brief use a (1D Quadratic Spline kernel)^D for density calculation
     *  This function considers only 1D Quadratic Spline kernel functions, which are then combined multiplicatively to form multidimensional weights (ref: Appendix A in Youdin & Johansen 2007). It will converts the mass of each particle to a density in the volume of one cell at first, and then sum them weighted by distances.
     *  RL: this kernel works better for small particle sub-sample (N_par ~ N_x * N_y).
     *  For the complete version of particle data, the number of particles is so large that the mass of one single particle is too small. Therefore, the density of any particle becomes very small and cannot exceed 160 even if you calculate 64 particles stacking with each other. What we need (for the full particle data) is a kernel that is independent of particle sampling and cell size. */
    template <class T>
    static double QseudoQuadraticSplinesKernel(const DataSet<T, D> &ds, const double radius, uint32_t self_id, uint32_t *indices, std::vector<std::pair<int, double>> &local_heaps) {
        double density = 0;
        for (unsigned int i = 0; i != progIO->numerical_parameters.num_neighbors_in_knn_search; i++) {
            // static member function, only have access to static data members, so use ds.tree
            density += ds.tree.particle_list[indices[i]].mass / progIO->numerical_parameters.cell_volume * ds.tree.QuadraticSpline(ds.tree.particle_list[indices[i]].pos-ds.tree.particle_list[self_id].pos);
        }
        return density;
    }

    /*! \fn template <class T> static double VerboseQseudoQuadraticSplinesKernel(const DataSet<T, D> &ds, const double radius, uint32_t self_id, uint32_t *indices, std::vector<std::pair<int, double>> &local_heaps)
     *  \brief use a (1D Quadratic Spline kernel)^D for density calculation and output the contribution from each particle
     *  See QseudoQuadraticSplinesKernel() for more details */
    template <class T>
    static double VerboseQseudoQuadraticSplinesKernel(const DataSet<T, D> &ds, const double radius, uint32_t self_id, uint32_t *indices, std::vector<std::pair<int, double>> &local_heaps) {
        double density = 0, tmp_density = 0, QSdist = 0;
        SmallVec<double, D> dist;
        for (unsigned int i = 0; i != progIO->numerical_parameters.num_neighbors_in_knn_search; i++) {
            // static member function, only have access to static data members, so use ds.tree
            dist = ds.tree.particle_list[indices[i]].pos-ds.tree.particle_list[self_id].pos;
            QSdist = ds.tree.QuadraticSpline(dist);
            tmp_density = ds.tree.particle_list[indices[i]].mass / progIO->numerical_parameters.cell_volume * QSdist;

            std::cout << "i=" << i << "dist=" << dist.Norm() << ", QSdist=" << QSdist << ", tmp_d=" << tmp_density << std::endl;

            density += tmp_density;
        }
        return density;
    }

    /*! \fn template <class T> static double PureSphericalDensityKernel(const DataSet<T, D> &ds, const double radius, uint32_t self_id, uint32_t *indices, std::vector<std::pair<int, double>> &local_heaps)
     *  \brief calculate density by dividing the total mass of certain neighbors by a spherical volume (radius = the distance to the furthest neighbor) */
    template <class T>
    static double PureSphericalDensityKernel(const DataSet<T, D> &ds, const double radius, uint32_t self_id, uint32_t *indices, std::vector<std::pair<int, double>> &local_heaps) {
        double density = 0;
        for (unsigned int i = 0; i != progIO->numerical_parameters.num_neighbors_in_knn_search; i++) {
            // static member function, only have access to static data members, so use ds.tree
            density += ds.tree.particle_list[indices[i]].mass;
        }
        density /= progIO->numerical_parameters.four_PI_over_three * std::pow(radius, 3);
        return density;
    }

    /*! \fn template <class T> static double MedianSphericalDensityKernel(const DataSet<T, D> &ds, const double radius, uint32_t self_id, uint32_t *indices, std::vector<std::pair<int, double>> &local_heaps)
     *  \brief calculate density by dividing the total mass of certain neighbors by a spherical volume (radius = the distance to the "furthest" neighbor within two times the median)
     *  The purpose of removing outskirt particles is to avoid incorrect density calculation due to unreasonable large volume. From experiments, the ratio of the distance to the outermost particle and the median distance is almost between 1 and 2, with a small portion exceeding 2 and a few extremes. */
    template <class T>
    static double MedianSphericalDensityKernel(const DataSet<T, D> &ds, const double radius, uint32_t self_id, uint32_t *indices, std::vector<std::pair<int, double>> &local_heaps) {
        int i = progIO->numerical_parameters.num_neighbors_in_knn_search;
        double density = 0;
        double twice_median_radii_squared = 4 * local_heaps[i/2].second; // don't forget: "second" stores d^2

        for (i = i-1; i >= 0; i--) {
            if (local_heaps[i].second > twice_median_radii_squared) {
                continue;
            } else {
                break;
            }
        }
        twice_median_radii_squared = local_heaps[i].second;
        for ( ; i >= 0; i--) {
            density += ds.tree.particle_list[indices[i]].mass;
        }
        density /= progIO->numerical_parameters.four_PI_over_three * std::pow(twice_median_radii_squared, 1.5); // the distance info in the local_heaps is squared
        return density;
    }

    /*! \fn template <class T> double InverseDistanceWeightingDensityKernel(const DataSet<T, D> &ds, const double radius, uint32_t self_id, uint32_t *indices, std::vector<std::pair<int, double>> &local_heaps)
     *  \brief calculate density by dividing the total weighted mass of certain neighbors by a spherical volume (radius = the distance to the "furthest" neighbor within two times the median)
     *  The purpose of removing outskirt particles is to avoid incorrect density calculation due to unreasonable large volume.
     *  RL: From experiments, Inverse Distance Weighting (IDW) density kernel almost does not change the results of median spherical density kernel. */
    template <class T>
    static double InverseDistanceWeightingDensityKernel(const DataSet<T, D> &ds, const double radius, uint32_t self_id, uint32_t *indices, std::vector<std::pair<int, double>> &local_heaps) {
        // Applying the modified "Modified Shepard's method"
        // see https://en.wikipedia.org/wiki/Inverse_distance_weighting
        int i = progIO->numerical_parameters.num_neighbors_in_knn_search;
        double cut_radius = 0, mass = 0;
        double twice_median_radii_squared = 4 * local_heaps[i/2].second;

        for (i = i-1; i >= 0; i--) {
            if (local_heaps[i].second > twice_median_radii_squared) {
                continue;
            } else {
                break;
            }
        }
        cut_radius = local_heaps[i].second;

        for ( ; i >= 0; i--) {
            mass += ds.tree.particle_list[indices[i]].mass * (twice_median_radii_squared - local_heaps[i].second) / twice_median_radii_squared;
        }

        return mass /= progIO->numerical_parameters.four_PI_over_three * std::pow(cut_radius, 1.5);
    }

    /*! \fn template <class T> void RemoveSmallMassAndLowPeak(DataSet<T, D> &ds)
     *  \brief Remove those groups that only have small masses or have relatively low peak densities */
    template <class T>
    void RemoveSmallMassAndLowPeak(DataSet<T, D> &ds) {
        std::vector<uint32_t> peaks_to_be_deleted;
        
        for (auto &it : ds.planetesimal_list.planetesimals) {
            if (it.second.total_mass < ds.planetesimal_list.clump_mass_threshold || particle_list[it.first].new_density < ds.planetesimal_list.peak_density_threshold) {
                peaks_to_be_deleted.push_back(it.first);
                continue;
            }
            // RL: though we need the original density also above the threshold, we need to check around to make sure since the peak id particle might fool us
            uint32_t idx_limit = it.second.indices.size() < progIO->numerical_parameters.num_neighbors_in_knn_search ? it.second.indices.size() : progIO->numerical_parameters.num_neighbors_in_knn_search;
            double peak_ori_density = 0;
            for (uint32_t idx = 0; idx < idx_limit; idx++) {
                peak_ori_density = std::max(peak_ori_density, ds.tree.particle_list[it.second.indices[idx]].ath_density);
            }
            if (peak_ori_density < ds.planetesimal_list.peak_density_threshold) {
                peaks_to_be_deleted.push_back(it.first);
            }
        }
        for (auto &it : peaks_to_be_deleted) {
            ds.planetesimal_list.planetesimals.erase(it);
        }
        peaks_to_be_deleted.resize(0);
    }

    /*! \fn template <class T> void OutputNaivePeakList(DataSet<T, D> &ds, std::string filename, boost::dynamic_bitset<> &mask)
     *  \brief Output the list of naive peaks after hopping step */
    template <class T>
    void OutputNaivePeakList(DataSet<T, D> &ds, const std::string &filename, boost::dynamic_bitset<> &mask) {
        std::ofstream tmp_file(filename, std::ofstream::out);
        if (!tmp_file.is_open()) {
            std::cout << "Fail to open "+filename << std::endl;
        }
        tmp_file << "#" << std::setw(23) << "x" << std::setw(24) << "y" << std::setw(24) << "z" << std::setw(24) << "dis_max" << std::setw(24) << "Npar" << std::setw(24) << "R_1/10" << std::setw(24) << "R_HalfM" << std::setw(24) << "R_moreM" << std::endl;
        tmp_file << std::scientific;
        int idx = 0;
        for (auto &it : ds.planetesimal_list.planetesimals) {
            if (mask[idx]) {
                for (int i_dim = 0; i_dim != 3; i_dim++) {
                    tmp_file << std::setprecision(16) << std::setw(24) << it.second.center_of_mass[i_dim];
                }
                tmp_file << std::setprecision(16) << std::setw(24) << it.second.particles.back().second << std::setw(24) << it.second.particles.size() << std::setw(24) << it.second.one10th_radius << std::setw(24) << it.second.inner_one10th_radius << std::setw(24) << it.second.outer_one10th_radius << std::endl;
            }
            idx++;
        }
        tmp_file.close();
    }

    /*! \fn template <class T> void FindPlanetesimals(DataSet<T, D> &ds, int loop_count)
     *  \brief choose density kernel (overloaded) to find planetesimals in particle data
     *  For simplicity, we only use spatial information to identify planetesimals at first.
     *  Also, we force all density kernels to take the same format of arguments. */
    template <class T>
    void FindPlanetesimals(DataSet<T, D> &ds, int loop_count) {
        auto &paras = progIO->numerical_parameters;
        uint32_t horizontal_resolution = paras.box_resolution[0] * paras.box_resolution[1];
        if (ds.tree.half_width < paras.max_half_width + paras.max_ghost_zone_width - 1e-16) {
            // if particles are confined at a small region
            horizontal_resolution = static_cast<uint32_t>(std::pow(ds.tree.half_width / paras.cell_length[0] * ((ds.tree.half_width - paras.max_ghost_zone_width) / ds.tree.half_width), 2.0));
        }
        if (ds.particle_set.num_particles >= 4 * horizontal_resolution) {
            if (!paras.fixed_num_neighbors_to_hop) {
                if (ds.particle_set.num_particles >= 1.67e7) { // 4096^2=16777216
                    paras.num_neighbors_to_hop = 64;
                } else if (ds.particle_set.num_particles >= 2.68e8) { // 16384^2=268435456
                    paras.num_neighbors_to_hop = 128;
                }
            } else {
                if ( (ds.particle_set.num_particles >= 1.67e7
                      && paras.num_neighbors_to_hop < 64)
                     || (ds.particle_set.num_particles >= 2.68e8
                         && paras.num_neighbors_to_hop < 128) ) {
                    progIO->log_info << "The number of particles retrieved from data is quite a lot, " << ds.particle_set.num_particles << ", you may want to set a larger num_neighbors_to_hop by specifying \"hop\" in the parameter input file (current value is " << paras.num_neighbors_to_hop << ")." << std::endl;
                    progIO->Output(std::cout, progIO->log_info, __normal_output, __master_only);
                }
            }
            enough_particle_resolution_flag = 1;
            FindPlanetesimals(ds, BHtree<dim>::MedianSphericalDensityKernel<float>, loop_count);
        } else if (ds.particle_set.num_particles >= horizontal_resolution) {
            enough_particle_resolution_flag = 0;
            progIO->log_info << "The number of particles retrieved from data files is merely " << ds.particle_set.num_particles / horizontal_resolution << " times the horizontal grid resolution. Consider using more particles in simulations or output more particle data for better results." << std::endl;
            progIO->Output(std::cout, progIO->log_info, __more_output, __master_only);
            FindPlanetesimals(ds, BHtree<dim>::QseudoQuadraticSplinesKernel<float>, loop_count);
        } else {
            enough_particle_resolution_flag = 0;
            progIO->log_info << "The number of particles retrieved from data files is only " << ds.particle_set.num_particles << ", a fraction of the horizontal grid resolution " << horizontal_resolution << ". Using more particles in simulations or output more particle data is strongly recommended for more reliable results." << std::endl;
            progIO->Output(std::cout, progIO->log_info, __normal_output, __all_processors);
            FindPlanetesimals(ds, BHtree<dim>::QseudoQuadraticSplinesKernel<float>, loop_count);
        }
    }

    /*! \fn template <class T, class F> void FindPlanetesimals(DataSet<T, D> &ds, F f, int loop_count)
     *  \brief find planetesimals in particle data
     *  \tparam F the function that used to calculate new particle density based on a given kernel
     *  For simplicity, we only use spatial information to identify planetesimals at first. */
    template <class T, class F>
    void FindPlanetesimals(DataSet<T, D> &ds, F DensityKernel, int loop_count) {

        ds.planetesimal_list.planetesimals.clear();
        ds.planetesimal_list.num_planetesimals = 0;
        ds.planetesimal_list.peaks_and_masses.resize(0);
        assert(ds.planetesimal_list.planetesimals.size() == 0);
        assert(ds.planetesimal_list.planetesimals.size() == 0);
        std::vector<uint32_t> peaks_to_be_deleted;

        /* Toomre Q ~ Omega^2 / PI / grav_constant / rho < 1
         * ==> rho_crit > Omega^2 / PI / grav_constant
         * so we are looking for regions/particles with at least rho_crit,
         * for now we use 2 * rho_crit as a threshold
         */
        ds.planetesimal_list.density_threshold = std::pow(progIO->numerical_parameters.Omega, 2.) / progIO->numerical_parameters.PI / progIO->numerical_parameters.grav_constant * 2.; // = 160 for the fiducial run
        // if the # of particles is low (m_par is large), we'd better use 9 * rho_0 / tilde_G
        if (progIO->numerical_parameters.mass_total_code_units / ds.particle_set.num_particles / progIO->numerical_parameters.cell_volume >  ds.planetesimal_list.density_threshold / 4) {
            ds.planetesimal_list.density_threshold *= 9. / 8.;
        }
        // if the hydro-resolution is low (coarse FFT-solver => particles not highly concentrated)
        // we need to avoid diffuse particle groups (diffuse like randomly distributed)
        // also we need the full Hill radius in merging b/c particles are more puffed in low resolution
        auto hydro_res_per_H = static_cast<unsigned int>((progIO->numerical_parameters.box_resolution / progIO->numerical_parameters.box_length).MinElement());
        if (hydro_res_per_H <= 1024) {
            ds.planetesimal_list.clump_diffuse_threshold = 0.35;
            ds.planetesimal_list.Hill_fraction_for_merge = 0.75;
        } else if (hydro_res_per_H <= 1536) {
            ds.planetesimal_list.clump_diffuse_threshold = 0.4;
            ds.planetesimal_list.Hill_fraction_for_merge = 0.5;
        } else if (hydro_res_per_H <= 2048) {
            ds.planetesimal_list.clump_diffuse_threshold = 0.5;
            ds.planetesimal_list.Hill_fraction_for_merge = 0.35;
        }
        // and we are looking for clumps with certain mass_crit & peak_rho_crit
        ds.planetesimal_list.clump_mass_threshold = progIO->numerical_parameters.min_trusted_mass_code_unit;
        ds.planetesimal_list.peak_density_threshold = 3 * ds.planetesimal_list.density_threshold; // from HOP's experience (ref: https://www.cfa.harvard.edu/~deisenst/hop/hop_doc.html)

        /////////////////////////////////////////////////////
        // 1, Calculate the particle densities
        double radius_Kth_NN = 0.;
        uint32_t *indices;
        //std::vector<double> other_density (num_particles, 0);
#ifdef OpenMP_ON
        omp_set_num_threads(progIO->numerical_parameters.num_avail_threads);
#pragma omp parallel private(radius_Kth_NN, indices)
        {
            indices = new uint32_t[progIO->numerical_parameters.num_neighbors_in_knn_search];
            std::vector<std::pair<int, double>> local_heaps;
#pragma omp for
            for (uint32_t i = 0; i < num_particles; i++) {
                KNN_Search_OpenMP(particle_list[i].pos, progIO->numerical_parameters.num_neighbors_in_knn_search, radius_Kth_NN, indices, local_heaps);
                particle_list[i].new_density = DensityKernel(ds, radius_Kth_NN, i, indices, local_heaps);
            }
            delete[] indices;
            indices = nullptr;
        }
#else // OpenMP_ON
        indices = new uint32_t[progIO->numerical_parameters.num_neighbors_in_knn_search];
        for (uint32_t i = 0; i != num_particles; i++) {
            KNN_Search(particle_list[i].pos, progIO->numerical_parameters.num_neighbors_in_knn_search, radius_Kth_NN, indices, true);
            particle_list[i].new_density = DensityKernel(ds, radius_Kth_NN, i, indices, heaps);
        }
        delete [] indices;
        indices = nullptr;
#endif // OpenMP_ON
        progIO->log_info << "Density calculation done, max_dpar = "
                         << std::max_element(particle_list, particle_list + num_particles,
                                             [](const InternalParticle &p1, const InternalParticle &p2) {
                                                 return p1.new_density < p2.new_density;
                                             }
                         )->new_density << ", dpar_threshold=" << ds.planetesimal_list.density_threshold << "; ";

        // RL: debug use -- print out the new_density for comparison (with the density output by Athena)
        /*
        std::ofstream dpar_file;
        dpar_file.open("dpar.txt");
        for (uint32_t i = 0; i != num_particles; i++) {
            dpar_file.unsetf(std::ios_base::floatfield);
            dpar_file << std::setw(9) << i;
            dpar_file << std::scientific << std::setw(16) << particle_list[i].new_density << std::setw(16) << particle_list[i].ath_density << std::endl;
        }
        dpar_file.close();
        //*/

        /////////////////////////////////////////////////////
        // 2, Identify the neighbor with the highest density (if same, then choose the one with a smaller original ID)
        auto *densest_neighbor_id_list = new uint32_t[num_particles]();
#ifdef OpenMP_ON
        omp_set_num_threads(progIO->numerical_parameters.num_avail_threads);
#pragma omp parallel private(radius_Kth_NN, indices)
        {
            indices = new uint32_t[progIO->numerical_parameters.num_neighbors_to_hop];
            std::vector<std::pair<int, double>> local_heaps;
#pragma omp for schedule(auto)
            for (uint32_t i = 0; i < num_particles; i++) {
                densest_neighbor_id_list[i] = i;
                if (particle_list[i].new_density < ds.planetesimal_list.density_threshold) {
                    continue;
                }
                KNN_Search_OpenMP(particle_list[i].pos, progIO->numerical_parameters.num_neighbors_to_hop, radius_Kth_NN, indices, local_heaps);
                double tmp_density = particle_list[i].new_density;
                for (unsigned int j = 0; j != progIO->numerical_parameters.num_neighbors_to_hop; j++) {
                    double density_j = particle_list[indices[j]].new_density;
                    if (tmp_density < density_j) {
                        densest_neighbor_id_list[i] = indices[j];
                        tmp_density = density_j;
                    } else if (tmp_density == density_j && densest_neighbor_id_list[i] > indices[j]) {
                        densest_neighbor_id_list[i] = indices[j];
                    }
                }
            }
            delete [] indices;
            indices = nullptr;
        }
#else // OpenMP_ON
        indices = new uint32_t[progIO->numerical_parameters.num_neighbors_to_hop];
        for (uint32_t i = 0; i < num_particles; i++) {
            densest_neighbor_id_list[i] = i;
            if (particle_list[i].new_density < ds.planetesimal_list.density_threshold) {
                continue;
            }
            KNN_Search(particle_list[i].pos, progIO->numerical_parameters.num_neighbors_to_hop, radius_Kth_NN, indices);
            double tmp_density = particle_list[i].new_density;
            for (unsigned int j = 0; j != progIO->numerical_parameters.num_neighbors_to_hop; j++) {
                double density_j = particle_list[indices[j]].new_density;
                if (tmp_density < density_j) {
                    densest_neighbor_id_list[i] = indices[j];
                    tmp_density = density_j;
                } else if (tmp_density == density_j && densest_neighbor_id_list[i] > indices[j]) {
                    densest_neighbor_id_list[i] = indices[j];
                }
            }
        }
        delete [] indices;
        indices = nullptr;
#endif // OpenMP_ON

        progIO->log_info << "Densest neighbor found, particle_list[0]'s densest_neighbor_id = " << densest_neighbor_id_list[0] << "; ";

        /////////////////////////////////////////////////////
        // 3, hop the densest neighbor all the way to the peak
        boost::dynamic_bitset<> mask(num_particles);
        mask.set(); // set all the bits to 1
        auto  *tmp_peak_id_list = new uint32_t[num_particles]();
        for (uint32_t i = 0; i != num_particles; i++) {
            if (!mask[i]) {
                continue;
            }
            if (particle_list[i].new_density < ds.planetesimal_list.density_threshold) {
                mask.flip(i);
                tmp_peak_id_list[i] = num_particles; // signal to be ignored
                continue;
            }
            std::vector<uint32_t> tmp_indices;
            tmp_peak_id_list[i] = i;
            uint32_t chain = i;
            tmp_indices.push_back(chain);
            mask.flip(chain);
            chain = densest_neighbor_id_list[i];
            // hop to denser region, but stop if we hit a particle that is already considered
            while (tmp_indices.back() != chain) {
                if (mask[chain]) {
                    tmp_indices.push_back(chain);
                    mask.flip(chain);
                    chain = densest_neighbor_id_list[chain];
                } else {
                    break;
                }
            }
            if (tmp_indices.back() == chain) { // if hop up to the peak
                for (auto it : tmp_indices) {
                    tmp_peak_id_list[it] = chain;
                }
            } else { // if not, which means breaking at some point
                for (auto it : tmp_indices) {
                    tmp_peak_id_list[it] = tmp_peak_id_list[chain];
                }
            }
            // check if we need to create a planetesimal object or just append particle list
            auto it = ds.planetesimal_list.planetesimals.emplace(tmp_peak_id_list[i], Planetesimal<D>());
            if (it.second) { // if insertion happens
                it.first->second.peak_index = tmp_peak_id_list[i];
            }
            it.first->second.indices.insert(it.first->second.indices.end(), tmp_indices.begin(), tmp_indices.end());
        }
        assert(mask.none()); // check if all the bits are 0
        mask.clear();
        delete [] tmp_peak_id_list;
        tmp_peak_id_list = nullptr;
        delete [] densest_neighbor_id_list;
        densest_neighbor_id_list = nullptr;
        for (auto &it : ds.planetesimal_list.planetesimals) {
            it.second.SortParticles(particle_list);
            if (ds.planetesimal_list.clump_diffuse_threshold < 0.525) {
                if (it.second.IsPositionDispersion2Large(ds.planetesimal_list.clump_diffuse_threshold)) {
                    peaks_to_be_deleted.push_back(it.first);
                }
            }
        }
        if (ds.planetesimal_list.clump_diffuse_threshold < 0.525) {
            for (auto it : peaks_to_be_deleted) {
                ds.planetesimal_list.planetesimals.erase(it);
            }
            peaks_to_be_deleted.resize(0);
        }
        ds.planetesimal_list.num_planetesimals = static_cast<uint32_t>(ds.planetesimal_list.planetesimals.size());

        /* With certain parameters, there will be A LOT clumps all over the computational domain. Hopping will also identify very loose associations, where particle densities are marginally higher than the threshold. Such associations are in fact sparse collections of a few particles. For efficiency, we remove them immediately.
         * RL update: For the complete version of particle data, the number of particles is so large that the hopping step need a large value of num_neighbors_to_hop to identify meaningful loose associations (like sub-sample). */
        if (enough_particle_resolution_flag == 0) {
            RemoveSmallMassAndLowPeak(ds);
        }
        ds.planetesimal_list.num_planetesimals = static_cast<uint32_t>(ds.planetesimal_list.planetesimals.size());
        progIO->log_info << "Hopping to peaks done, naive peak amount = " << ds.planetesimal_list.planetesimals.size() << "; ";

        // RL: debug use -- print out particle groups' properties after hopping and before merging
        /*
        std::ofstream tmp_file("naive_peak_list.txt", std::ofstream::out);
        if (!tmp_file.is_open()) {
            std::cout << "Fail to open naive_peak_list.txt" << std::endl;
        }
        tmp_file << "#" << std::setw(23) << "x" << std::setw(24) << "y" << std::setw(24) << "z" << std::setw(24) << "dis_max" << std::setw(24) << "Npar" << std::setw(24) << "R_1/10" << std::setw(24) << "R_HalfM" << std::setw(24) << "R_moreM" << std::endl;
        tmp_file << std::scientific;
        for (auto &it : ds.planetesimal_list.planetesimals) {
            for (int i_dim = 0; i_dim != 3; i_dim++) {
                tmp_file << std::setprecision(16) << std::setw(24) << it.second.center_of_mass[i_dim];
            }
            tmp_file << std::setprecision(16) << std::setw(24) << it.second.particles.back().second << std::setw(24) << it.second.particles.size() << std::setw(24) << it.second.one10th_radius << std::setw(24) << it.second.inner_one10th_radius << std::setw(24) << it.second.outer_one10th_radius << std::endl;
        }
        tmp_file.close();
        //*/

        /* RL: Originally, the next step should be find the boundary particles and catalog the highest density boundary found between each pair of groups. After that, the code determines whether each naive clumps are viable or unviable depending on a peak_threshold and then merge two viable groups if they touch with each other and their boundary density > merge_threshold. An unviable clump is merged to the viable group with which it shares the highest boundary density. At last, all particles whose density are less than a outer_threshold are excluded from groups.
         * But here I found the naive peaks found so far already obviously stand out of the background in the form of individual clump-clusters. Thus, from here, I'm going to merge those naive clumps if: (a) one contains another smaller one; (b) one touches another one and they are gravitationally bound with the other. After that, we'll go through each outer particle in those naive planetesimals to check if they are gravitationally bound.
         */

        /////////////////////////////////////////////////////
        // 4, merge bound clumps
        int merge_happened_flag = 1;
        double max_radius = 0.;
        uint32_t *nearby_mask;
        unsigned int merging_count = 0, delete_count = 0, predator_count = 0;
        using pl_iterator = decltype(ds.planetesimal_list.planetesimals.begin());
        std::vector<pl_iterator> nearby_pi;
        std::vector<std::pair<pl_iterator, pl_iterator>> merging_pairs;
        std::vector<std::pair<pl_iterator, std::vector<pl_iterator>>> combined_merging_pairs;
        nearby_pi.resize(ds.planetesimal_list.planetesimals.size());
        merging_pairs.resize(ds.planetesimal_list.num_planetesimals);
        combined_merging_pairs.resize(ds.planetesimal_list.num_planetesimals);
        peaks_to_be_deleted.resize(ds.planetesimal_list.num_planetesimals);

        while (merge_happened_flag) {
            // Build Clump Tree every time is necessary
            ds.planetesimal_list.BuildClumpTree(root_center, half_width, max_radius);
            merge_happened_flag = 0;
            merging_count = 0;
            predator_count = 0;
            nearby_mask = new uint32_t[ds.planetesimal_list.num_planetesimals]();
#ifdef OpenMP_ON
            omp_set_num_threads(progIO->numerical_parameters.num_avail_threads);
#pragma omp parallel
#endif
            {
                auto tmp_p1 = ds.planetesimal_list.planetesimals.begin();
                auto tmp_p2 = tmp_p1;
                double r_p1 = 0, r_p2 = 0, center_dist = 0;
                uint32_t nearby_count = 0;
                auto *nearby_indices = new uint32_t[ds.planetesimal_list.planetesimals.size()];

                for (tmp_p1 = ds.planetesimal_list.planetesimals.begin(); tmp_p1 != ds.planetesimal_list.planetesimals.end(); tmp_p1++) {
                    r_p1 = tmp_p1->second.one10th_radius;
                    if (tmp_p1->second.mask) {
                        // BallSearch does not sort the results, don't assume 1st one is itself
                        ds.planetesimal_list.clump_tree.BallSearch(tmp_p1->second.center_of_mass, tmp_p1->second.one10th_radius+max_radius, nearby_indices, nearby_count);
#ifdef OpenMP_ON
#pragma omp for
#endif
                        for (size_t idx = 0; idx < nearby_count; idx++) {
                            uint32_t tmp_p2_id = ds.planetesimal_list.clump_tree.particle_list[nearby_indices[idx]].ath_density; // this density stores peak_id
                            if (tmp_p2_id == tmp_p1->first) {
                                continue;
                            }
                            tmp_p2 = ds.planetesimal_list.planetesimals.find(tmp_p2_id);
                            if (tmp_p2 != ds.planetesimal_list.planetesimals.end()) {

                                if (tmp_p2->second.mask) {
                                    r_p2 = tmp_p2->second.one10th_radius;
                                    center_dist = (tmp_p1->second.center_of_mass - tmp_p2->second.center_of_mass).Norm();
                                    if (center_dist < std::max(r_p1, r_p2)) {
                                        // one contains the center of another, just merge
                                        nearby_mask[idx] = 1;
                                        nearby_pi[idx] = tmp_p2;
                                        tmp_p2->second.mask = false;
                                    } else if (center_dist < (r_p1 + r_p2) && ds.planetesimal_list.IsGravitationallyBound(tmp_p1->second, tmp_p2->second)) {
                                        // they intersect and bound
                                        if (!ds.planetesimal_list.IsSaddlePointDeepEnough(ds, DensityKernel, tmp_p1->second, tmp_p2->second)) {
                                            nearby_mask[idx] = 1;
                                            nearby_pi[idx] = tmp_p2;
                                            tmp_p2->second.mask = false;
                                        }
                                    }
                                }
                            } else {
                                progIO->error_message << "Error: Cannot find a clump with peak_id=" << tmp_p2_id << ". This should not happen. Please report a bug. Proceed for now." << std::endl;
                                progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                            }
                        }
                    }
#ifdef OpenMP_ON
#pragma omp barrier
#pragma omp single
#endif
                    {
                        for (size_t idx = 0; idx < nearby_count; idx++) {
                            if (nearby_mask[idx]) {
                                merging_pairs[merging_count].first = tmp_p1;
                                merging_pairs[merging_count].second = nearby_pi[idx];
                                merging_count++;
                                peaks_to_be_deleted[delete_count] = nearby_pi[idx]->first;
                                delete_count++;
                                nearby_mask[idx] = 0;
                            }
                        }
                        if (merging_count > 0) {
                            merge_happened_flag = 1;
                        }
                    }
                }
                delete [] nearby_indices;
                nearby_indices = nullptr;
#ifdef OpenMP_ON
#pragma omp single
#endif
                {
                    std::vector<decltype(ds.planetesimal_list.planetesimals.begin())> preys;
                    auto idx_limit = merging_count-1;
                    for (unsigned int idx = 0; idx != merging_count; ) {
                        preys.resize(0);
                        preys.push_back(merging_pairs[idx].second);
                        while (idx < idx_limit && merging_pairs[idx].first->first == merging_pairs[idx+1].first->first) {
                            preys.push_back(merging_pairs[idx+1].second);
                            idx++;
                        }
                        combined_merging_pairs[predator_count] = std::make_pair(merging_pairs[idx].first, preys);
                        predator_count++;
                        idx++;
                    }
                }
#ifdef OpenMP_ON
#pragma omp for
#endif
                for (uint32_t i = 0; i < predator_count; i++) {
                    for (auto &it : combined_merging_pairs[i].second) {
                        combined_merging_pairs[i].first->second.MergeAnotherPlanetesimal(it->second, particle_list);
                    }
                }
            }
            delete [] nearby_mask;
            nearby_mask = nullptr;
            for (auto &it : ds.planetesimal_list.planetesimals) {
                if (it.second.mask) {
                    if (it.second.IsPositionDispersion2Large()) {
                        it.second.mask = false;
                        peaks_to_be_deleted[delete_count] = it.first;
                        delete_count++;
                    }
                }
            }
            ds.planetesimal_list.num_planetesimals = ds.planetesimal_list.planetesimals.size() - delete_count;
            /* RL: debug use
            std::cout << "peaks to be deleted now equals " << delete_count << std::endl;
            std::cout << "merging performed " << merging_count << " times" << std::endl;
            //*/
        }

        for (auto it : peaks_to_be_deleted) {
            ds.planetesimal_list.planetesimals.erase(it);
        }
        merging_pairs.resize(0);
        combined_merging_pairs.resize(0);
        peaks_to_be_deleted.resize(0);
        nearby_pi.resize(0);

        /*
         * RL: now we merge primitive clumps according to a fraction of their Hill radii
         */
        merge_happened_flag = 1;
        merging_count = 0; delete_count = 0; predator_count = 0; max_radius = 0;
        for (auto &it : ds.planetesimal_list.planetesimals) {
            it.second.CalculateHillRadius();
        }
        nearby_pi.resize(ds.planetesimal_list.planetesimals.size());
        merging_pairs.resize(ds.planetesimal_list.num_planetesimals);
        combined_merging_pairs.resize(ds.planetesimal_list.num_planetesimals);
        peaks_to_be_deleted.resize(ds.planetesimal_list.num_planetesimals);

        while (merge_happened_flag) {
            // Build Clump Tree every time is necessary
            ds.planetesimal_list.BuildClumpTree(root_center, half_width, max_radius, true);
            merge_happened_flag = 0;
            merging_count = 0;
            predator_count = 0;
            nearby_mask = new uint32_t[ds.planetesimal_list.num_planetesimals]();
#ifdef OpenMP_ON
            omp_set_num_threads(progIO->numerical_parameters.num_avail_threads);
#pragma omp parallel
#endif
            {
                auto tmp_p1 = ds.planetesimal_list.planetesimals.begin();
                auto tmp_p2 = tmp_p1;
                double r_p1 = 0, r_p2 = 0, center_dist = 0;
                uint32_t nearby_count = 0;
                auto *nearby_indices = new uint32_t[ds.planetesimal_list.planetesimals.size()];

                for (tmp_p1 = ds.planetesimal_list.planetesimals.begin(); tmp_p1 != ds.planetesimal_list.planetesimals.end(); tmp_p1++) {
                    r_p1 = tmp_p1->second.Hill_radius * ds.planetesimal_list.Hill_fraction_for_merge;
                    if (tmp_p1->second.mask) {
                        // BallSearch does not sort the results, don't assume 1st one is itself
                        ds.planetesimal_list.clump_tree.BallSearch(tmp_p1->second.center_of_mass, tmp_p1->second.Hill_radius * ds.planetesimal_list.Hill_fraction_for_merge + max_radius, nearby_indices, nearby_count);
#ifdef OpenMP_ON
#pragma omp for
#endif
                        for (size_t idx = 0; idx < nearby_count; idx++) {
                            uint32_t tmp_p2_id = ds.planetesimal_list.clump_tree.particle_list[nearby_indices[idx]].ath_density; // this density stores peak_id
                            if (tmp_p2_id == tmp_p1->first) {
                                continue;
                            }
                            tmp_p2 = ds.planetesimal_list.planetesimals.find(tmp_p2_id);
                            if (tmp_p2 != ds.planetesimal_list.planetesimals.end()) {

                                if (tmp_p2->second.mask) {
                                    r_p2 = tmp_p2->second.Hill_radius * ds.planetesimal_list.Hill_fraction_for_merge;
                                    center_dist = (tmp_p1->second.center_of_mass - tmp_p2->second.center_of_mass).Norm();
                                    if (center_dist < std::max(r_p1, r_p2)) {
                                        // one contains the center of another, just merge
                                        nearby_mask[idx] = 1;
                                        nearby_pi[idx] = tmp_p2;
                                        tmp_p2->second.mask = false;
                                    } else if (center_dist < (r_p1 + r_p2) && ds.planetesimal_list.IsGravitationallyBound(tmp_p1->second, tmp_p2->second)) {
                                        // they intersect and bound
                                        if (!ds.planetesimal_list.IsHillSaddlePointDeepEnough(ds, DensityKernel, tmp_p1->second, tmp_p2->second)) {
                                            nearby_mask[idx] = 1;
                                            nearby_pi[idx] = tmp_p2;
                                            tmp_p2->second.mask = false;
                                        }
                                    }
                                }
                            } else {
                                progIO->error_message << "Error: Cannot find a clump with peak_id=" << tmp_p2_id << ". This should not happen. Please report a bug. Proceed for now." << std::endl;
                                progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
                            }
                        }
                    }
#ifdef OpenMP_ON
#pragma omp barrier
#pragma omp single
#endif
                    {
                        for (size_t idx = 0; idx < nearby_count; idx++) {
                            if (nearby_mask[idx]) {
                                merging_pairs[merging_count].first = tmp_p1;
                                merging_pairs[merging_count].second = nearby_pi[idx];
                                merging_count++;
                                peaks_to_be_deleted[delete_count] = nearby_pi[idx]->first;
                                delete_count++;
                                nearby_mask[idx] = 0;
                            }
                        }
                        if (merging_count > 0) {
                            merge_happened_flag = 1;
                        }
                    }
                }
                delete [] nearby_indices;
                nearby_indices = nullptr;
#ifdef OpenMP_ON
#pragma omp single
#endif
                {
                    std::vector<decltype(ds.planetesimal_list.planetesimals.begin())> preys;
                    auto idx_limit = merging_count-1;
                    for (unsigned int idx = 0; idx != merging_count; ) {
                        preys.resize(0);
                        preys.push_back(merging_pairs[idx].second);
                        while (idx < idx_limit && merging_pairs[idx].first->first == merging_pairs[idx+1].first->first) {
                            preys.push_back(merging_pairs[idx+1].second);
                            idx++;
                        }
                        combined_merging_pairs[predator_count] = std::make_pair(merging_pairs[idx].first, preys);
                        predator_count++;
                        idx++;
                    }
                }
#ifdef OpenMP_ON
#pragma omp for
#endif
                for (uint32_t i = 0; i < predator_count; i++) {
                    for (auto &it : combined_merging_pairs[i].second) {
                        combined_merging_pairs[i].first->second.MergeAnotherPlanetesimal(it->second, particle_list);
                        combined_merging_pairs[i].first->second.CalculateHillRadius();
                    }
                }
            }
            delete [] nearby_mask;
            nearby_mask = nullptr;
            for (auto &it : ds.planetesimal_list.planetesimals) {
                if (it.second.mask) {
                    if (it.second.IsPositionDispersion2Large()) {
                        it.second.mask = false;
                        peaks_to_be_deleted[delete_count] = it.first;
                        delete_count++;
                    }
                }
            }
            ds.planetesimal_list.num_planetesimals = ds.planetesimal_list.planetesimals.size() - delete_count;
            /* RL: debug use
            std::cout << "peaks to be deleted now equals " << delete_count << std::endl;
            std::cout << "merging performed " << merging_count << " times" << std::endl;
            //*/
        }

        for (auto it : peaks_to_be_deleted) {
            ds.planetesimal_list.planetesimals.erase(it);
        }
        merging_pairs.resize(0);
        combined_merging_pairs.resize(0);
        peaks_to_be_deleted.resize(0);
        nearby_pi.resize(0);

        // put real peak indices into keys
        // \todo: if this happens, I also need to change the potential_subpeak_indices
        std::vector<uint32_t> real_peak_indices;
        for (auto &it : ds.planetesimal_list.planetesimals) {
            uint32_t tmp_peak_id = it.first;
            for (auto it_other_peak : it.second.potential_subpeak_indices) {
                if (particle_list[tmp_peak_id].new_density < particle_list[it_other_peak].new_density) {
                    tmp_peak_id = it_other_peak;
                }
            }
            if (tmp_peak_id != it.first) {
                peaks_to_be_deleted.push_back(it.first);
                real_peak_indices.push_back(tmp_peak_id);
            }
        }
        for (unsigned int i = 0; i != peaks_to_be_deleted.size(); i++) {
            auto it = ds.planetesimal_list.planetesimals.emplace(real_peak_indices[i], Planetesimal<D>());
            if (it.second) { // if insertion happens successfully, which it should
                auto original_item = ds.planetesimal_list.planetesimals.find(peaks_to_be_deleted[i]);
                std::swap(original_item->second, it.first->second);
                it.first->second.potential_subpeak_indices.push_back(original_item->first);
                ds.planetesimal_list.planetesimals.erase(original_item);
            } else {
                progIO->error_message << "Real peak index replacing failed. " << std::endl;
                progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
            }
        }
        peaks_to_be_deleted.resize(0);
        real_peak_indices.resize(0); //*/

        progIO->log_info << "Merging done, now " << ds.planetesimal_list.planetesimals.size() << " left; ";

        /////////////////////////////////////////////////////
        // 5, remove unbound particles and search potential bound particles within Hill Radius
        // \todo: we may need to collect these unbound particles --> we may find some sub-clumps
        for (auto &it : ds.planetesimal_list.planetesimals) {
            it.second.RemoveUnboundParticles(particle_list);
            it.second.CalculateHillRadius();
        }
        progIO->log_info << "Unbinding particles done. ";

        // mark particles that are already in clumps
        for (auto &it : ds.planetesimal_list.planetesimals) {
            for (auto index : it.second.indices) {
                particle_list[index].in_clump_flag = true;
            }
        }
        for (auto &it : ds.planetesimal_list.planetesimals) {
            size_t tmp_num_particles;
            do {
                tmp_num_particles = it.second.indices.size();
                it.second.SearchBoundParticlesWithinHillRadius(ds.tree, ds.planetesimal_list.density_threshold);
                it.second.CalculateHillRadius();
            } while (it.second.indices.size() > tmp_num_particles);
        }
        progIO->log_info << "Rebinding particles done. ";

        /////////////////////////////////////////////////////
        // 6, remove those groups that only have small masses or have relatively low peak densities
        RemoveSmallMassAndLowPeak(ds);
        for (auto &it : ds.planetesimal_list.planetesimals) { // also remove losse streams
            if (it.second.mask) {
                if (it.second.IsPositionDispersion2Large(std::min(ds.planetesimal_list.clump_diffuse_threshold+0.1, 0.475))) {
                    peaks_to_be_deleted.push_back(it.first);
                }
            }
        }
        for (auto &it : peaks_to_be_deleted) {
            ds.planetesimal_list.planetesimals.erase(it);
        }
        peaks_to_be_deleted.resize(0);
        progIO->log_info << "Remove low&small peaks, now " << ds.planetesimal_list.planetesimals.size() << " left; ";

        /////////////////////////////////////////////////////
        // 7, even though we remove unbound particles, there are still plenty of particles that are marginally attached to our planetesimals (outside the naive Hill radius, which is in fact larger than the really radius since we include more mass). Let's try to remove them (thinking about removing puffed-envelopes).
        int ripping_outerior_count = 0;
        for (auto &it : ds.planetesimal_list.planetesimals) {
            while (it.second.Hill_radius < it.second.particles.back().second) {
                auto low = std::lower_bound(it.second.particles.begin(), it.second.particles.end(), std::pair<uint32_t, double>(0, it.second.Hill_radius), less_second<uint32_t, double>());
                ripping_outerior_count++;
                it.second.particles.resize(low - it.second.particles.begin());
                it.second.indices.resize(it.second.particles.size());
                if (it.second.indices.size() == 0) {
                    peaks_to_be_deleted.push_back(it.first);
                    break;
                }
                it.second.CalculateKinematicProperties(particle_list);
                it.second.RemoveUnboundParticles(particle_list);
                it.second.CalculateHillRadius();
            }
        }
        for (auto it : peaks_to_be_deleted) {
            ds.planetesimal_list.planetesimals.erase(it);
        }
        peaks_to_be_deleted.resize(0);
        progIO->log_info << "Ripping Envelope (" << ripping_outerior_count << "); ";

        /////////////////////////////////////////////////////
        // 8, again remove those groups that only have small masses relatively low peak densities
        RemoveSmallMassAndLowPeak(ds);
        progIO->log_info << "Remove low&small peaks again, now " << ds.planetesimal_list.planetesimals.size() << " left; ";

        /////////////////////////////////////////////////////
        // 9, remove those planetesimals that locate at ghost zones
        for (auto &it : ds.planetesimal_list.planetesimals) {
            if (!it.second.center_of_mass.InRange(progIO->numerical_parameters.box_min, progIO->numerical_parameters.box_max)) {
                peaks_to_be_deleted.push_back(it.first);
            }
        }
        for (auto &it : peaks_to_be_deleted) {
            ds.planetesimal_list.planetesimals.erase(it);
        }
        peaks_to_be_deleted.resize(0);
        progIO->log_info << "Erase ghost planetesimals, now " << ds.planetesimal_list.planetesimals.size() << " left; ";
        ds.planetesimal_list.num_planetesimals = static_cast<uint32_t >(ds.planetesimal_list.planetesimals.size());
        progIO->out_content << "Finish clump finding for t = " << ds.particle_set.time << ", ";

        /////////////////////////////////////////////////////
        // 10, perform more scientific calculations
        for (auto &it : ds.planetesimal_list.planetesimals) {
            it.second.CalculateAngularMomentum(ds.tree);
        }

        // RL: output the cumulative mass/Jz function
        /*
        std::ofstream cumulative_file("cumulative_Jz.txt", std::ofstream::out);
        if (!cumulative_file.is_open()) {
            std::cout << "Fail to open cumulative_Jz.txt" << std::endl;
        }
        for (auto &it : ds.planetesimal_list.planetesimals) {
            it.second.CalculateCumulativeAngularMomentum(ds.tree, cumulative_file);
        }
        cumulative_file.close();
        //*/

        /////////////////////////////////////////////////////
        // Before ending, output some info
        if (ds.planetesimal_list.num_planetesimals > 0) {
            // construct "peaks_and_masses" in order to sort clumps by mass
            for (auto &it : ds.planetesimal_list.planetesimals) {
                ds.planetesimal_list.peaks_and_masses.push_back(std::pair<uint32_t, double>(it.first, it.second.total_mass));
            }
            std::sort(ds.planetesimal_list.peaks_and_masses.begin(), ds.planetesimal_list.peaks_and_masses.end(), [](const std::pair<uint32_t, double> &a, const std::pair<uint32_t, double> &b) {
                if (a.second == b.second) {
                    return a.first < b.first; // for clumps with the same # of particles, ranking by peak_id
                }
                return a.second < b.second;
            });

            if (progIO->numerical_parameters.num_peaks > 0 && ds.planetesimal_list.num_planetesimals > progIO->numerical_parameters.num_peaks) {
                size_t num_to_be_deleted = ds.planetesimal_list.num_planetesimals - progIO->numerical_parameters.num_peaks;
                for (auto it = ds.planetesimal_list.peaks_and_masses.begin(); it != ds.planetesimal_list.peaks_and_masses.begin() + num_to_be_deleted; it++) {
                    peaks_to_be_deleted.push_back(it->first);
                }
                std::vector<typename decltype(ds.planetesimal_list.peaks_and_masses)::value_type>(
                        ds.planetesimal_list.peaks_and_masses.begin() + num_to_be_deleted,
                        ds.planetesimal_list.peaks_and_masses.end()).swap(ds.planetesimal_list.peaks_and_masses);
                for (auto it : peaks_to_be_deleted) {
                    ds.planetesimal_list.planetesimals.erase(it);
                }
                peaks_to_be_deleted.resize(0);
                ds.planetesimal_list.num_planetesimals = static_cast<uint32_t >(ds.planetesimal_list.planetesimals.size());

                progIO->log_info << "Remove low&small peaks once more due to the input max num_peaks limit, now " << ds.planetesimal_list.planetesimals.size() << " left; ";
            }

            double Mp_tot = std::accumulate(ds.planetesimal_list.peaks_and_masses.begin(), ds.planetesimal_list.peaks_and_masses.end(), 0., [](const double &a, const std::pair<uint32_t, double> &b) {
                return a + b.second;
            });
            progIO->out_content << "found " << ds.planetesimal_list.num_planetesimals << " clumps; " << " Mp_max = " << ds.planetesimal_list.peaks_and_masses.back().second << ", Mp_tot = " << Mp_tot << "(" << std::fixed << Mp_tot/progIO->numerical_parameters.mass_total_code_units*100 << "%) in code units.";
        } else {
            progIO->out_content << "found zero clumps";
        }
        if (loop_count >= 0) { // RL: input negative loop_count for test purpose
            ds.planetesimal_list.OutputPlanetesimalsInfo(loop_count, ds);
        }

        // \todo: think about small-passing-by-clumps --> they are very weakly-bound with the big one while embeded in it
        // \todo: for the small-passing-by-clumps, I found the angles between their velocities and the center-of-mass velocity are fairly large, but are not the largest among all the particles. However, their |P_grav + E_k|/E_k are the lowest (between 0.2 to 0.3) among all the particles. This may be a hint to distinguish them.
        // \todo: according to the definition of Hill Radius, for those small clumps with super dense core (multiple particles at the same location), their Hill Radii are in fact pretty large (~ cell length), which means more particles belong to such small clumps. And those particles may have only low densities and are excluded from the begining.

        progIO->log_info << std::endl;
        progIO->Output(std::clog, progIO->log_info, __more_output, __all_processors);
        progIO->out_content << std::endl;
        progIO->Output(std::cout, progIO->out_content, __normal_output, __all_processors);
    }
    
};


/**************************************/
/********** Planetesimal Part *********/
/**************************************/

/*! \class template <int D> Planetesimal
 *  \brief data for one planetesimal structure */
template <int D>
class Planetesimal {
private:

public:
    /*! \var uint32_t peak_index
     *  \brief the index of the particle with the peak density */
    uint32_t peak_index {0};

    /*! \var std::vector<uint32_t> potential_subpeak_indices;
     *  \brief the indices of particles that may lead a sub-peak */
    std::vector<uint32_t> potential_subpeak_indices;

    /*! \var std::vector<uint32_t> indices
     *  \brief indices of particles that belongs to this planetesimals */
    std::vector<uint32_t> indices;

    /*! \var bool mask
     *  \brief its mask for selection process */
    bool mask {true};

    /*! \var double mass
     *  \brief total mass */
    double total_mass {0};

    /*! \var SmallVec<double, D> center_of_mass
     *  \brief center of mass */
    SmallVec<double, D> center_of_mass {0};

    /*! \var SmallVec<double, D> vel_com
     *  \brief velocity of the center of mass */
    SmallVec<double, D> vel_com {0};

    /*! \var std::vector<std::pair<uint32_t, double>> particles;
     *  \brief pairs of particle indices and their distances to the center of mass
     *  \todo: see if it is possible to use distances^2 to save computational cost */
    std::vector<std::pair<uint32_t, double>> particles;

    /* RL: After the hopping step, we need a criterion to determine whether or not to merge two particle groups (since there are usually many hopping chains within one dense clump). Again, two groups will be merged if one totally contains the other or one intersects with the other and they are gravitationally bound. This method requires the knowledge of "radius." The naivest choice is the distance from the COM to the farthest particle, which works fine for the case where planetesimals form sparely. But in order to generalize this code to situations where planetesimals form closely, we now move to a new definition of radius, one tenth radius (see definitions below)
     */
    
    /*! \var double outer_one10th_radius
     *  \brief the radius of the outermost particle that has a density > 1/10 the peak density */
    double outer_one10th_radius {0};
    
    /*! \var double inner_one10th_radius
     *  \brief the radius of the innermost particle that has a density < 1/10 the peak density */
    double inner_one10th_radius {0};
    
    /*! \var double one10th_radius
     *  \brief the average radius = (inner_one10th_radius + outer_one10th_radius) / 2. */
    double one10th_radius {0};
    
    /* RL: On the one hand, very-newly-formed clumps, early-stage-merging clumps, and clumps starting accreting smaller companions may deviate from spherical shapes a lot. This "inner_one10th_radius" will avoid over-estiamting their radii and prevent unnecessary "group-merging" during clump-finding. On the other hand, late-stage-mering clumps and clumps with tidally-disrupted companions/streams only deviate a little from spherical shapes. This "outer_one10th_radius" will identified them and their sub-structures as one individual clump and perform group-merging during clump-finding. Therefore, using the average radius "one10th_radius" enables this program to find planetesimals more accurately.
     * RL: Note that outer_one10th_radius is not necessary larger than inner_one10th_radius (e.g., two adjacent indices but with outer_one10th_radius being the smaller one). Furthermore, one clump identified by (outer_)one10th_radius only may correspond to multiple clumps identified by inner_one10th_radius only, which may provide a good practice for categorizing sub-clumps.
     */
    
    /*! \var double half_mass_radius
     *  \brief the radius that contains half of the total mass */
    //double half_mass_radius;

    /*! \var double two_sigma_mass_radius
     *  \brief the radius that contains 90% of the total mass */
    double two_sigma_mass_radius {0};

    std::vector<uint32_t> preys;
    
    /*! \var double Hill_radius
     *  \brief Hill radius */
    double Hill_radius {0};

    /*! \var double J
     *  \brief angular momentum vector */
    SmallVec<double, D> J {0};

    /*! \var double accumulated_J_in_quarter_Hill_radius
     *  \brief accumulated angular momentum vector within 0.25 * Hill_radius */
    SmallVec<double, D> accumulated_J_in_quarter_Hill_radius {0};

    /*! \fn void CalculateKinematicProperties(typename BHtree<D>::InternalParticle *particle_list)
     *  \brief calculate the total_mass, center_of_mass, and vel_com */
    void CalculateKinematicProperties(typename BHtree<D>::InternalParticle *particle_list) {
        total_mass = 0;
        center_of_mass = SmallVec<double, D>(0);
        vel_com = SmallVec<double, D>(0);
        for (auto it : indices) {
            total_mass += particle_list[it].mass;
            center_of_mass += particle_list[it].pos * particle_list[it].mass;
            vel_com += particle_list[it].vel * particle_list[it].mass;
        }
        center_of_mass /= total_mass;
        vel_com /= total_mass;
    }

    /*! \fn void SortParticles(typename BHtree<D>::InternalParticle *particle_list)
     *  \brief sort the particle indices by their distances to the center of mass */
    void SortParticles(typename BHtree<D>::InternalParticle *particle_list) {
        CalculateKinematicProperties(particle_list);
        for (auto it : indices) {
            particles.push_back(std::pair<uint32_t, double>(it, (particle_list[it].pos-center_of_mass).Norm()));
        }
        std::sort(particles.begin(), particles.end(), [](const std::pair<uint32_t, double> &a, const std::pair<uint32_t, double> &b) {
            return a.second < b.second;
        });
        
        auto tmp_it = indices.begin();
        for (auto it : particles) {
            *tmp_it = it.first;
            tmp_it++;
        }

        // calculate one tenth radius
        auto par_it = particles.begin();
        auto par_rit = particles.rbegin();
        //auto tmp_one10th_peak_density = particle_list[peak_index].new_density / 10.;
        auto tmp_one10th_peak_density = particle_list[indices[indices.size()/2]].new_density / 10.;
        //
        for (; par_it != particles.end(); ++par_it) {
            if (particle_list[par_it->first].new_density < tmp_one10th_peak_density) {
                break;
            }
        }
        if (par_it != particles.end()) {
            inner_one10th_radius = par_it->second;
        } else {
            inner_one10th_radius = particles.back().second;
        }
        for (; par_rit != particles.rend(); ++par_rit) {
            if (particle_list[par_rit->first].new_density > tmp_one10th_peak_density) {
                break;
            }
        }
        outer_one10th_radius = par_rit->second;
        one10th_radius = (inner_one10th_radius + outer_one10th_radius) / 2.; //*/

        // calculate half mass radius
        par_rit = particles.rbegin();
        double tmp_mass = 0., tenth_mass = total_mass * (1.0-0.95);
        for (; par_rit != particles.rend(); ++par_rit) {
            tmp_mass += particle_list[par_rit->first].mass;
            if (tmp_mass > tenth_mass) {
                break;
            }
        }
        two_sigma_mass_radius = par_rit->second;
        /*
        double half_mass = total_mass / 2.;
        for (; par_rit != particles.rend(); ++par_rit) {
            tmp_mass += particle_list[par_rit->first].mass;
            if (tmp_mass > half_mass) {
                break;
            }
        }
        half_mass_radius = par_rit->second; //*/

        // RL: here we use inner/outer as a way to identify elliptical/stream-like structures
        // \todo: use adaptive # of neighbors for hopping to find streams
        if ((inner_one10th_radius > 0 && outer_one10th_radius / inner_one10th_radius > 2.5)) {
            one10th_radius = std::min(outer_one10th_radius, two_sigma_mass_radius);
        } else {
            one10th_radius = std::min(one10th_radius, two_sigma_mass_radius);
        }
    }

    /*! \fn bool IsPositionDispersion2Large(double tolerance)
     *  \brief By examining the ratio of position dispersion over max(dist_to_COM), we may be able to remove some diffuse particle groups. From experiments, I found a concentrated particle group usually has such a ratio below 0.4, a stream-like particle group could have 0.5-ish, a randomly-distributed group will give ~0.6. */
    bool IsPositionDispersion2Large(double tolerance=0.55, double fraction=1.0) {
        // to avoid some outliers (a few particles really faraway)
        size_t num_poi = particles.size(); // poi = particle of interest
        /* RL: I've tried many different ways to eliminate extremes
         * none of them really works
        num_poi = std::round(particles.size() * fraction);
        if (particles.size() > 50) {
            num_poi = particles.size() - 5;
        } else {
            num_poi = particles.size();
        } //*/
        return std::accumulate(particles.begin(), particles.begin()+num_poi, 0.0, [](const double &a, const std::pair<uint32_t, double> &b) { return a + b.second*b.second; }) / num_poi > std::pow(tolerance * particles[num_poi-1].second, 2);
    }

    /*! \fn void MergeAnotherPlanetesimal(Planetesimal<D> carnivore, typename BHtree<D>::InternalParticle *particle_list)
     *  \brief take another planetesimal's data */
    void MergeAnotherPlanetesimal(Planetesimal<D> &carnivore, typename BHtree<D>::InternalParticle *particle_list) {
        if (particle_list[peak_index].new_density < particle_list[carnivore.peak_index].new_density) {
            potential_subpeak_indices.push_back(peak_index);
            peak_index = carnivore.peak_index;
        } else {
            potential_subpeak_indices.push_back(carnivore.peak_index);
        }
        double tmp_total_mass = total_mass + carnivore.total_mass;
        center_of_mass = (center_of_mass * total_mass + carnivore.center_of_mass * carnivore.total_mass) / tmp_total_mass;
        vel_com = (vel_com * total_mass + carnivore.vel_com * carnivore.total_mass) / tmp_total_mass;
        total_mass = tmp_total_mass;

        particles.insert(particles.end(), carnivore.particles.begin(), carnivore.particles.end());
        for (auto &it : particles) {
            it.second = (particle_list[it.first].pos - center_of_mass).Norm();
        }
        std::sort(particles.begin(), particles.end(), [](const std::pair<uint32_t, double> &a, const std::pair<uint32_t, double> &b) {
            return a.second < b.second;
        });

        // copy indices (do we really need it?)
        indices.resize(particles.size());
        auto tmp_it = indices.begin();
        for (auto it : particles) {
            *tmp_it = it.first;
            tmp_it++;
        }
        
        // recalculate one10th_radius
        auto par_it = particles.begin();
        auto par_rit = particles.rbegin();
        //auto tmp_one10th_peak_density = particle_list[peak_index].new_density / 10.;
        auto tmp_one10th_peak_density = particle_list[indices[indices.size()/2]].new_density / 10.;
        //
        for (; par_it != particles.end(); ++par_it) {
            if (particle_list[par_it->first].new_density < tmp_one10th_peak_density) {
                break;
            }
        }
        if (par_it != particles.end()) {
            inner_one10th_radius = par_it->second;
        } else {
            inner_one10th_radius = particles.back().second;
        }
        for (; par_rit != particles.rend(); ++par_rit) {
            if (particle_list[par_rit->first].new_density > tmp_one10th_peak_density) {
                break;
            }
        }
        outer_one10th_radius = par_rit->second;
        one10th_radius = (inner_one10th_radius + outer_one10th_radius) / 2.; //*/

        // calculate half mass radius
        par_rit = particles.rbegin();
        double tmp_mass = 0., tenth_mass = total_mass * (1.0-0.95);
        for (; par_rit != particles.rend(); ++par_rit) {
            tmp_mass += particle_list[par_rit->first].mass;
            if (tmp_mass > tenth_mass) {
                break;
            }
        }
        two_sigma_mass_radius = par_rit->second;
        /*
        double half_mass = total_mass / 2.;
        for (; par_rit != particles.rend(); ++par_rit) {
            tmp_mass += particle_list[par_rit->first].mass;
            if (tmp_mass > half_mass) {
                break;
            }
        }
        half_mass_radius = par_rit->second; //*/

        // RL: here we use inner/outer as a way to identify elliptical/stream-like structures
        // \todo: use adaptive # of neighbors for hopping to find streams
        if ((inner_one10th_radius > 0 && outer_one10th_radius / inner_one10th_radius > 2.5)) {
            one10th_radius = std::min(outer_one10th_radius, two_sigma_mass_radius);
        } else {
            one10th_radius = std::min(one10th_radius, two_sigma_mass_radius);
        }
    }

    /*! \fn void RemoveUnboundParticles(typename BHtree<D>::InternalParticle *particle_list)
     *  \brief remove unbound particles */
    void RemoveUnboundParticles(typename BHtree<D>::InternalParticle *particle_list) {
        auto i = particles.size()-1;
        double tmp_total_mass;
        SmallVec<double, D> tmp_center_of_mass, tmp_vel_com;
        double total_energy_over_mass_product;

        /* RL: test use
        double last_tmp_total_mass;
        tmp_total_mass = total_mass;
        last_tmp_total_mass = total_mass;
        while (i != 1) {
            tmp_total_mass -= particle_list[particles[i].first].mass; // only internal particles
            tmp_center_of_mass = (tmp_center_of_mass * last_tmp_total_mass - particle_list[particles[i].first].pos * particle_list[particles[i].first].mass) / tmp_total_mass;
            tmp_vel_com = (tmp_vel_com * last_tmp_total_mass - particle_list[particles[i].first].vel * particle_list[particles[i].first].mass) / tmp_total_mass;
            last_tmp_total_mass = tmp_total_mass;
        //*/

        //
        while (i != 0) {
            tmp_total_mass = total_mass - particle_list[particles[i].first].mass;
            tmp_center_of_mass = (center_of_mass * total_mass - particle_list[particles[i].first].pos * particle_list[particles[i].first].mass) / tmp_total_mass;
            tmp_vel_com = (vel_com * total_mass - particle_list[particles[i].first].vel * particle_list[particles[i].first].mass) / tmp_total_mass;
        //*/

            /* you may tweak them to achieve P_grav + 2 * E_k  0
            double P_grav, E_k;
            P_grav = -progIO->numerical_parameters.grav_constant / (tmp_center_of_mass - particle_list[particles[i].first].pos).Norm();
            E_k = 0.5 / total_mass * (tmp_vel_com - particle_list[particles[i].first].vel).Norm2();
             * we omit tmp_total_mass * particle_list[particles[i].first].mass
             * since P_grav and E_k all have it
             * ***** ***** ***** *****
             * I also tried to include "the tidal energy" and "vertically gravitational potential" like below
            double P_tidal, P_vg;
            P_tidal = + 1.5 / total_mass * pow(progIO->numerical_parameters.Omega, 2) * pow(tmp_center_of_mass[0] - particle_list[particles[i].first].pos[0], 2);
            P_vg = - 0.5 / total_mass * pow(progIO->numerical_parameters.Omega, 2) * pow(tmp_center_of_mass[2] - particle_list[particles[i].first].pos[2], 2);
             * However, such a tidal energy does not differentiate the prograde and retrograde orbits where the latter one is usually more stable in numerical studies.
             * And such a vertically gravitational potential may classify two clumps with a large vertical separation as bound binaries even though they don't feel each other.
             * Therefore, we should still check those standard energies.
             * An promising approach is that we check separation / Hill radius, and then check relative velocity / Hill velocity, where Hill velocity = Hill radius * Omega_K. The second step is equivalent to checking binding energy. If we plot Delta r / R_Hill versus Delta v / V_Hill and draw a line following 1/sqrt(dr), we can get a rough idea of how bound those possible binaries are.
             */

            total_energy_over_mass_product =
                    + 0.5 / total_mass * (tmp_vel_com - progIO->numerical_parameters.shear_vector*(tmp_center_of_mass[0]-particle_list[particles[i].first].pos[0]) - particle_list[particles[i].first].vel).Norm2()
                    - progIO->numerical_parameters.grav_constant / (tmp_center_of_mass - particle_list[particles[i].first].pos).Norm();
            if (total_energy_over_mass_product > 0) { // unbound
                std::swap(indices[i], indices.back());
                indices.pop_back();
                std::swap(particles[i], particles.back());
                particles.pop_back();
                vel_com = tmp_vel_com;
                center_of_mass = tmp_center_of_mass;
                total_mass = tmp_total_mass;
            }
            i--;
        }

        std::vector<std::pair<uint32_t, double>> tmp;
        particles.swap(tmp);
        SortParticles(particle_list);
    }

    /*! \fn void SearchBoundParticlesWithinHillRadius(BHtree<D> &tree, double density_threshold)
     *  \brief search potential bound particles within Hill radius */
    void SearchBoundParticlesWithinHillRadius(BHtree<D> &tree, double density_threshold) {
        // \todo should consider only particles within half the Hill radius
        uint32_t nearby_count = 0, idx = 0;
        tree.RecursiveBallSearchCount(center_of_mass, tree.root, Hill_radius, nearby_count);
        auto *nearby_indices = new uint32_t[nearby_count];
        tree.BallSearch(center_of_mass, Hill_radius, nearby_indices, nearby_count);

        double total_energy_over_mass_product;

        for (uint32_t i = 0; i != nearby_count; i++) {
            idx = nearby_indices[i];
            if (tree.particle_list[idx].in_clump_flag || tree.particle_list[idx].new_density < density_threshold) {
                continue;
            } else {
                auto tmp_total_mass = total_mass + tree.particle_list[idx].mass;
                // again we omit total_mass * tree.particle_list[idx].mass since P_grav and E_k all have it
                total_energy_over_mass_product = + 0.5 / tmp_total_mass * (vel_com - progIO->numerical_parameters.shear_vector*(center_of_mass[0]-tree.particle_list[idx].pos[0]) - tree.particle_list[idx].vel).Norm2()
                                                 - progIO->numerical_parameters.grav_constant / (center_of_mass - tree.particle_list[idx].pos).Norm();
                if (total_energy_over_mass_product < 0) { // bound
                    indices.push_back(idx);
                    total_mass = tmp_total_mass;
                    center_of_mass = (center_of_mass * total_mass + tree.particle_list[idx].mass * tree.particle_list[idx].pos) / total_mass;
                    vel_com = (total_mass * vel_com + tree.particle_list[idx].mass * tree.particle_list[idx].vel) / total_mass;
                    tree.particle_list[idx].in_clump_flag = true;
                }
            }
        }
        if (indices.size() > particles.size()) {
            std::vector<std::pair<uint32_t, double>> tmp;
            particles.swap(tmp);
            SortParticles(tree.particle_list);
        }
        delete[] nearby_indices;
    }

    /*! \fn void CalculateHillRadius()
     *  \brief calculate the Hill radius */
    void CalculateHillRadius() {
        Hill_radius = pow(total_mass * progIO->numerical_parameters.grav_constant / 3. / progIO->numerical_parameters.Omega / progIO->numerical_parameters.Omega, 1./3.);
    }

    /*! \fn void CalculateAngularMomentum(BHtree<D> &tree)
     *  \brief calculate the angular momentum in Hill units */
    void CalculateAngularMomentum(BHtree<D> &tree) {
        SmallVec<double, D> tmp_j {0};
        SmallVec<double, D> tmp_dr {0};
        SmallVec<double, D> tmp_dv {0};
        double quarter_Hill_radius = 0.25 * Hill_radius;
        double Hill_units {total_mass * Hill_radius * Hill_radius * progIO->numerical_parameters.Omega};
        SmallVec<double, D> shear_vector (0., progIO->numerical_parameters.q * progIO->numerical_parameters.Omega, 0.);

        for (auto it : particles) {
            tmp_dr = tree.particle_list[it.first].pos - center_of_mass;
            tmp_dv = tree.particle_list[it.first].vel - shear_vector * tree.particle_list[it.first].pos[0]
                     - (                vel_com - shear_vector * center_of_mass[0]);
            tmp_j = tmp_dr.Cross(tmp_dv);
            tmp_j += progIO->numerical_parameters.Omega
                     * SmallVec<double, 3>(-tmp_dr[0]*tmp_dr[2],
                                           -tmp_dr[1]*tmp_dr[2],
                                           (tmp_dr[0]*tmp_dr[0] + tmp_dr[1]*tmp_dr[1]));
            J += tree.particle_list[it.first].mass * tmp_j;
            if (it.second < quarter_Hill_radius) {
                accumulated_J_in_quarter_Hill_radius += tree.particle_list[it.first].mass * tmp_j;
            }
        }

        // normalizing angular momentum in units of (Mass * R_Hill^2 * Omega)
        J /= Hill_units;
        accumulated_J_in_quarter_Hill_radius /= Hill_units;
    }

    /*! \fn void void CalculateCumulativeAngularMomentum(BHtree<D> &tree, std::ofstream &f)
     *  \brief calculate the cumulative angular momentum inside out in Hill units */
    void CalculateCumulativeAngularMomentum(BHtree<D> &tree, uint32_t id_peak, std::ofstream &f) {
        SmallVec<double, D> tmp_j {0};
        SmallVec<double, D> tmp_dr {0};
        SmallVec<double, D> tmp_dv {0};
        double quarter_Hill_radius = 0.25 * Hill_radius;
        double Hill_units {total_mass * Hill_radius * Hill_radius * progIO->numerical_parameters.Omega};
        SmallVec<double, D> shear_vector (0., progIO->numerical_parameters.q * progIO->numerical_parameters.Omega, 0.);
        std::vector<double> accumulated_m (indices.size(), 0);
        std::vector<double> accumulated_Jz (indices.size(), 0);
        std::vector<std::pair<double, double>> xyJz (indices.size(), std::pair<double, double>(0, 0));

        size_t idx = 0;
        for (auto it : indices) {
            tmp_dr = tree.particle_list[it].pos - center_of_mass;
            tmp_dv = tree.particle_list[it].vel - shear_vector * tree.particle_list[it].pos[0]
                     - (vel_com - shear_vector * center_of_mass[0]);
            tmp_j = tmp_dr.Cross(tmp_dv);

            tmp_j[0] -= progIO->numerical_parameters.Omega * tmp_dr[0] * tmp_dr[2];
            tmp_j[1] -= progIO->numerical_parameters.Omega * tmp_dr[1] * tmp_dr[2];
            tmp_j[2] += progIO->numerical_parameters.Omega * (tmp_dr[0]*tmp_dr[0] + tmp_dr[1]*tmp_dr[1]);

            if (particles[idx].second < quarter_Hill_radius) {
                accumulated_J_in_quarter_Hill_radius += tree.particle_list[it].mass * tmp_j;
            }

            J += tree.particle_list[it].mass * tmp_j;
            accumulated_Jz[idx] = tree.particle_list[it].mass * tmp_j[2];
            accumulated_m[idx] = tree.particle_list[it].mass;
            if (idx > 0) {
                accumulated_Jz[idx] += accumulated_Jz[idx-1];
                accumulated_m[idx] += accumulated_m[idx-1];
            }
            xyJz[idx].first = std::sqrt(tmp_dr[0]*tmp_dr[0]+tmp_dr[1]*tmp_dr[1]);
            xyJz[idx].second = tree.particle_list[it].mass * tmp_j[2];
            idx++;
        }
        // resort Jz following the order of cylindrical r
        std::sort(xyJz.begin(), xyJz.end(), [](const std::pair<double, double> &a, const std::pair<double, double> &b) {
            return a.first < b.first;
        });

        //f << std::defaultfloat << indices.size() << std::endl; // std::defaultfloat not supported in GCC 4.9.4
        f.unsetf(std::ios_base::floatfield);
        f << id_peak << ' ' << indices.size() << std::endl;
        idx = 0;
        for (auto it : particles) {
            if (idx > 0) {
                xyJz[idx].second += xyJz[idx-1].second;
            }
            f << std::scientific << std::setprecision(12) << std::setw(20) << it.second/Hill_radius << std::setw(20) << accumulated_m[idx]/total_mass << std::setw(20) << accumulated_Jz[idx]/Hill_units << std::setw(20) << xyJz[idx].first/Hill_radius  << std::setw(20) << xyJz[idx].second/Hill_units << std::endl;
            idx++;
        }

        // normalizing angular momentum in units of (Mass * R_Hill^2 * Omega)
        J /= Hill_units;
        accumulated_J_in_quarter_Hill_radius /= Hill_units;
    }
};

/*! \class template <int D> PlanetesimalList
 *  \brief store a list of planetesimal structures
 * @tparam D dimension of data */
template <int D>
class PlanetesimalList {
private:

public:
    /*! \var uint32_t num_planetesimals
     *  \breif the number of planetesimals found */
    uint32_t num_planetesimals {0};

    /*! \var double density_threshold
     *  \brief the density criterion for choosing particles (~ rho_crit for outer rims of planetesimals) */
    double density_threshold {0};

    /*! \var double clump_mass_threshold
     *  \brief the mass criterion for removing weak clumps */
    double clump_mass_threshold {0};

    /*! \var double reliable_fraction_for_diffuse
     *  \brief the reliable fraction of inner particles used in IsPositionDispersion2Large(); lower hydro-resolution needs lower value
    double reliable_fraction_for_diffuse {1.0}; //*/

    /*! \var double clump_diffuse_threshold
     *  \brief the diffuse threshold for a particle group; input for IsPositionDispersion2Large() */
    double clump_diffuse_threshold {0.55};

    /*! \var double Hill_fraction_for_merge
     *  \brief the fraction of Hill radius used when further merging primitive clumps */
    double Hill_fraction_for_merge {0.25};

    /*! \var double peak_density_threshold
     *  \brief the density criterion for peak densities of planetesimals (rho_crit for the particle with the highest density) */
    double peak_density_threshold {0};

    /*! \var std::vector<std::pair<uint32_t, double>> peaks_and_masses
     *  \brief a ordered vector of pairs recording the peak_indices (key to the map below) and masses */
    std::vector<std::pair<uint32_t, double>> peaks_and_masses;

    /*! \var std::map<uint32_t, Planetesimal>
     *  \brief a map pairing the peak particle indices with planetesimal structures
     *  This STL container, unordered_map, is suitable for quick-access by key values (with its internal hash table), which is very useful during grouping particles by those peak-density-particles since many partcle-hopping-chains belongs to one ultimate peak-density-particle. In addition, it uses non-contiguous memory which fits our needs because the object Planetesimal requires a lot of changes in size during unbinding processes.
     *  RL: the order of the iterations through an unordered_map is not guaranteed to be the same. So don't count on it during the merging process. It could be "A" eats "B" or "B" eats "A", which depends on who comes first. I switched to ordered_map for better cross-platform/compiler consistency. Keep in mind that iterating through an ordered_map always follows the non-descending order of keys. */
    std::map<uint32_t, Planetesimal<D>> planetesimals;

    /*! \var ParticleSet<D> clump_set
     *  \brief like particle_set, used in building clump_tree */
    ParticleSet<D> clump_set;

    /*! \var BHtree<D> clump_tree
     *  \brief put clumps into tree structures for easier manipulations */
    BHtree<D> clump_tree;

    /*! \fn bool IsGravitationallyBound(Planetesimal<D> &p1, Planetesimal<D> &p2)
     *  \brief determine if two planeteismals structures are gravitationally bound */
    bool IsGravitationallyBound(const Planetesimal<D> &p1, const Planetesimal<D> &p2) {
        // first calculate the mutual gravitational potential energy total kinematic energy, mutual tidal energy and vertically gravitational potential
        // omitting p1.total_mass * p2.total_mass since both P_grav and E_k have it
        double total_mass = p1.total_mass + p2.total_mass;
        double P_grav = - progIO->numerical_parameters.grav_constant / (p1.center_of_mass - p2.center_of_mass).Norm();
        double E_k = 0.5 / total_mass * (p1.vel_com - progIO->numerical_parameters.shear_vector * (p1.center_of_mass[0] - p2.center_of_mass[0]) - p2.vel_com).Norm2();
        return P_grav + E_k < 0.;
    }

    /*! \fn bool template <class T> IsSaddlePointDeepEnough(const DataSet<T, D> &ds, Planetesimal<D> &p1, Planetesimal<D> &p2)
     *  \brief determine if two planetesimals structures have a deep saddle point in between
     *  If so, these two clumps won't be merged. The criteria for a deep saddle point is under 2.5 times the density threshold (from HOP's experience, ref: https://www.cfa.harvard.edu/~deisenst/hop/hop_doc.html) */
    template <class T, class F>
    bool IsSaddlePointDeepEnough(DataSet<T, D> &ds, F DensityKernel, const Planetesimal<D> &p1, const Planetesimal<D> &p2, double saddle_threshold=2.5) {
        auto r12 = p2.center_of_mass - p1.center_of_mass;
        r12 *= p1.one10th_radius / (p1.one10th_radius + p2.one10th_radius);
        auto possible_saddle_point = p1.center_of_mass + r12;
        
        double radius_Kth_NN = 0;
        auto *indices = new uint32_t[progIO->numerical_parameters.num_neighbors_in_knn_search];
        std::vector<std::pair<int, double>> local_heaps;
        ds.tree.KNN_Search_OpenMP(possible_saddle_point, progIO->numerical_parameters.num_neighbors_in_knn_search, radius_Kth_NN, indices, local_heaps, std::min(p1.one10th_radius, p2.one10th_radius));
        double saddle_density = DensityKernel(ds, radius_Kth_NN, 0, indices, local_heaps);
        delete [] indices;
        indices = nullptr;
        return saddle_density < saddle_threshold * ds.planetesimal_list.density_threshold;
    }

    template <class T, class F>
    bool IsHillSaddlePointDeepEnough(DataSet<T, D> &ds, F DensityKernel, const Planetesimal<D> &p1, const Planetesimal<D> &p2, double saddle_threshold=2.5) {
        auto r12 = p2.center_of_mass - p1.center_of_mass;
        r12 *= p1.Hill_radius / (p1.Hill_radius + p2.Hill_radius);
        auto possible_saddle_point = p1.center_of_mass + r12;

        double radius_Kth_NN = 0;
        auto *indices = new uint32_t[progIO->numerical_parameters.num_neighbors_in_knn_search];
        std::vector<std::pair<int, double>> local_heaps;
        ds.tree.KNN_Search_OpenMP(possible_saddle_point, progIO->numerical_parameters.num_neighbors_in_knn_search, radius_Kth_NN, indices, local_heaps, std::min(p1.Hill_radius, p2.Hill_radius) / Hill_fraction_for_merge);
        double saddle_density = DensityKernel(ds, radius_Kth_NN, 0, indices, local_heaps);
        delete [] indices;
        indices = nullptr;
        return saddle_density < saddle_threshold * ds.planetesimal_list.density_threshold;
    }

    template <class T>
    bool IsPhaseSpaceDistanceWithinTenSigma(DataSet<T, D> &ds, const Planetesimal<D> &p1, const Planetesimal<D> &p2) {
        const Planetesimal<D> &small_p = (p1.total_mass > p2.total_mass)? p2 : p1;
        const Planetesimal<D> &large_p = (p1.total_mass > p2.total_mass)? p1 : p2;
        double sigma2_pos = 0, sigma2_vel = 0, small_n = small_p.particles.size();
        double phase_dist2 = 0;
        typename BHtree<D>::InternalParticle *p;
        for (auto &par: small_p.particles) {
            sigma2_pos += par.second*par.second;
            p = &ds.tree.particle_list[par.first];
            sigma2_vel += (p->vel - small_p.vel_com - progIO->numerical_parameters.shear_vector * (p->pos[0] - small_p.center_of_mass[0]) ).Norm2();
            //sigma2_vel += (p->vel - small_p.vel_com).Norm2();
        }
        sigma2_pos /= small_n;
        sigma2_vel /= small_n;

        phase_dist2 = (large_p.center_of_mass - small_p.center_of_mass).Norm2() / (sigma2_pos / small_n)
                     +(large_p.vel_com - small_p.vel_com - progIO->numerical_parameters.shear_vector
                     * (large_p.center_of_mass[0] - small_p.center_of_mass[0])).Norm2() / (sigma2_vel / small_n);
        return phase_dist2 < 200.; // (10*sqrt(2))^2
    }

    /*! \fn void WriteBasicResults(int loop_count)
     *  \brief write basic peak-finding results into one file */
    void WriteBasicResults(int loop_count) {
        if (loop_count == mpi->loop_begin) {
            // write header
            progIO->out_content << std::setw(progIO->width) << std::setfill(' ') << "#time" << std::setw(progIO->width) << "N_peak" << std::setw(progIO->width) << std::setfill(' ') << "M_p,max,code" << std::setw(progIO->width) << std::setfill(' ') << "M_p,tot,code" << std::setw(progIO->width) << std::setfill(' ') << "M_p,max,frac" << std::setw(progIO->width) << std::setfill(' ') << "M_p,tot,frac" << std::setw(progIO->width) << std::setfill(' ') << "M_p,max,Ceres" << std::setw(progIO->width) << std::setfill(' ') << "M_p,tot,Ceres" << std::endl;
            mpi->WriteSingleFile(mpi->result_files[mpi->file_pos[progIO->file_name.planetesimals_file]], progIO->out_content, __master_only);
        }

        progIO->out_content <<  std::setw(progIO->width) << std::setfill(' ') << std::scientific << progIO->physical_quantities[loop_count].time;
        progIO->out_content << std::setw(progIO->width) << std::setfill(' ') << num_planetesimals;
        for (auto &it : planetesimals) {
            progIO->physical_quantities[loop_count].max_planetesimal_mass = std::max(progIO->physical_quantities[loop_count].max_planetesimal_mass, it.second.total_mass);
            progIO->physical_quantities[loop_count].total_planetesimal_mass += it.second.total_mass;
        }
        progIO->out_content << std::setw(progIO->width) << std::setfill(' ') << progIO->physical_quantities[loop_count].max_planetesimal_mass << std::setw(progIO->width) << std::setfill(' ') << progIO->physical_quantities[loop_count].total_planetesimal_mass;
        progIO->out_content << std::setw(progIO->width) << std::setfill(' ') << progIO->physical_quantities[loop_count].max_planetesimal_mass / progIO->numerical_parameters.mass_total_code_units << std::setw(progIO->width) << std::setfill(' ') << progIO->physical_quantities[loop_count].total_planetesimal_mass / progIO->numerical_parameters.mass_total_code_units;
        progIO->out_content << std::setw(progIO->width) << std::setfill(' ') << progIO->physical_quantities[loop_count].max_planetesimal_mass * progIO->numerical_parameters.mass_physical / progIO->numerical_parameters.mass_ceres << std::setw(progIO->width) << std::setfill(' ') << progIO->physical_quantities[loop_count].total_planetesimal_mass * progIO->numerical_parameters.mass_physical / progIO->numerical_parameters.mass_ceres;
        progIO->out_content << std::endl;
        mpi->WriteSingleFile(mpi->result_files[mpi->file_pos[progIO->file_name.planetesimals_file]], progIO->out_content, __all_processors);
    }

    /*! \fn void OutputPlanetesimalsInfo(int loop_count, DataSet<T, D> &ds)
     *  \brief output the kinematic info of planetesimals */
    template <class T>
    void OutputPlanetesimalsInfo(int loop_count, DataSet<T, D> &ds) {
        std::ofstream file_planetesimals;
        std::ostringstream tmp_ss;
        tmp_ss << std::setprecision(3) << std::fixed << std::setw(7) << std::setfill('0') << progIO->physical_quantities[loop_count].time;
        std::string tmp_file_name;
        if (progIO->file_name.output_file_path.find_last_of('/') != std::string::npos) {
            tmp_file_name = progIO->file_name.output_file_path.substr(0, progIO->file_name.output_file_path.find_last_of('/')) + std::string("/peaks_at_") + tmp_ss.str() + std::string(".txt");
        } else {
            // operates under current directory
            tmp_file_name = std::string("peaks_at_") + tmp_ss.str() + std::string(".txt");
        }

        file_planetesimals.open(tmp_file_name, std::ofstream::out);
        if (!(file_planetesimals.is_open())) {
            progIO->error_message << "Error: Failed to open file " << tmp_file_name << " due to " << std::strerror(errno) << std::endl;
            progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
        }

        std::stringstream tmp_ss_id;
        if (progIO->flags.save_clumps_flag) {
            char mkdir_cmd[500] = "mkdir -p ParList.";
            tmp_ss_id << std::setw(4) << std::setfill('0') << loop_count * progIO->interval + progIO->start_num;
            std::strcat(mkdir_cmd, tmp_ss_id.str().c_str());
            if (std::system(mkdir_cmd) == -1) {
                progIO->error_message << "Error: Failed to execute: " << mkdir_cmd << std::endl;
                progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
            }
        }

        /* RL: old output with many diagnostic info
        file_planetesimals << std::setw(24) << "# center_of_mass[x]" << std::setw(24) << "center_of_mass[y]" << std::setw(24) << "center_of_mass[z]" << std::setw(10) << "peak_id" << std::setw(24) << "<-its_rho" << std::setw(10) << "Npar" << std::setw(24) << "total_mass" << std::setw(24) << "Hill_radius" << std::setw(24) << "max(dist2COM)" << std::setw(24) << "total_mass/volume_by_max_dist2COM" << std::endl;
        for (auto peak : peaks_and_masses) {
            auto it = planetesimals.find(peak.first);
            file_planetesimals << std::setprecision(16) << std::setw(24) << it->second.center_of_mass[0] << std::setw(24) << it->second.center_of_mass[1] << std::setw(24) << it->second.center_of_mass[2] << std::setw(10) << it->first << std::setw(24) << ds.tree.particle_list[it->first].new_density;
            auto tmp_num_particles = it->second.indices.size();
            for (auto item : it->second.indices) {
                auto tmp_it = ds.tree.sink_particle_indices.find(item);
                if (tmp_it != ds.tree.sink_particle_indices.end()) {
                    tmp_num_particles += tmp_it->second.size()-1;
                }
            }
            file_planetesimals << std::setw(10) << tmp_num_particles << std::setw(24) << it->second.total_mass << std::setw(24) << it->second.Hill_radius << std::setw(24) << it->second.particles.back().second << std::setw(24) << it->second.total_mass / (progIO->numerical_parameters.four_PI_over_three * pow(it->second.particles.back().second, 3));
            file_planetesimals << std::endl;
        }
        //*/

        file_planetesimals << std::setw(12) << "#peak_id"
                           << std::setw(12) << "Npar"
                           << std::setw(24) << "total_mass"
                           << std::setw(24) << "Hill_radius"
                           << std::setw(24) << "center_of_mass[x]"
                           << std::setw(24) << "center_of_mass[y]"
                           << std::setw(24) << "center_of_mass[z]"
                           << std::setw(24) << "vel_COM[x]"
                           << std::setw(24) << "vel_COM[y]"
                           << std::setw(24) << "vel_COM[z]"
                           << std::setw(24) << "J[x]"
                           << std::setw(24) << "J[y]"
                           << std::setw(24) << "J[z]/(M R_H^2 Omega)"
                           << std::setw(24) << "geo_mean_offset[x]"
                           << std::setw(24) << "geo_mean_offset[y]"
                           << std::setw(24) << "geo_mean_offset[z]"
                           << std::setw(24) << "median_offset[x]"
                           << std::setw(24) << "median_offset[y]"
                           << std::setw(24) << "median_offset[z]"
                           //<< std::setw(24) << "accum_J_RH/4[x]"
                           //<< std::setw(24) << "accum_J_RH/4[y]"
                           //<< std::setw(24) << "accum_J_RH/4[z]"
                           << std::endl;
        for (auto peak : peaks_and_masses) {
            auto it = planetesimals.find(peak.first);
            auto tmp_num_particles = it->second.indices.size();
            sn::dvec geo_mean_offset, median_offset;
            //double ath_density = 0, peak_density = ds.tree.particle_list[it->first].new_density;
            for (auto item : it->second.indices) {
                auto tmp_it = ds.tree.sink_particle_indices.find(item);
                if (tmp_it != ds.tree.sink_particle_indices.end()) {
                    tmp_num_particles += tmp_it->second.size()-1;
                }
                //ath_density = MaxOf(ds.tree.particle_list[item].ath_density, ath_density);
            }
            std::vector<std::vector<double>> offset;
            offset.resize(3);
            for (auto &item : offset) {
                item.resize(tmp_num_particles);
            }
            uint32_t idx = 0;
            for (auto item : it->second.indices) {
                auto tmp_it = ds.tree.sink_particle_indices.find(item);
                if (tmp_it != ds.tree.sink_particle_indices.end()) {
                    for (auto &tmp_sink_it : tmp_it->second) {
                        for (size_t d = 0; d < 3; d++) {
                            offset[d][idx] = tmp_sink_it.pos[d] - it->second.center_of_mass[d];
                        }
                        idx++;
                    }
                } else {
                    for (size_t d = 0; d < 3; d++) {
                        offset[d][idx] = ds.tree.particle_list[item].pos[d] - it->second.center_of_mass[d];
                    }
                    idx++;
                }
            }

            // perform angle correction based on J before calculating the offset
            auto tmp_J = sn::dvec(it->second.J[0], it->second.J[1], it->second.J[2]);
            double tmp_theta = std::acos(it->second.J[2] / tmp_J.Norm());
            double tmp_phi = std::atan2(it->second.J[1], it->second.J[0]);

            double rot_z[3][3] = {{std::cos(tmp_phi), -std::sin(tmp_phi), 0},
                                  {std::sin(tmp_phi),  std::cos(tmp_phi), 0},
                                  {                0,                  0, 1}};
            double rot_y[3][3] = {{ std::cos(tmp_theta), 0, std::sin(tmp_theta)},
                                  {                   0, 1,                   0},
                                  {-std::sin(tmp_theta), 0, std::cos(tmp_theta)}};
            double rot[3][3];
            for (size_t d1 = 0; d1 < 3; d1++) {
                for (size_t d2 = 0; d2 < 3; d2++) {
                    rot[d1][d2] = 0;
                    for (size_t d = 0; d < 3; d++) {
                        rot[d1][d2] += rot_z[d1][d] * rot_y[d][d2];
                    }
                }
            }

            for (idx = 0; idx < tmp_num_particles; idx++) {
                double tmp_offset[3] = {offset[0][idx], offset[1][idx], offset[2][idx]};
                for (size_t d = 0; d < 3; d++) {
                    offset[d][idx] = 0;
                    for (size_t d1 = 0; d1 < 3; d1++) {
                        offset[d][idx] += tmp_offset[d1] * rot[d1][d];
                    }
                }
            }

            auto half_size_offset = tmp_num_particles / 2;
            bool is_even = !(tmp_num_particles & 1);
            for (size_t d = 0; d < 3; d++) {
                for (auto &item : offset[d]) item = std::abs(item);
                geo_mean_offset[d] = std::pow(10.0, std::accumulate(offset[d].begin(), offset[d].end(), 0., [](const double &a, const double &b) {
                    if (b < 1e-32) {
                        return a - 32;
                    } else {
                        return a + std::log10(b);
                    }
                }) / tmp_num_particles);
                std::nth_element(offset[d].begin(), offset[d].begin() + half_size_offset, offset[d].end());
                median_offset[d] = offset[d][half_size_offset];
                if (is_even) {
                    median_offset[d] = (median_offset[d] + *std::max_element(offset[d].begin(), offset[d].begin() + half_size_offset)) / 2.0;
                }
            }
            /* RL: output the cumulative angular momentum inside out
            std::ofstream cumulative_file("cumulative_Jz.txt", std::ofstream::out | std::ofstream::app);
            if (!cumulative_file.is_open()) {
                std::cout << "Fail to open cumulative_Jz.txt" << std::endl;
            }
            it->second.CalculateCumulativeAngularMomentum(ds.tree, it->first, cumulative_file);
            cumulative_file.close();
            //*/

            if (progIO->flags.save_clumps_flag) {
                OutputSinglePlanetesimal(std::string("ParList.")+tmp_ss_id.str()+std::string("/")+std::to_string(it->first)+std::string(".txt"), it->first, ds);
            }
            file_planetesimals << std::setw(12) << it->first << std::setw(12) << tmp_num_particles
                               << std::setprecision(16)
                               << std::setw(24) << it->second.total_mass
                               << std::setw(24) << it->second.Hill_radius
                               << std::setw(24) << it->second.center_of_mass[0]
                               << std::setw(24) << it->second.center_of_mass[1]
                               << std::setw(24) << it->second.center_of_mass[2]
                               << std::setw(24) << it->second.vel_com[0]
                               << std::setw(24) << it->second.vel_com[1]
                               << std::setw(24) << it->second.vel_com[2]
                               << std::setw(24) << it->second.J[0]
                               << std::setw(24) << it->second.J[1]
                               << std::setw(24) << it->second.J[2]
                               << std::setw(24) << geo_mean_offset[0]
                               << std::setw(24) << geo_mean_offset[1]
                               << std::setw(24) << geo_mean_offset[2]
                               << std::setw(24) << median_offset[0]
                               << std::setw(24) << median_offset[1]
                               << std::setw(24) << median_offset[2]
                               //<< std::setw(24) << peak_density
                               //<< std::setw(24) << ath_density
                               //<< std::setw(24) << it->second.accumulated_J_in_quarter_Hill_radius[0]
                               //<< std::setw(24) << it->second.accumulated_J_in_quarter_Hill_radius[1]
                               //<< std::setw(24) << it->second.accumulated_J_in_quarter_Hill_radius[2]
                               << std::endl;
            // RL: output the particle list of every clump
            //OutputSinglePlanetesimal("clump_parlist/"+std::to_string(it->first)+".txt", it->first, ds);
        }
        file_planetesimals.close(); //*/


        // todo: after revealing sub-clumps, we need to output them hierarchically
        /* RL: output all clumps and their particle lists in one file
        std::ofstream file_parlist;
        tmp_file_name = progIO->file_name.output_file_path.substr(0, progIO->file_name.output_file_path.find_last_of('/'))+std::string("/parlist_")+tmp_ss.str()+std::string(".txt");
        file_parlist.open(tmp_file_name);
        if (!(file_parlist.is_open())) {
            progIO->error_message << "Error: Failed to open file_parlist due to " << std::strerror(errno) << std::endl;
            progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
        }

        file_parlist << planetesimals.size() << std::endl; // (1), # of planetesimals
        for (auto &it : planetesimals) {
            auto tmp_num_particles = it.second.indices.size();
            for (auto item : it.second.indices) {
                auto tmp_it = ds.tree.sink_particle_indices.find(item);
                if (tmp_it != ds.tree.sink_particle_indices.end()) {
                    tmp_num_particles += tmp_it->second.size()-1;
                }
            }
            file_parlist << tmp_num_particles << std::endl; // (2) # of particles
            for (auto item : it.second.indices) {
                auto tmp_it = ds.tree.sink_particle_indices.find(item);
                if (tmp_it != ds.tree.sink_particle_indices.end()) {
                    for (auto sink_par : tmp_it->second) {
                        file_parlist << sink_par.original_id << std::endl;
                    }
                } else {
                    file_parlist << ds.tree.particle_list[item].original_id << std::endl;
                }
            }
        }
        file_parlist.close();
        //*/
    }

    /*! \fn template<class T> void OutputSinglePlanetesimal(std::string file_name, uint32_t peak_id, DataSet<T, D> &ds)
     *  \brief output all particles in one planetesimal */
    template<class T>
    void OutputSinglePlanetesimal(const std::string &file_name, uint32_t peak_id, DataSet<T, D> &ds, size_t precision=16) {
        auto width = precision + 8;
        std::ofstream file_single_clump;
        auto search_it = planetesimals.find(peak_id);
        if (search_it != planetesimals.end()) {
            file_single_clump.open(file_name);
            if (!(file_single_clump.is_open())) {
                progIO->error_message << "Error: Failed to open file: " << file_name << ". But we proceed." << std::endl;
                progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
            }
            uint32_t skip_for_sub_sampling = 0;
            for (auto it : search_it->second.particles) {
                if (skip_for_sub_sampling > 0) {
                    skip_for_sub_sampling--;
                    continue;
                } else {
                    skip_for_sub_sampling = progIO->save_clump_sampling_rate - 1;
                }
                file_single_clump.unsetf(std::ios_base::floatfield);
                file_single_clump << std::setw(precision) << ds.tree.particle_list[it.first].original_id;
                file_single_clump << std::scientific;
                for (int i = 0; i != D; i++) {
                    file_single_clump << std::setprecision(precision) << std::setw(width) << ds.tree.particle_list[it.first].pos[i];
                }
                for (int i = 0; i != D; i++) {
                    file_single_clump << std::setprecision(precision) << std::setw(width) << ds.tree.particle_list[it.first].vel[i];
                }
                file_single_clump << std::setprecision(precision) << std::setw(width) << ds.tree.particle_list[it.first].mass << std::endl;
            }
            file_single_clump.close();
        }
    }

    /*! \fn void OutputParticlesByIndices(const std::string &file_name, const std::vector<uint32_t> &indices, const ParticleSet<D> &particle_set)
     *  \brief output particles by a given index list (indices in particle_set); this only works if particle_set.particles is sorted and all the particles have contiguous id numbers (so id = index) */
    void OutputParticlesByIndices(const std::string &file_name, const std::vector<uint32_t> &indices, const ParticleSet<D> &particle_set) {
        std::ofstream file_particles;
        file_particles.open(file_name);
        if (!(file_particles.is_open())) {
            progIO->error_message << "Error: Failed to open file: " << file_name << ". But we proceed. " << std::endl;
            progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
        }

        for (auto it : indices) {
            file_particles.unsetf(std::ios_base::floatfield);
            file_particles << std::setw(16) << particle_set[it].id;
            file_particles << std::scientific;
            for (int i = 0; i != D; i++) {
                file_particles << std::setprecision(16) << std::setw(24) << particle_set[it].pos[i];
            }
            for (int i = 0; i != D; i++) {
                file_particles << std::setprecision(16) << std::setw(24) << particle_set[it].vel[i];
            }
            file_particles <<  std::setprecision(16) << std::setw(24) << progIO->numerical_parameters.mass_per_particle[particle_set[it].property_index] << std::endl;
        }
        file_particles.close();
    }

    /*! \fn void BuildClumpTree(sn::dvec root_center, double half_width, double &max_radius)
     *  \brief put planetesimals into tree structures */
    void BuildClumpTree(sn::dvec &root_center, double half_width, double &max_radius, bool Hill=false) {
        clump_set.Reset();
        clump_set.num_total_particles = static_cast<uint32_t>(num_planetesimals);
        clump_set.num_particles = static_cast<uint32_t>(num_planetesimals);
        clump_set.AllocateSpace(clump_set.num_total_particles);

        uint32_t tmp_id = 0;
        Particle<D> *p;
        for (auto &it : planetesimals) {
            if (it.second.mask) {
                p = &clump_set[tmp_id];
                p->pos = it.second.center_of_mass;
                p->vel = it.second.vel_com;
                // Previously, p->id = it.first. However, using it.first will cause EXC_BAD_ACCESS in MakeSinkParticle() due to the particle_set[it_par]
                p->id = tmp_id;
                p->density = it.first; // now we use density to store the peak id
                p->property_index = 0; // has no meaning, N.B., tmp_ds.tree does not have correct mass data
                tmp_id++;
                if (Hill) {
                    max_radius = MaxOf(it.second.Hill_radius * Hill_fraction_for_merge, max_radius);
                } else {
                    max_radius = MaxOf(it.second.one10th_radius, max_radius);
                }
            }
        }

        // build a clump tree in an overlapping spatial volume
        clump_tree.root_center = root_center;
        clump_tree.half_width = half_width;
        clump_tree.BuildTree(progIO->numerical_parameters, clump_set, true, false);
        if (clump_tree.sink_particle_indices.size() > 0) {
            progIO->log_info << "(Warning: got " << clump_tree.sink_particle_indices.size() << " sink clumps while building clump tree. Proceed for now.) ";
        }
    }

};


/**********************************/
/********** DataSet Part **********/
/**********************************/

/*! \class template <class T, int D> DataSet
 *  \brief data set that encapsulates other data classes
 *  \tparam T type of data
 *  \tparam D dimension of data */
template <class T, int D>
class DataSet {
private:
    
public:
    
    /*! \var VtkData<T, D> vtk_data
     *  \brief data from vtk files */
    VtkData<T, D> vtk_data;
    
    /*! \var ParticleSet<D> particle_set
     *  \brief data about particles
     *  RL: I found even for 2D simulation, LIS output 3D data (with certain vars=0) */
    ParticleSet<3> particle_set;
    
    /*! \var BHtree<dim> tree
     *  \brief tree that holds particle data */
    BHtree<D> tree;
    
    /*! \var PlanetesimalList<D> planetesimal_list
     *  \brief TBD */
    PlanetesimalList<D> planetesimal_list;
    
};


#endif /* tree_hpp */

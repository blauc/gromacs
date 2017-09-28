#include <vector>
#include "finite3dlatticeindices.h"

namespace gmx
{

template <class T> class GridDataAccess
{
    public:
        typedef typename std::vector<T>::iterator t_iterator;
        GridDataAccess(const IVec extend, const std::vector<T> &data)
            : data_ {const_cast<std::vector<T> &>(data)},
        latticeIndices3d_ {extend} {}; // TODO: this is very dirty. Think of
                                       // proper const implementation
        GridDataAccess(IVec extend, std::vector<T> &data)
            : data_ {data}, latticeIndices3d_ {extend} {};
        /*! \brief
         * Return the raw 1d grid data
         *
         */
        std::vector<T> &data() { return data_; };
        const std::vector<T> &data() const { return data_; };
        t_iterator begin(){ return data_.begin(); };
        t_iterator end(){return data_.end(); };
        /*! \brief
         * Access a section of the volume density data via start and end iterator.
         */
        t_iterator sectionBegin(int z) const { return zy_column_begin(z, 0); }
        /*! \brief
         * Access column of the volume density data through start and end iterator.
         */
        std::array<t_iterator, 2> zy_column(int z, int y) const
        {
            return {{zy_column_begin(z, y), zy_column_begin(z, y + 1)}};
        }
        /*! \brief
         * Access the start of a column through iterator.
         * No bounds checking for fast data access
         */
        t_iterator zy_column_begin(int z, int y) const
        {
            return std::begin(data_) + latticeIndices3d_.getNumLatticePointsXY() * z +
                   latticeIndices3d_.getExtend()[XX] * y;
        };

        t_iterator next_column(t_iterator x) const
        {
            return x + latticeIndices3d_.getExtend()[XX];
        }
        t_iterator next_slice(t_iterator x) const
        {
            return x + latticeIndices3d_.getNumLatticePointsXY();
        };
        t_iterator previousSection(t_iterator x) const
        {
            return x - latticeIndices3d_.getNumLatticePointsXY();
        };

        /*! \brief
         * Directly access an index element.
         * \throws std::out_of_range if element is out of array bounds
         */
        T &at(const IVec index)
        {
            return data_.at(latticeIndices3d_.getLinearIndexFromLatticeIndex(index));
        };

        const T &get(const IVec index) const
        {
            return data_.at(latticeIndices3d_.getLinearIndexFromLatticeIndex(index));
        };

        const Finite3DLatticeIndices &indices()
        {
            return latticeIndices3d_;
        }

    private:
        std::vector<T>        &data_;
        Finite3DLatticeIndices latticeIndices3d_;
};
}

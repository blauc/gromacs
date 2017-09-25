#ifndef GMX_MATH_FINITE3DLATTICEINDICES
#define GMX_MATH_FINITE3DLATTICEINDICES

#include "gromacs/math/vectypes.h"

namespace gmx
{


typedef BasicVector<int> IVec; //!< RVec equivalent for int

class Finite3DLatticeIndices
{
    public:
        Finite3DLatticeIndices();
        Finite3DLatticeIndices(IVec extend);
        Finite3DLatticeIndices(const Finite3DLatticeIndices &other);
        int num_gridpoints() const;  //!< returns pre-evaluated extend[0]*extend[1]*extend[2]
        int numGridPointsXY() const; //!< returns pre-evaluated extend[0]*extend[1]
        IVec extend() const;         //!< return the extend of the grid
                                     /*! \brief
                                      * The extend of the grid.
                                      *
                                      * Grid indices will alwasy run from (0,0,0) to extend =
                                      * (extend[XX]-1,extend[YY]-1,extend[ZZ]-1)
                                      */
        void set_extend(IVec extend);

        /*! \brief
         * multiply grid extend by factor */
        void multiplyExtend(const RVec factor);

        /*! \brief
         * Unique one-dimensional grid index  = x + extend[XX] * y + extend[XX] *
         * extend[YY] * z.
         */
        int ndx3d_to_ndx1d(IVec i_xyz) const;

        /*! \brief
         * Inverse for ndx3d_to_ndx1d ;
         */
        IVec ndx1d_to_ndx3d(int i) const;

        bool inGrid(IVec gridIndex) const;
        bool inGrid(int gridIndex, int dimension) const;

    private:
        IVec   extend_;
        int    numGridPointsXY_;
        int    numGridPoints_;
};

} //gmx

#endif

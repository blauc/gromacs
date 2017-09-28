#ifndef GMX_MATH_FINITE3DLATTICEINDICES
#define GMX_MATH_FINITE3DLATTICEINDICES

#include "gromacs/math/vectypes.h"

namespace gmx
{

typedef BasicVector<int> IVec; //!< RVec equivalent for int

class Finite3DLatticeIndices
{
    public:
        Finite3DLatticeIndices() = default;
        Finite3DLatticeIndices(const IVec &extend);
        Finite3DLatticeIndices(const Finite3DLatticeIndices &other);
        int getNumLatticePoints() const;   //!< returns extend[0]*extend[1]*extend[2]
        int getNumLatticePointsXY() const; //!< returns pre-evaluated extend[0]*extend[1]
        IVec getExtend() const;            //!< return the extend of the lattice

        /*! \brief
         * The extend of the lattice.
         *
         * Lattice indices will alwasy run from (0,0,0) to extend =
         * (extend[XX]-1,extend[YY]-1,extend[ZZ]-1)
         */
        void setExtend(const IVec extend);

        /*! \brief
         * multiply lattice extend by factor */
        void multiplyExtend(const RVec factor);

        /*! \brief
         * Unique one-dimensional lattice index  = x + extend[XX] * y + extend[XX] *
         * extend[YY] * z.
         */
        int getLinearIndexFromLatticeIndex(const IVec &latticeIndex) const;

        /*! \brief
         * Inverse for ndx3d_to_ndx1d ;
         */
        IVec getLatticeIndexFromLinearIndex(int linearIndex) const;

        bool inLattice(const IVec &latticeIndex) const;
        bool inLattice(int latticeIndex, int dimension) const;

    private:
        IVec   extend_;
};

} //gmx

#endif

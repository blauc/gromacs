#include "externalpotential.h"
#include "externalpotentialregistration.h"

#include "gromacs/fileio/filenm.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/real.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/mdlib/groupcoord.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/fileio/oenv.h"

#include <string>
#include <memory>

ExternalPotential::ExternalPotential(){};

void ExternalPotential::dd_make_local_groups( gmx_domdec_t *dd)
{
    data_->make_local_group_indices(dd);
};

real ExternalPotential::summed_potential(t_commrec *cr)
{
    return data_->summed_potential_on_ranks(cr);
};

rvec* ExternalPotential::x_assembled(gmx_int64_t step, t_commrec *cr, rvec x[], matrix box)
{
    return data_->x_assembled(step, cr, x, box);
};

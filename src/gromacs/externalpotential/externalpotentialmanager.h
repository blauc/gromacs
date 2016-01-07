/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifndef _externalpotentialmanager_h_
#define _externalpotentialmanager_h_


#include "gromacs/externalpotential/externalpotential.h"
#include "gromacs/externalpotential/modules/template/template.h"
#include "gromacs/fileio/readinp.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/basedefinitions.h"


struct ext_pot;
struct t_blocka;
struct t_commrec;
struct t_inputrec;
struct gmx_ext_pot;
struct gmx_output_env_t;
struct gmx_mtop_t;
struct gmx_domdec_t;
struct gmx_wallcycle;


#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>


namespace externalpotential
{

struct Registry{
    struct Info{
        decltype(&Template::create) factorymethod;
        std::string description;
    };
    static const std::map<std::string, Info>
    methods;
};

namespace stringutils
{
/*! \brief
 * Helper function, that splits a string at a token
 * \todo: allow multiple tokens
 */
std::vector<std::string> split_string_at_token(const char * names, char token);

/*! \brief
 * Return vector with non-empty, non-whitespace containing string.
 *
 * E.g. "a, \t \n fs  g \n" is returned as {"a,","fs","g"}
 */
std::vector<std::string> split_string_at_whitespace(std::string s);

std::vector<std::string> split_input_at_token(int * ninp_p, t_inpfile ** inp_p, std::string mdp_line_lhs, std::string default_rhs, char token);

/*! \brief
 * An exeption throwing search_string, analogue to the one in readir.cpp */
int search_string_(const char *s, int ng, char *gn[]);

void pr_str(FILE *fp, int indent, const char * title, char * s);
};

namespace inputrecordutils
{
/*! \brief
 * Set the filenames and indexgroup names for all external potentials in grompp.
 *
 * \param[in] ninp_p mdp number of input parameters for STYPE-macro magic
 * \param[in] inp_p mdp input parameters for STYPE-macro magic
 * \param[in] ep all external potential data that goes into the input-record
 * \result groupnames for all types of external potentials
 * todo: enforce input consistency (eg. by writing checksums of input files to the .tpr file / inputrecord)
 */
char ** set_external_potential(int *ninp_p, t_inpfile **inp_p, struct ext_pot * ep );

/*! \brief
 * Generate an index group per applied external potential during pre-processing in grompp.
 *
 * \param[in] ext_pot all external potential data that goes into the input-record
 * \param[in] group_name the index group names of the external potentials
 * \param[in] groups the indexgroups as read by grompp
 * \param[in] all_group_names the group names as read by grompp
 */
void make_groups(struct ext_pot *ep, char **group_name, t_blocka *groups, char **all_group_names);

void broadcast_inputrecord_data(const t_commrec * cr, struct ext_pot * external_potential);
void pr_externalpotential(FILE *fp, int indent, t_ext_pot * pot);
};

/*! \brief
 * Manage iterations over external potentials.
 */
class Manager
{
    public:
        /*! \brief
         * Initialize the external potentials during the run (see runner.cpp).
         *
         * Check input consistency, test if files are read/writable,
         * setup coordinate communication and index groups
         *
         * \todo: Parse info from input-files in xml/json format?.
         * \todo: Check input file consistency with checksums.
         */
        Manager( FILE *fplog, t_inputrec *ir, gmx_mtop_t *mtop, rvec *x, matrix box, t_commrec *cr, const gmx_output_env_t *oenv, unsigned long Flags, gmx_bool bVerbose );

        /*! \brief
         * Trigger calculation of external potentials in the respective external potential module classes.
         * \todo: this calculation could be triggered in parallel for all external potentials
         *
         * \param[in] cr Communication record
         * \param[in] ir Input record
         * \result potential_ and force_ will be updated in all experimental input modules, if applied at this step.
         */
        void do_potential( t_commrec *cr, t_inputrec *ir, matrix box, rvec x[], real t, gmx_int64_t step, gmx_wallcycle* wcycle, gmx_bool bNS);

        /*! \brief
         * Add the forces from the external potentials to the overall force in this simulation step.
         * \todo: implement lambda-value reading for weighting the external potentials
         * \todo: implement exponential averaging for forces; make that default
         *
         * \param[in,out] f The updated forces
         * \param[in,out] vir The updated virial
         * \result contribution to the total potential from external potentials, updated force and virial
         */
        real add_forces( rvec f[], tensor vir, t_commrec *cr, gmx_int64_t step);

        /*! \brief
         * Keep the local indices in all applied potentials up-to-date after domain decomposition.
         *
         * Is called whenever the system is partioned.
         * \param[in] dd the domain decomposition data
         * \param[in,out] the external potentials
         *
         * \result updated num_atoms_loc_, ind_loc_ and nalloc_loc_,
         */
        void dd_make_local_groups( gmx_domdec_t *dd);



    private:

        std::vector<real> calculate_weights();

        int nr_colons_in_string_(std::string s);

        /*! \brief
         * Throw errors if input is inconsistent
         */
        void throw_at_input_inconsistency_(t_commrec * cr, t_inputrec * ir, std::string input_file, std::string output_file, int current);

        std::vector<std::unique_ptr<ExternalPotential> > potentials_;
        std::vector<real> V_external_;

};
} // namespace externalpotential
#endif

#ifndef _externalpotentialmanager_h_
#define _externalpotentialmanager_h_


#include "gromacs/utility/basedefinitions.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/gmxlib/readinp.h"
#include "gromacs/mdtypes/inputrec.h"
#include "externalpotential.h"
#include "modules/densityfitting/densityfitting.h"
#include "modules/pulling/pulling.h"



struct ext_pot;
struct t_blocka;
struct t_commrec;
struct t_inputrec;
struct gmx_ext_pot;
struct gmx_output_env_t;
struct gmx_mtop_t;
struct gmx_domdec_t;
struct gmx_wallcycle;


#include <map>
#include <vector>
#include <string>
#include <memory>
#include <functional>


namespace externalpotential{
    typedef std::function<std::unique_ptr<ExternalPotential>(struct ext_pot_ir, t_commrec *, t_inputrec *, const gmx_mtop_t*, rvec [], matrix, FILE *, FILE *, FILE *, gmx_bool , const gmx_output_env_t *, unsigned long)> creatorFunction ;


/*! \brief
 * Manage iterations over external potentials.
 */
class Manager
{
    public:

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

        /*! \brief
         * Trigger calculation of external potentials in the respective external potential module classes.
         * \todo: this calculation could be triggered in parallel for all external potentials
         *
         * \param[in] cr Communication record
         * \todo: create own communicators for all external potentials
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
        real add_forces( rvec f[], tensor vir, t_commrec *cr, gmx_int64_t step, real t);

        /*! \brief
         * Initialize the external potentials during the run (see runner.cpp).
         *
         * Check input consistency, test if files are read/writable,
         * setup coordinate communication and index groups
         *
         * \todo: Parse info from input-files in xml/json format?.
         * \todo: Check input file consistency with checksums.
         */
        void init_external_potentials( FILE *fplog, t_inputrec *ir, gmx_mtop_t *mtop, rvec *x, matrix box, t_commrec *cr, const gmx_output_env_t *oenv, unsigned long Flags, gmx_bool bVerbose );

        /*! \brief
         * Keep the local indices in all applied potentials up-to-date after domain decomposition.
         *
         * Is called whenever the system is partioned.
         * \param[in] dd the domain decomposition data
         * \param[in,out] the external potentials
         *
         * \result updated num_atoms_loc_, ind_loc_ and nalloc_loc_,
         */
        void dd_make_local_groups( gmx_domdec_t *dd, struct ext_pot *external_potential );

        void broadcast_inputrecord_data(const t_commrec * cr, struct ext_pot * external_potential);

    private:

        /*! \brief
         * Checks if the number of inputfile arguments matches the number of outputilfe matches groupnames
         * (Arguments can be empty, but should be present) */
        bool matching_number_of_inputs_(std::string filenames, std::string output, std::string groupnames);
        int nr_colons_in_string_(std::string s);

        /*! \brief
         * Throw errors if input is inconsistent
         */
        void throw_at_input_inconsistency_(t_commrec * cr, t_inputrec * ir, std::string input_file, std::string output_file, int current);

        std::vector<std::string> split_input_at_token_(int * ninp_p, t_inpfile ** inp_p, std::string mdp_line_lhs, std::string default_rhs,char token);

        /*! \brief
         * Helper function, that splits a string at a token
         * \todo: allow multiple tokens
         */
        std::vector<std::string> split_string_at_token_(const char * names, char token);

        /*! \brief
         * An exeption throwing search_string, analogue to the one in readir.cpp */
        int search_string_(const char *s, int ng, char *gn[]);

        /*! \brief
         * Return vector with non-empty, non-whitespace containing string.
         *
         * E.g. "a, \t \n fs  g \n" is returned as {"a,","fs","g"}
         */
        std::vector<std::string> split_string_at_whitespace_(std::string s);
        std::vector<std::pair<std::string,decltype(&DensityFitting::create)>> methods_={{"density-fitting",&DensityFitting::create}};

        std::vector<std::unique_ptr<ExternalPotential> > potentials_;

};
} // namespace externalpotential
#endif

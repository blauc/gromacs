#include "externalpotentialutil.h"




// void ExternalPotentialUtil::set_external_potential(t_ext_pot external_potential,
//   char **external_potential_filenames, int n_external_potentials)
// {
//     snew(external_potential,1);
//     external_potential->number_of_potentials = n_external_potentials;
//
//     /* set the filenames */
//
//     snew(external_potential->filenames,n_external_potentials);
//     for (int i = 0; i < n_external_potentials; i++) {
//         snew(external_potential->filenames[i], STRLEN);
//         external_potential->filenames[i] = external_potential_filenames[i];
//     }
//
//     return;
// };
//

// void ExternalPotentialUtil::do_external_potentials(t_commrec *cr, t_inputrec *ir, matrix box, rvec x[],
    // real t, gmx_int64_t step, gmx_wallcycle_t wcycle, gmx_bool  bNS)
// {
//     for (std::list<ExternalPotential*>::iterator it =
//         ir->external_potential->extpot->potentials.begin();
//         it != ir->external_potential->extpot->potentials.end(); it++)
//     {
//         (*it)->do_potential(cr, ir, box, x, t, step, wcycle, bNS);
//     }
//     return;
// };

// real ExternalPotentialUtil::add_ext_forces(t_gmx_ext_pot *extpot, rvec f[], t_commrec *cr, gmx_int64_t step, real t){
//
//     real result=0;
//     for (std::list<ExternalPotential*>::iterator it =
//         extpot->potentials.begin();
//         it != extpot->potentials.end(); it++)
//     {
//          result+= (*it)->potential();
//         (*it)->add_forces(f, cr, step, t);
//     }
//     return result;
// };


// void ExternalPotentialUtil::init_external_potentials(FILE *fplog, t_inputrec *ir, int nfile,
//         const t_filenm fnm[], gmx_mtop_t *mtop, rvec *x, matrix box,
//         t_commrec *cr, const output_env_t oenv, unsigned long Flags,
//         gmx_bool bVerbose)
// {
//
//     FILE* current_potential_file;
//     t_external_potential_type current_type;
//     char * current_type_string;
//     ExternalPotential* external_potential_buf;
//
//     for (size_t i = 0; i < external_potential->number_of_potentials; i++)
//     {
//
//         current_potential_file = gmx_ffopen(external_potential->filename[i],"r");
//         if (current_potential_file != NULL) {
//             if (fread(&current_type_string, sizeof(current_type_string), 8, current_potential_file)!=8)
//             {
//                 gmx_error("input","Cannot read external potential format string.");
//             }
//         }
//
//         external_potential_buf = external_potential_string_to_init(current_type_string,
//           fplog, inputrec, nfile, fnm, mtop, state->x, box, cr, oenv, Flags, bVerbose);
//         if(external_potential_buf == NULL)
//         {
//             gmx_error("input","Failed to initialise external potential.");
//         }
//         else
//         {
//             ir->external_potential->external_potentals.push_back(external_potential_buf);
//         }
//
//         gmx_ffclose(current_potential_file);
//     }
//
//     return;
// };

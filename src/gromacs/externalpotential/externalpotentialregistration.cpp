//
// void ExternalPotentialRegistration::init(
//     const char *typestring,
//     FILE * inputfile,
//     FILE *fplog,
//     t_inputrec *ir,
//     int nfile,
//     const t_filenm fnm[],
//     gmx_mtop_t *mtop,
//     rvec *x,
//     matrix box,
//     t_commrec *cr,
//     const output_env_t oenv,
//     unsigned long Flags,
    // gmx_bool bVerbose){
    //
    //     char error_msg[STRLEN];
    //     this->registered_potential_ = NULL;
    //     this->success_ = false;
    //     this->known_method_=false;
    //
    //     if (gmx_strncasecmp_min("CRYOEM3D",typestring,STRLENTYPE))
    //     {
    //         init_instance<DensityFitting>(inputfile, fplog, ir,
    //             nfile, fnm, mtop, x, box,cr, oenv, Flags, bVerbose);
    //     };
    //
    //     if (gmx_strncasecmp_min("PULLING ",typestring,STRLENTYPE))
    //     {
    //         init_instance<Pulling>(inputfile, fplog, ir,
    //             nfile, fnm, mtop, x, box,cr, oenv, Flags, bVerbose);
    //     };
    //
    //     if(!this->known_method_)
    //     {
    //         sprintf(error_msg,"Cannot recognise ""%s"" method in first eight characters of the given input file.", typestring);
    //
    //     }
    //
    //     if(!this->success_)
    //     {
    //         sprintf(error_msg,"Recognised ""%s"", but failed to initialise external potential.", typestring);
    //         gmx_error("input",error_msg);
    //     }
    //
    //
    // };

// ExternalPotential* ExternalPotentialRegistration::registered_potential(){
//     ExternalPotential* result;
//     result=registered_potential_;
//     this->success_ = false;
//     this->known_method_ = false;
//     this->registered_potential_ = NULL;
//     return result;
// };

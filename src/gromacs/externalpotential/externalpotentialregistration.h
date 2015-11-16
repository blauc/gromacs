/*! \file
 * \brief
 * If you provided a new external potential module, this is the place to let gromacs know
 */
#ifndef _externalpotentialregistration_h_
#define _externalpotentialregistration_h_

class ExternalPotential;
class ExternalPotentialData;
class ExternalPotentialInfo;

#include <string>
#include <vector>
#include <memory>

typedef std::shared_ptr<ExternalPotentialData> ExternalPotentialDataPointer;

/*! \brief
 * This class registers the external potentials by generating a pointer to each impemented external potential.
 *
 * The external potentials are invoked via the create-Routine of the ModuleInfo-classes.
 */
class ExternalPotentialRegistration {
    public:
        ExternalPotentialRegistration();
        ExternalPotential* init(size_t method, ExternalPotentialDataPointer data);
        size_t number_methods();
        size_t method_id(std::string methodstring);
        std::string name(size_t method_nr);
        std::vector<std::string> names();
    private:
        std::vector<ExternalPotentialInfo*> methods_;
        template <class T> void register_(){
            methods_.push_back( (ExternalPotentialInfo*) (new T) );
        };
};

#endif

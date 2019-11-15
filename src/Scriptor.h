//
// Created by parnet on 08.06.19.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_OUTPUT_H
#define UG_PLUGIN_XBRAIDFORUG4_OUTPUT_H

#include "lib_disc/io/vtkoutput.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/interpolate.h"
#include "lib_disc/spatial_disc/domain_disc.h"
#include "lib_algebra/vector_interface/vec_functions.h"

template<typename TDomain, typename TAlgebra>
class Scriptor {
public:
    typedef ug::GridFunction <TDomain, TAlgebra> TGridFunction;
    typedef SmartPtr <TGridFunction> SPGridFunction;

    Scriptor() = default;

    virtual ~Scriptor() {} ;

    virtual bool write(TGridFunction *u, int index, double time) = 0;

    virtual bool write(TGridFunction *u, int index, double time, int iteration, int level) {};

    virtual void print(const char *filename, TGridFunction &u, int index, double time) {};
};

/**
 * Compatibility class to simulate ug::VTKOutput
 * @tparam TDomain
 * @tparam TAlgebra
 */
template<typename TDomain, typename TAlgebra>
class VTKScriptor : public Scriptor<TDomain, TAlgebra> {
public:
    typedef ug::GridFunction <TDomain, TAlgebra> TGridFunction;
    typedef SmartPtr<TGridFunction> SPGridFunction;
    typedef SmartPtr<ug::VTKOutput<TDomain::dim> > SPVTKOutput;
    typedef Scriptor<TDomain, TAlgebra> TScriptor;

    SPVTKOutput m_out;
    const char *m_filename;

    VTKScriptor(SPVTKOutput p_out, const char *filename) : Scriptor<TDomain, TAlgebra>() {
        m_out = p_out;
        m_filename = filename;
    }

    ~VTKScriptor() {}

    bool write(TGridFunction *u, int index, double time) override {
        this->m_out->print(m_filename, *u, index, time);
        return true;
    };

    bool write(TGridFunction *u, int index, double time, int iteration, int level) override {
        std::stringstream ss;
        ss << m_filename << "_k" << iteration << "_l" << level;
        this->m_out->print(ss.str().c_str(), *u, index, time);
        return true;
    };

    void writeTimePVD(TGridFunction *u) {
        this->m_out->write_time_pvd(this->m_filename, *u);
    };

    void write_time_pvd(const char *filename, TGridFunction &u) {
        this->writeTimePVD(&u);
    };
};

template<typename TDomain, typename TAlgebra>
class VTKIterationScriptor : public Scriptor<TDomain, TAlgebra> {
public:
    typedef ug::GridFunction <TDomain, TAlgebra> TGridFunction;
    typedef SmartPtr<TGridFunction> SPGridFunction;
    typedef SmartPtr<ug::VTKOutput<TDomain::dim> > SPVTKOutput;

    SPVTKOutput m_out;
    const char *m_filename;

    VTKIterationScriptor(SPVTKOutput p_out, const char *filename) : Scriptor<TDomain, TAlgebra>() {
        m_out = p_out;
        m_filename = filename;
    }

    ~VTKIterationScriptor()= default;

    bool write(TGridFunction *u, int index, double time) override {
        this->m_out->print(m_filename, *u, index, time);
        return true;
    };

    bool write(TGridFunction *u, int index, double time, int iteration, int level) override {
        std::stringstream ss;
        ss << m_filename << "_l" << level << "_k" << index ;
        this->m_out->print(ss.str().c_str(), *u, iteration, time);
        return true;
    };

    void writeTimePVD(TGridFunction *u) {
        this->m_out->write_time_pvd(this->m_filename, *u);
    };

    void write_time_pvd(const char *filename, TGridFunction &u) {
        this->writeTimePVD(&u);
    };
};

template<typename TDomain, typename TAlgebra>
class VTKModScriptor : public Scriptor<TDomain, TAlgebra> {
public:
    typedef ug::GridFunction <TDomain, TAlgebra> TGridFunction;
    typedef SmartPtr<TGridFunction> SPGridFunction;
    typedef SmartPtr<ug::VTKOutput<TDomain::dim> > SPVTKOutput;

    SPVTKOutput m_out;
    const char *m_filename;
    int modal = 1;

    VTKModScriptor(SPVTKOutput p_out, const char *filename) : Scriptor<TDomain, TAlgebra>() {
        m_out = p_out;
        m_filename = filename;
    }

    ~VTKModScriptor() {}

    void writeTimePVD(TGridFunction *u) {
        this->m_out->write_time_pvd(this->m_filename, *u);
    }

    void write_time_pvd(const char *filename, TGridFunction &u) {
        this->writeTimePVD(&u);
    }

    void setModal(int p_modal) {
        this->modal = p_modal;
    }

    bool write(TGridFunction *u, int index, double time) override {
        if (index % modal == 0) {
            this->m_out->print(m_filename, *u, index, time);
        }
        return true;
    };

    virtual bool write(TGridFunction *u, int index, double time, int iteration, int level) override {
        if (index % this->modal == 0) {
            std::stringstream ss;
            ss << m_filename << "_k" << iteration << "_l" << level;
            this->m_out->print(ss.str().c_str(), *u, index, time);
        }
        return true;
    };

    void print(const char *filename, TGridFunction &u, int index, double time) override {
        TGridFunction *uptr = &u;
        this->write(uptr, index, time);
    }
};

template<typename TDomain, typename TAlgebra>
class EvalScriptor : public Scriptor<TDomain, TAlgebra> {
public:
    typedef ug::GridFunction <TDomain, TAlgebra> TGridFunction;
    typedef SmartPtr<TGridFunction> SPGridFunction;
    typedef SmartPtr<ug::VTKOutput<TDomain::dim> > SPVTKOutput;
    typedef SmartPtr<ug::IDomainDiscretization<TAlgebra> > SPDomainDisc;

    typedef SmartPtr<ug::UserData<double, TGridFunction::dim> > SPData;
    SPData m_data;

    const char *m_cmp;

    bool m_relative = false;
    std::ofstream *outfile;
    SPDomainDisc m_domainDisc;

    EvalScriptor() : Scriptor<TDomain, TAlgebra>() {}

    ~EvalScriptor() {
        outfile->close();
        delete outfile;
    }

    void setDomain(SPDomainDisc p_domainDisc) {
        m_domainDisc = p_domainDisc;
    }

    void setRelative(bool val){
        this->m_relative = val;
    }

    void setFile(const char *filename) {
        outfile = new std::ofstream();
        outfile->open(filename);
        (*this->outfile) << "iteration; level; index; time; || u - v ||;" ;
        if(this->m_relative){
            (*this->outfile) << "  || u ||; relative;";
        }
        (*this->outfile) << std::endl;
    }

    void setFile(std::ofstream &file) {
        this->outfile = &file;
    }

    void writeTimePVD(TGridFunction *u) {};

    bool write(TGridFunction *u, int index, double time) override {
	    return this->write(u,index,time,-1,-1);
    };

    virtual bool write(TGridFunction *u, int index, double time, int iteration, int level) override  {
        SPGridFunction vec = u->clone_without_values();

        Interpolate(this->m_data, vec, this->m_cmp, NULL, time); // todo eliminate NULL ?
        m_domainDisc->adjust_solution(*vec.get(), time);

        VecAdd(-1.0, *vec.get(), 1.0, *u);
	if(iteration == -1){
		(*this->outfile) << std::setw(4) << ""<<";"<< std::setw(4) << "" <<";";
	} else {
		(*this->outfile) << std::setw(4) << iteration << ";"<< std::setw(4)  << level << ";";
	}

	    double vecnorm = vec->norm();
        (*this->outfile) << std::setw(7) << index <<";"
                         << std::setw(12) << time << ";"
                         << std::setw(12) << vecnorm << ";";
	    if(this->m_relative){
            SPGridFunction uvec = u->clone();
            double unorm = uvec->norm();
            (*this->outfile) << std::setw(12) << unorm << ";"
                             << std::setw(12) << (vecnorm / unorm) << ";";
        }

            (*this->outfile) << std::endl;

        return true;
    };

    void print(const char *filename, TGridFunction &u, int index, double time) override {
        TGridFunction *uptr = &u;
        this->write(uptr, index, time);
    }

    void write_time_pvd(const char *filename, TGridFunction &u) {
    }


    void setGeneratorComponent(const char *cmp) {
        this->m_cmp = cmp;
    }


    void setVectorGenerator(SPData p_data) {
        this->m_data = p_data;
    }

    void setDomainDisc(SPDomainDisc p_domainDisc) {
        m_domainDisc = p_domainDisc;
    }

#ifdef UG_FOR_LUA

    void setVectorGenerator(const char *fctName) {
        setVectorGenerator(ug::LuaUserDataFactory<double, TDomain::dim>::create(fctName));
    }

    void setVectorGenerator(ug::LuaFunctionHandle fct) {
        setVectorGenerator(make_sp(new ug::LuaUserData<double, TDomain::dim>(fct)));
    }

#endif
};

template<typename TDomain, typename TAlgebra>
class MATLABScriptor : public Scriptor<TDomain, TAlgebra> {
public:
    typedef ug::GridFunction <TDomain, TAlgebra> TGridFunction;
    typedef SmartPtr<TGridFunction> SPGridFunction;
    typedef SmartPtr<ug::VTKOutput<TDomain::dim> > SPVTKOutput;

    SPVTKOutput m_out;
    std::ofstream &debugwriter;

    MATLABScriptor(std::ofstream &p_debugwriter) : Scriptor<TDomain, TAlgebra>(), debugwriter(p_debugwriter) {
    }

    ~MATLABScriptor() {

    }

    void setFile(std::ofstream &file) {
        this->debugwriter = file;
    }


    bool write(TGridFunction *u, int index, double time) override {

        int sz = u->size();
        this->debugwriter << "u_" << index << " = [";
        for (int i = 0; i < sz - 1; i++) {
            this->debugwriter << std::setw(7) << (*u)[i] << ", ";
        }
        this->debugwriter << std::setw(7) << (*u)[sz - 1];
        this->debugwriter << "]" << std::endl;

        return true;
    };


    virtual bool write(TGridFunction *u, int index, double time, int iteration, int level) override  {
        this->write(u, index, time);
        return true;
    };


    void print(const char *filename, TGridFunction &u, int index, double time) override  {
        TGridFunction *uptr = &u;
        this->write(uptr, index, time);
    }

};

template<typename TDomain, typename TAlgebra>
class MultiScriptor : public Scriptor<TDomain, TAlgebra> {

public:
    typedef typename ug::GridFunction<TDomain, TAlgebra> TGridFunction;
    typedef SmartPtr<TGridFunction> SPGridFunction;
    typedef SmartPtr<ug::VTKOutput<TDomain::dim> > SPVTKOutput;
    typedef Scriptor<TDomain, TAlgebra> TScriptor;



    MultiScriptor(){}

    MultiScriptor(const char *filename) : Scriptor<TDomain, TAlgebra>() {
        //lst = std::vector<SmartPtr<Scriptor<TDomain, TAlgebra>>>();
    }

    ~MultiScriptor() {

    }

    void addScriptor(SmartPtr<TScriptor> scriptor) {
        // this->lst.emplace_back(scriptor);
    	this->lst.push_back(scriptor);
    }

    void writeTimePVD(TGridFunction *u) {};

    void write_time_pvd(const char *filename, TGridFunction &u) {};

    bool write(TGridFunction *u, int index, double time) override {
        for (size_t k = 0; k < this->lst.size(); k++) {
            this->lst[k]->write(u, index, time);
        }
        return true;
    };


    virtual bool write(TGridFunction *u, int index, double time, int iteration, int level) override {
        for (size_t k = 0; k < this->lst.size(); k++) {
            this->lst[k]->write(u, index, time, iteration, level);
        }
        return true;
    };

    void print(const char *filename, TGridFunction &u, int index, double time) override {
        for (size_t k = 0; k < this->lst.size(); k++) {
            this->lst[k]->print(filename, u, index, time);
        }
    }
protected:
    std::vector<SmartPtr<TScriptor> > lst;
};

std::string ptr2str(void *c) {
    std::stringstream ss;
    ss << c;
    return ss.str();
}

#endif //UG_PLUGIN_XBRAIDFORUG4_OUTPUT_H

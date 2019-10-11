//
// Created by mparnet on 19.07.19.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_FIELDGENERATOR_H
#define UG_PLUGIN_XBRAIDFORUG4_FIELDGENERATOR_H

#include "lib_disc/spatial_disc/user_data/std_glob_pos_data.h"
#include "common/math/math_vector_matrix/math_vector.h"

template<typename TDatatype, int dim>
class SinSourceOneCube : public ug::StdGlobPosData<SinSourceOneCube<TDatatype, dim>, TDatatype, dim> {
public:
    SinSourceOneCube() = default;

    ~SinSourceOneCube() = default;

    std::string config_string() const {
        return "-sin(pi*x)*sin(pi*y)*sin(pi*z)*(sin(t)-dim*alpha*pi*pi*cos(t))";
    };

    TDatatype alpha;

    void setAlpha(TDatatype alpha) {
        this->alpha = alpha;
    };

    inline void evaluate(TDatatype &value, const ug::MathVector<dim> &x, number time, int si) const {
        TDatatype prod = 1;
        for (size_t i = 0; i < dim; i++) {
            prod *= sin(M_PI * x[i]);
        }
        value = -1 * prod * (sin(time) - dim * alpha * M_PI * M_PI * cos(time));
    }
};


template<typename TDatatype, int dim>
class SinAnalyticSolutionOneCube
        : public ug::StdGlobPosData<SinAnalyticSolutionOneCube<TDatatype, dim>, TDatatype, dim> {
public:
    SinAnalyticSolutionOneCube() = default;

    ~SinAnalyticSolutionOneCube() = default;

    std::string config_string() const {
        return "sin(pi*x)*sin(pi*y)*sin(pi*z)*cos(t)";
    };

    TDatatype alpha;

    void setAlpha(TDatatype alpha) {
        this->alpha = alpha;
    };

    inline void evaluate(TDatatype &value, const ug::MathVector<dim> &x, number time, int si) const {
        TDatatype prod = 1;
        for (size_t i = 0; i < dim; i++) {
            prod *= sin(M_PI * x[i]);
        }
        value = prod * cos(time);
    }
};


template<typename TDatatype, int dim>
class SinSourcePiCube : public ug::StdGlobPosData<SinSourcePiCube<TDatatype, dim>, TDatatype, dim> {
public:
    SinSourcePiCube() = default;

    ~SinSourcePiCube() = default;

    std::string config_string() const {
        return "-sin(x)*sin(y)*sin(z)*(sin(t)-dim*alpha*cos(t))";
    };

    TDatatype alpha;

    void setAlpha(TDatatype alpha) {
        this->alpha = alpha;
    };

    inline void evaluate(TDatatype &value, const ug::MathVector<dim> &x, number time, int si) const {
        TDatatype prod = 1;
        for (size_t i = 0; i < dim; i++) {
            prod *= sin(x[i]);
        }
        value = -1 * prod * (sin(time) - dim * alpha * cos(time));
    }
};

template<typename TDatatype, int dim>
class SinAnalyticSolutionPiCube : public ug::StdGlobPosData<SinAnalyticSolutionPiCube<TDatatype, dim>, TDatatype, dim> {
public:
    SinAnalyticSolutionPiCube() = default;

    ~SinAnalyticSolutionPiCube() = default;

    std::string config_string() const {
        return "sin(x)*sin(y)*sin(z)*cos(t)";
    };

    TDatatype alpha;

    void setAlpha(TDatatype alpha) {
        this->alpha = alpha;
    };

    inline void evaluate(TDatatype &value, const ug::MathVector<dim> &x, number time, int si) const {
        TDatatype prod = 1;
        for (size_t i = 0; i < dim; i++) {
            prod *= sin(x[i]);
        }
        value = prod * cos(time);
    }
};


#endif //UG_PLUGIN_XBRAIDFORUG4_FIELDGENERATOR_H

/**
 * @file    m4dGeodesicDP54.h
 * @author  Thomas Mueller
 *
 *   Calculate geodesics with a Dormand-Prince [RK6(5)] method.
 *
 *        The algorithm is taken from  Andreas Guthmann,
 *        "Einfuehrung in die Himmelsmechanik und Ephemeridenrechnung",
 *        Spektrum Verlag (2000), 2. Auflage, Seite 276
 *
 *        Order of the method is p=5;
 *
 *  This file is part of libMotion4D.
 */
#ifndef M4D_GEODESIC_DP54_H
#define M4D_GEODESIC_DP54_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include <m4dGlobalDefs.h>
#include <metric/m4dMetric.h>
#include <motion/m4dGeodesic.h>

#ifdef USE_DP_INT

namespace m4d {

// ---------------------------------------------------
//    class definition:   GeodesicDP54
// ---------------------------------------------------
class API_M4D_EXPORT GeodesicDP54 : public Geodesic
{
public:
    GeodesicDP54(Metric* metric, enum_geodesic_type type = enum_geodesic_lightlike);
    virtual ~GeodesicDP54();

    // --------- public methods -----------
public:
    void setStepSizeControlled(bool control = true);

    virtual enum_break_condition calculateGeodesic(const vec4 initPos, const vec4 initDir, const int maxNumPoints,
        std::vector<vec4>& points, std::vector<vec4>& dirs, std::vector<double>& lambda);

    virtual enum_break_condition calculateGeodesic(
        const vec4 initPos, const vec4 initDir, const int maxNumPoints, vec4*& points, vec4*& dirs, int& numPoints);

    virtual enum_break_condition calculateGeodesicData(const vec4 initPos, const vec4 initDir, const int maxNumPoints,
        std::vector<vec4>& points, std::vector<vec4>& dirs, std::vector<double>& epsilons, std::vector<double>& lambda);

    virtual enum_break_condition calcParTransport(const vec4 initPos, const vec4 initDir, const vec4 e0, const vec4 e1,
        const vec4 e2, const vec4 e3, const int maxNumPoints, std::vector<vec4>& points, std::vector<vec4>& dirs,
        std::vector<double>& lambda, std::vector<vec4>& base0, std::vector<vec4>& base1, std::vector<vec4>& base2,
        std::vector<vec4>& base3);

    virtual enum_break_condition calcSachsJacobi(const vec4 initPos, const vec4 initCoordDir, const vec3 localNullDir,
        const vec3 locX, const vec3 locY, const vec3 locZ, const vec4 e0, const vec4 e1, const vec4 e2, const vec4 e3,
        const enum_nat_tetrad_type tetrad_type, const int maxNumPoints, std::vector<vec4>& points,
        std::vector<vec4>& dirs, std::vector<double>& lambda, std::vector<vec4>& sachs0, std::vector<vec4>& sachs1,
        std::vector<vec5>& jacobi, vec5& maxJacobi);

    virtual enum_break_condition calcSachsJacobi(const vec4 initPos, const vec4 initCoordDir, const vec3 localNullDir,
        const vec3 locX, const vec3 locY, const vec3 locZ, const vec4 b0, const vec4 b1, const vec4 b2, const vec4 b3,
        const enum_nat_tetrad_type tetrad_type, const int maxNumPoints, vec4*& points, vec4*& dirs, double*& lambda,
        vec4*& sachs0, vec4*& sachs1, vec5*& jacobi, vec5& maxJacobi, int& numPoints);

    // void  nextStep              ( double* yo, double* yn, double* yerr, double h );
    // void  nextStepPar           ( double* yo, double* yn, double* yerr, double h );
    // void  nextStepSachsJacobi   ( double* yo, double* yn, double* yerr, double h );

    virtual bool nextStep(int& status);
    virtual bool nextStepPar(int& status);
    virtual bool nextStepSachsJacobi(int& status);

    virtual void print(FILE* fptr = stderr);

    // -------- protected attribute ---------
protected:
    bool mCalcWithParTransport;
    int mNumCoords;

    double a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, a71, a72, a73, a74, a75, a76;
    double b1, b2, b3, b4, b5, b6, b7, db1, db2, db3, db4, db5, db6, db7;

    int mOrder;
    double stepSigma;
    double stepFac;

    double yn[DEF_MAX_YS], yerr[DEF_MAX_YS], h;
};

} // end namespace m4d

#endif // USE_DP_INT

#endif

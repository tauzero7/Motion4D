// -------------------------------------------------------------------------------
/*
    m4dObject.cpp

  Copyright (c) 2009-2014  Thomas Mueller, Frank Grave


   This file is part of the m4d-library.

   The m4d-library is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The m4d-library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with the m4d-library.  If not, see <http://www.gnu.org/licenses/>.

*/
// -------------------------------------------------------------------------------

#include "m4dObject.h"

namespace m4d {

/**
 * @brief Standard constructor.
 */
Object::Object() :
    currMetric(NULL),
    geodSolver(NULL) 
{
    resetAll();
}

/**
 * @brief Standard destructor.
 */
Object::~Object() {
    clearAll();
}


bool Object::setMetric(const char *metricName) {
    if (currMetric != nullptr) {
        delete currMetric;
    }
    currMetric = metricDB.getMetric(metricName);
    return (currMetric == nullptr ? false : true);
}


bool Object::setMetricParam(const char *paramName, double value) {
    if (currMetric == nullptr) {
        fprintf(stderr,"Object::setMetricParam() ... metric is missing!\n");
        return false;
    }
    return currMetric->setParam(paramName, value);
}


bool Object::setSolver(const char *solverName) {
    if (currMetric == nullptr) {
        fprintf(stderr,"Object::setSolver() ... metric is missing!\n");
        return false;
    }
    if (geodSolver != nullptr) {
        delete geodSolver;
    }
    geodSolver = solverDB.getIntegrator(currMetric, solverName);
    return (geodSolver == nullptr ? false : true);
}


bool Object::setSolverParam(const char *paramName, bool val) {
    if (geodSolver == nullptr) {
        fprintf(stderr,"Object::setSolverParam() ... solver is missing!\n");
        return false;
    }
    //M4D_MSG((val ? "yes" : "no"));
    return geodSolver->setParam(paramName, val);
}


bool Object::setSolverParam(const char *paramName, double value) {
    if (strcmp(paramName, "eps_a") == 0) {
        this->epsAbs = value;
    }
    else if (strcmp(paramName, "eps_r") == 0) {
        this->epsRel = value;
    }
    else if (strcmp(paramName, "stepsize") == 0) {
        this->stepsize = value;
    }

    if (geodSolver == nullptr) {
        fprintf(stderr,"Object::setSolverParam() ... solver is missing!\n");
        return false;
    }    
    //M4D_MSG(std::to_string(value).c_str());
    return geodSolver->setParam(paramName, value);
}


bool Object::setSolverParam(const char *paramName, double v0, double v1, double v2, double v3) {
    if (geodSolver == nullptr) {
        fprintf(stderr,"Object::setSolverParam() ... solver is missing!\n");
        return false;
    }
    return geodSolver->setParam(paramName, v0, v1, v2, v3);
}

bool Object::setParam(const char* paramName, int val) {
    if (paramName != nullptr && strcmp(paramName, "maxNumPoints") == 0) {
        maxNumPoints = static_cast<unsigned int>(val);
    }
    return true;
}


bool Object::setInitialPosition(const double* x) {
    if (x == nullptr) {
        return false;
    }
    return setInitialPosition(x[0], x[1], x[2], x[3]);
}


bool Object::setInitialPosition(double x0, double x1, double x2, double x3) {
    if (currMetric == nullptr) {
        fprintf(stderr,"Object::setInitialPosition() ... metric is missing!\n");
        return false;
    }
    this->startPos = vec4(x0, x1, x2, x3);

    // test if initial position is valid
    currMetric->calculateMetric(this->startPos);
    return (currMetric->breakCondition(this->startPos) == enum_break_none ? true : false);
}


bool Object::setInitialDirection(const double *v) {
    if (v == nullptr) {
        return false;
    }
    return setInitialDirection(v[0], v[1], v[2], v[3]);
}


bool Object::setInitialDirection(double v0, double v1, double v2, double v3) {
    this->coordDir = vec4(v0, v1, v2, v3);
    return true;
}


bool Object::setInitialLocalNullDirection(enum_time_direction tdir, const double *l,
                                        enum_nat_tetrad_type nattype)
{
    if (l == nullptr) {
        return false;
    }
    return setInitialLocalNullDirection(tdir, l[0], l[1], l[2], nattype);
}

bool Object::setInitialLocalNullDirection(enum_time_direction tdir,
                                          double l0, double l1, double l2,
                                          enum_nat_tetrad_type natType)
{
    if (currMetric == nullptr) {
        fprintf(stderr,"Object::setInitialPosition() ... metric is missing!\n");
    }

    double tf = (tdir == enum_time_backward ? -1.0 : 1.0);

    vec4 locDir, coDir;
    vec3 dir = vec3(l0, l1, l2).getNormalized();
    locDir = vec4(tf, dir[0], dir[1], dir[2]);
    this->startDir = dir;
    currMetric->localToCoord(this->startPos, locDir, coDir, natType);
    this->coordDir = coDir;
    this->type = enum_geodesic_lightlike;
    return true;
}


bool Object::setInitialLocalTimeDirection(enum_time_direction tdir,
                                    double l0, double l1, double l2, double beta,
                                    enum_nat_tetrad_type natType)
{
    if (currMetric == nullptr) {
        fprintf(stderr,"Object::setInitialPosition() ... metric is missing!\n");
    }

    double tf = (tdir == enum_time_backward ? -1.0 : 1.0);
    if (fabs(beta) >= 1.0) {
        fprintf(stderr,"Object::setInitialPosition() ... velocity too fast!\n");
        return false;
    }
    double gam = 1.0/sqrt(1.0 - beta*beta);

    vec4 locDir, coDir;
    vec3 dir = vec3(l0, l1, l2).getNormalized();

    locDir = vec4(tf*gam, gam*beta*dir[0], gam*beta*dir[1], gam*beta*dir[2]);
    currMetric->localToCoord(this->startPos, locDir, coDir, natType);
    this->coordDir = coDir;
    this->type = enum_geodesic_timelike;
    return true;
}


enum_break_condition Object::calculateGeodesic() {
    if (geodSolver == nullptr) {
        fprintf(stderr,"Object::calculateGeodesic() ... solver is missing!\n");
        return enum_break_other;
    }
    clearAll();
    geodSolver->setGeodesicType(this->type);   
    geodSolver->setAffineParamStep(this->stepsize);
    return geodSolver->calculateGeodesic(this->startPos, this->coordDir, this->maxNumPoints,
                                         this->points, this->dirs, this->lambda);
}


enum_break_condition  Object::calcSachsJacobi() {
    if (geodSolver == nullptr) {
        fprintf(stderr,"Object::calcSachsJacobi() ... solver is missing!\n");
        return enum_break_other;
    }
    clearAll();
    geodSolver->setGeodesicType(enum_geodesic_lightlike_sachs);

    vec3 localNullDir = this->startDir;

    vec3 locX = vec3(1.0, 0.0, 0.0);
    vec3 locY = vec3(0.0, 1.0, 0.0);
    vec3 locZ = vec3(0.0, 0.0, 1.0);

    return geodSolver->calcSachsJacobi(this->startPos, this->coordDir, 
        localNullDir, locX, locY, locZ, 
        this->base[0], this->base[1], this->base[2], this->base[3],
        this->tetradType, this->maxNumPoints,
        this->points, this->dirs, this->lambda, this->sachs1, this->sachs2,
        this->jacobi, this->maxJacobi);
}


void Object::printStatus() {
    if (currMetric != nullptr) {
        currMetric->printF();
    }
    if (geodSolver != nullptr) {
        geodSolver->printF();
    }

    fprintf(stderr,"\n");
    fprintf(stderr,"initial position:  %f %f %f %f\n",
            this->startPos[0], this->startPos[1], this->startPos[2], this->startPos[3]);
    fprintf(stderr,"initial direction: %f %f %f %f\n",
            this->coordDir[0], this->coordDir[1], this->coordDir[2], this->coordDir[3]);

    fprintf(stderr,"\n");
    fprintf(stderr,"#points: %d\n", (int)this->points.size());
}


unsigned int Object::getNumPoints() {
    return static_cast<unsigned int>(this->points.size());
}

vec4 Object::getPosition(unsigned int num) {
    if (num >= this->points.size()) {
        return vec4(0);
    }
    return this->points[num];
}


double Object::getAffineParam(unsigned int num) {
    if (num >= this->lambda.size()) {
        return 0.0;
    }
    return this->lambda[num];
}


/**
 * @brief Object::clearAll
 */
void Object::clearAll() {
    if (!points.empty()) {
        points.clear();
    }
    if (!dirs.empty()) {
        dirs.clear();
    }
    if (!lambda.empty()) {
        lambda.clear();
    }
    if (!sachs1.empty()) {
        sachs1.clear();
    }
    if (!sachs2.empty()) {
        sachs2.clear();
    }
    if (!jacobi.empty()) {
        jacobi.clear();
    }

    for(unsigned int i=0; i < 4; i++) {
        if (!trans_lt[i].empty()) {
            trans_lt[i].clear();
        }
    }
}

/**
 * @brief Reset all parameters.
 */
void Object::resetAll() {
    if (currMetric != NULL) {
        delete currMetric;
        currMetric = NULL;
    }

    if (geodSolver != NULL) {
        delete geodSolver;
    }
    geodSolver = NULL;

    timeDirection = 1;
    tetradType    = enum_nat_tetrad_default;
    maxNumPoints  = 3000;

    stepsizeControlled = false;
    stepsize = 0.01;
    max_stepsize = DEF_MAX_STEPSIZE;
    min_stepsize = DEF_MIN_STEPSIZE;

    epsAbs    = 1.0e-6;
    epsRel    = 0.0;
    epsConstr = DEF_CONSTRAINT_EPSILON;
    epsResize = DEF_RESIZE_EPSILON;

    axes_orient = 0;
    ksi =  0.0;
    chi = 90.0;
    vel =  0.99;

    isBaseInCoords = false;
    base[0] = vec4(1, 0, 0, 0);
    base[1] = vec4(0, 1, 0, 0);
    base[2] = vec4(0, 0, 1, 0);
    base[3] = vec4(0, 0, 0, 1);


    boost_ksi  =  0.0;
    boost_chi  = 90.0;
    boost_beta =  0.0;
    lorentz.setIdent();

    geodSolverType = gsIrk4;
    type =  enum_geodesic_lightlike;

    clearAll();

    // geometric units
    speed_of_light  = 1.0;
    grav_constant   = 1.0;
    dielectric_perm = 1.0;
}


/**
 * @brief Get parameter value of Object.
 * @param paramName
 * @param paramValue
 * @return true : if parameter was found.\n
 *         false : parameter was not found.
 */
bool Object::getParam(const char* paramName, int &paramValue) {
    if (strcmp(paramName,"axes_orient") == 0) {
        paramValue = axes_orient;
        return true;
    }
    else if (strcmp(paramName,"timeDirection") == 0) {
        paramValue = timeDirection;
        return true;
    }
    else if (strcmp(paramName,"maxNumPoints") == 0) {
        paramValue = (int)maxNumPoints;
        return true;
    }
    return false;
}


/**
 * @brief Get parameter value of Object.
 * @param paramName
 * @param paramValue
 * @return true : if parameter was found.\n
 *         false : parameter was not found.
 */
bool Object::getParam(const char* paramName, double &paramValue) {
    if (strcmp(paramName,"stepsize")==0) {
        paramValue = stepsize;
        return true;
    }
    else if (strcmp(paramName,"max_stepsize")==0) {
        paramValue = max_stepsize;
        return true;
    }
    else if (strcmp(paramName,"min_stepsize")==0) {
        paramValue = min_stepsize;
        return true;
    }
    else if (strcmp(paramName,"epsAbs")==0) {
        paramValue = epsAbs;
        return true;
    }
    else if (strcmp(paramName,"epsRel")==0) {
        paramValue = epsRel;
        return true;
    }
    else if (strcmp(paramName,"epsConstr")==0) {
        paramValue = epsConstr;
        return true;
    }
    else if (strcmp(paramName,"epsResize")==0) {
        paramValue = epsResize;
        return true;
    }
    else if (strcmp(paramName,"ksi")==0) {
        paramValue = ksi;
        return true;
    }
    else if (strcmp(paramName,"chi")==0) {
        paramValue = chi;
        return true;
    }
    else if (strcmp(paramName,"vel")==0) {
        paramValue = vel;
        return true;
    }
    else if (strcmp(paramName,"boost_ksi")==0) {
        paramValue = boost_ksi;
        return true;
    }
    else if (strcmp(paramName,"boost_chi")==0) {
        paramValue = boost_chi;
        return true;
    }
    else if (strcmp(paramName,"boost_beta")==0) {
        paramValue = boost_beta;
        return true;
    }
    return false;
}


/**
 * @brief Get parameter value of Object.
 * @param paramName
 * @param paramValue
 * @return true : if parameter was found.\n
 *         false : parameter was not found.
 */
bool Object::getParam(const char* paramName, m4d::vec3 &paramValue) {
    if (strcmp(paramName,"startDir")==0) {
        paramValue = startDir;
        return true;
    }
    return false;
}


/**
 * @brief Get parameter value of Object.
 * @param paramName
 * @param paramValue
 * @return true : if parameter was found.\n
 *         false : parameter was not found.
 */
bool Object::getParam(const char* paramName, m4d::vec4 &paramValue) {
    if (strcmp(paramName,"startPos")==0) {
        paramValue = startPos;
        return true;
    }
    else if (strcmp(paramName,"base0")==0) {
        paramValue = base[0];
        return true;
    }
    else if (strcmp(paramName,"base1")==0) {
        paramValue = base[1];
        return true;
    }
    else if (strcmp(paramName,"base2")==0) {
        paramValue = base[2];
        return true;
    }
    else if (strcmp(paramName,"base3")==0) {
        paramValue = base[3];
        return true;
    }
    else if (strcmp(paramName,"coordDir")==0) {
        paramValue = coordDir;
        return true;
    }
    return false;
}


/**
 * @brief Set Lorentz transformation.
 *
 *  \param  chi : angle in deg.
 *  \param  ksi : angle in deg.
 *  \param  beta : velocity (v/c).
 */
bool Object::setLorentzTransf(const double chi, const double ksi, const double beta) {
    boost_ksi  = ksi;
    boost_chi  = chi;
    boost_beta = beta;

    if (fabs(beta) < 1.0) {
        resetLorentzTransf();
        double n[3] = {sin(chi * DEG_TO_RAD)*cos(ksi * DEG_TO_RAD),
                       sin(chi * DEG_TO_RAD)*sin(ksi * DEG_TO_RAD),
                       cos(chi * DEG_TO_RAD)
                      };
        double gamma = 1.0 / sqrt(1.0 - beta * beta);

        lorentz.setElem(0, 0, gamma);
        for (int row = 1; row < 4; row++) {
            for (int col = 1; col < 4; col++) {
                lorentz.setElem(row, col, (gamma - 1.0)*n[row - 1]*n[col - 1] + M4D_DELTA(row, col));
            }
            lorentz.setElem(0, row, beta * gamma * n[row - 1]);
            lorentz.setElem(row, 0, beta * gamma * n[row - 1]);
        }
        return true;
    }
    return false;
}


bool Object::setLorentzTransf(const m4d::vec3 beta) {
    double bn = beta.getNorm();
    if (bn >= 1.0) {
        return false;
    }
    m4d::vec3 n = beta.getNormalized();

    resetLorentzTransf();
    double gamma = 1.0 / sqrt(1.0 - bn * bn);

    lorentz.setElem(0, 0, gamma);
    for (int row = 1; row < 4; row++) {
        for (int col = 1; col < 4; col++) {
            lorentz.setElem(row, col, (gamma - 1.0) * n[row - 1] * n[col - 1] + M4D_DELTA(row, col));
        }
        lorentz.setElem(0, row, bn * gamma * n[row - 1]);
        lorentz.setElem(row, 0, bn * gamma * n[row - 1]);
    }
    return true;
}

/**
 * @brief Reset Lorentz transformation.
 */
void Object::resetLorentzTransf() {
    lorentz.setIdent();
}

/**
 * @brief Load settings.
 *
 *  \param filename : name of setting file.
 *  \param printset : print setting.
 *  \return true : success.
 *  \return false : error occured.
 */
bool Object::loadSettings(const char* filename, bool printset) {
    setlocale(LC_NUMERIC, "C");
    std::vector<std::vector<std::string> >  tokens;
    m4d::tokenizeFile(std::string(filename), tokens);

    bool ok = true;
    for (unsigned int i = 0; i < tokens.size(); i++) {
        if (tokens[i].size() == 0) {
            continue;
        }

        std::string baseString = tokens[i][0];
        if (baseString.compare("METRIC") == 0 && tokens[i].size() > 1) {
            currMetric = metricDB.getMetric(tokens[i][1].c_str());
            if (currMetric == NULL) {
                return false;
            }
            ok &= true;
        }
        else if (baseString.compare("PARAM") == 0 && tokens[i].size() > 3) {
            //int    pnum  = atoi(tokens[i][1].c_str());
            std::string pname = tokens[i][2];
            double val   = atof(tokens[i][3].c_str());
            currMetric->setParam(pname.c_str(), val);
        }
        else if (baseString.compare("INIT_POS") == 0 && tokens[i].size() > 4) {
            for (int j = 0; j < 4; j++) {
                startPos[j] = atof(tokens[i][j + 1].c_str());
            }
        }
        else if (baseString.compare("INIT_DIR") == 0 && tokens[i].size() > 3) {
            for (int j = 0; j < 3; j++) {
                startDir[j] = atof(tokens[i][j + 1].c_str());
            }
        }
        else if (baseString.compare("INIT_ANGLE_VEL") == 0 && tokens[i].size() > 3) {
            ksi = atof(tokens[i][1].c_str());
            chi = atof(tokens[i][2].c_str());
            vel = atof(tokens[i][3].c_str());

            startDir[0] = sin(chi * DEG_TO_RAD) * cos(ksi * DEG_TO_RAD);
            startDir[1] = sin(chi * DEG_TO_RAD) * sin(ksi * DEG_TO_RAD);
            startDir[2] = cos(chi * DEG_TO_RAD);
        }
        else if (baseString.compare("TIME_DIR") == 0 && tokens[i].size() > 1) {
            timeDirection = atoi(tokens[i][1].c_str());
        }
        else if (baseString.compare("AXES_ORIENT") == 0 && tokens[i].size() > 1) {
            axes_orient = atoi(tokens[i][1].c_str());
        }
        else if (baseString.compare("GEOD_SOLVER_TYPE") == 0 && tokens[i].size() > 1) {
            geodSolverType = enum_integrator(atoi(tokens[i][1].c_str()));
            if (currMetric != NULL) {
                geodSolver = solverDB.getIntegrator(currMetric, geodSolverType);
            }
        }
        else if (baseString.compare("GEODESIC_TYPE") == 0 && tokens[i].size() > 1) {
            type = enum_geodesic_lightlike;
            unsigned int j = 0;
            while (j < NUM_ENUM_GEODESIC_TYPE && (tokens[i][1].compare(stl_geodesic_type[j]) != 0)) {
                j++;
            }
            if (j < NUM_ENUM_GEODESIC_TYPE) {
                type = (enum_geodesic_type)j;
            } else {
                fprintf(stderr, "m4dObject::loadSettings() ... geodesic type not recognized! Please check ini-file.");
            }
        }
        else if (baseString.compare("STEPSIZE_CTRL") == 0 && tokens[i].size() > 1) {
            stepsizeControlled = (atoi(tokens[i][1].c_str())) == 1 ? true : false;
        }
        else if (baseString.compare("STEPSIZE") == 0 && tokens[i].size() > 1) {
            stepsize = atof(tokens[i][1].c_str());
        }
        /*
        else if (baseString.compare("STEPSIZE_MAX_MIN")==0 && tokens[i].size()>2)
        {
        max_stepsize = atof(tokens[i][1].c_str());
        min_stepsize = atof(tokens[i][2].c_str());
        }
        */
        else if (baseString.compare("STEPSIZE_MAX") == 0 && tokens[i].size() > 1) {
            max_stepsize = atof(tokens[i][1].c_str());
            //min_stepsize = atof(tokens[i][2].c_str());
        }
        else if (baseString.compare("EPSILONS") == 0 && tokens[i].size() > 2) {
            epsAbs = atof(tokens[i][1].c_str());
            epsRel = atof(tokens[i][2].c_str());
        }
        else if (baseString.compare("CONSTR_EPSILON") == 0 && tokens[i].size() > 1) {
            epsConstr = atof(tokens[i][1].c_str());
        }
        // else if (baseString.compare("RESIZE_EPSILON")==0 && tokens[i].size()>1)
        //   epsResize = atof(tokens[i][1].c_str());
        else if (baseString.compare("MAX_NUM_POINTS") == 0 && tokens[i].size() > 1) {
            maxNumPoints = atoi(tokens[i][1].c_str());
        }
        else if (baseString.compare("TETRAD_TYPE") == 0 && tokens[i].size() > 1) {
            tetradType = enum_nat_tetrad_type(atoi(tokens[i][1].c_str()));
        }
        else if (baseString.compare("BASE_0") == 0 && tokens[i].size() > 4) {
            for (int j = 0; j < 4; j++) {
                base[0][j] = atof(tokens[i][j + 1].c_str());
            }
        }
        else if (baseString.compare("BASE_1") == 0 && tokens[i].size() > 4) {
            for (int j = 0; j < 4; j++) {
                base[1][j] = atof(tokens[i][j + 1].c_str());
            }
        }
        else if (baseString.compare("BASE_2") == 0 && tokens[i].size() > 4) {
            for (int j = 0; j < 4; j++) {
                base[2][j] = atof(tokens[i][j + 1].c_str());
            }
        }
        else if (baseString.compare("BASE_3") == 0 && tokens[i].size() > 4) {
            for (int j = 0; j < 4; j++) {
                base[3][j] = atof(tokens[i][j + 1].c_str());
            }
        }
        else if (baseString.compare("BOOST") == 0 && tokens[i].size() > 3) {
            boost_ksi = atof(tokens[i][1].c_str());
            boost_chi = atof(tokens[i][2].c_str());
            boost_beta = atof(tokens[i][3].c_str());
            setLorentzTransf(boost_chi, boost_ksi, boost_beta);
        }
        else if (baseString.compare("SPEED_OF_LIGHT") == 0 && tokens[i].size() > 1) {
            speed_of_light = atof(tokens[i][1].c_str());
        }
        else if (baseString.compare("GRAV_CONSTANT") == 0 && tokens[i].size() > 1) {
            grav_constant = atof(tokens[i][1].c_str());
        }
        else if (baseString.compare("DIELECTRIC_PERM") == 0 && tokens[i].size() > 1) {
            dielectric_perm = atof(tokens[i][1].c_str());
        }
    }

    if (printset && ok) {
        printSettings();
        currMetric->printF();
    }

    if (geodSolver != NULL) {
        geodSolver->setGeodesicType(type);
        geodSolver->setEpsilons(epsAbs, epsRel);
        geodSolver->setStepSizeControlled(stepsizeControlled);
        geodSolver->setAffineParamStep(stepsize);
    }
    return ok;
}

/*!  Save settings.
 *
 *  \param filename : name of the settings file.
 *  \param dat      : current date.
 *  \return true : success.
 *  \return false : error occured.
 */
bool Object::saveSettings(const char* filename, const char *dat) {
    if (currMetric == NULL) {
        return false;
    }

    FILE* fptr;
#ifdef _WIN32
    fopen_s(&fptr, filename, "w");
#else
    fptr = fopen(filename, "w");
#endif
    if (fptr == NULL) {
        fprintf(stderr, "Cannot open %s for output!\n", filename);
        return false;
    }

    char* oldlocale = setlocale(LC_NUMERIC, "C");

    fprintf(fptr, "# ----------------------------------------------------------\n");
    fprintf(fptr, "# Settings file : %s\n", filename);
    if (dat != nullptr) {
        fprintf(fptr, "#          date : %s\n", dat);
    }
    fprintf(fptr, "# ----------------------------------------------------------\n");

    fprintf(fptr, "METRIC       %s\n", currMetric->getMetricName());
    double val;

    char* pname = nullptr;
    for (int i = 0; i < currMetric->getNumParams(); i++) {
        currMetric->getParam(i, pname, val);
        fprintf(fptr, "PARAM  %d  %10s  %16.12f\n", i, pname, val);
    }
    if (pname != nullptr) {
        delete [] pname;
    }
    
    fprintf(fptr, "INIT_POS         %18.14f %18.14f %18.14f %18.14f\n", startPos[0], startPos[1], startPos[2], startPos[3]);
    fprintf(fptr, "INIT_DIR         %18.14f %18.14f %18.14f\n", startDir[0], startDir[1], startDir[2]);
    fprintf(fptr, "INIT_ANGLE_VEL   %18.14f %18.14f %18.14f\n", ksi, chi, vel);
    fprintf(fptr, "TIME_DIR          %d\n", timeDirection);
    fprintf(fptr, "AXES_ORIENT       %d\n", axes_orient);
    fprintf(fptr, "GEOD_SOLVER_TYPE  %d\n", int(geodSolverType));
    fprintf(fptr, "GEODESIC_TYPE     %s\n", stl_geodesic_type[type]);
    if (stepsizeControlled) {
        fprintf(fptr, "STEPSIZE_CTRL     1\n");
    } else {
        fprintf(fptr, "STEPSIZE_CTRL     0\n");
    }
    fprintf(fptr, "STEPSIZE          %16.12e\n", stepsize);
    //fprintf(fptr,"STEPSIZE_MAX_MIN  %16.12e %16.12e\n",max_stepsize,min_stepsize);
    fprintf(fptr, "STEPSIZE_MAX      %16.12e \n", max_stepsize);
    fprintf(fptr, "EPSILONS          %16.12e %16.12e\n", epsAbs, epsRel);
    fprintf(fptr, "CONSTR_EPSILON    %16.12e\n", epsConstr);
    // fprintf(fptr,"RESIZE_EPSILON    %16.12e\n",epsResize);
    fprintf(fptr, "MAX_NUM_POINTS    %d\n", maxNumPoints);
    fprintf(fptr, "TETRAD_TYPE       %d\n", int(tetradType));
    for (int i = 0; i < 4; i++) {
        fprintf(fptr, "BASE_%d          %16.12f %16.12f %16.12f %16.12f\n", i, base[i][0], base[i][1], base[i][2], base[i][3]);
    }
    fprintf(fptr, "BOOST           %16.12f %16.12f %16.12f\n", boost_ksi, boost_chi, boost_beta);
    fprintf(fptr, "SPEED_OF_LIGHT    %12.6f\n", speed_of_light);
    fprintf(fptr, "GRAV_CONSTANT     %12.6e\n", grav_constant);
    fprintf(fptr, "DIELECTRIC_PERM   %12.6e\n", dielectric_perm);
    fclose(fptr);

    setlocale(LC_NUMERIC, oldlocale);
    return true;
}

/*! Print settings to fptr.
 *
 *  \param fptr : pointer to file.
 */
void  Object::printSettings(FILE* fptr) {
    if (currMetric == NULL || geodSolver == NULL) {
        fprintf(fptr, "Error in Object::printSettings() ... no metric or solver set!\n");
        return;
    }

    char* oldlocale = setlocale(LC_NUMERIC, "C");

    fprintf(fptr, "\n--------------------- parameter settings ---------------------\n");
    fprintf(fptr, "METRIC       %s\n", currMetric->getMetricName());
    double val = 0.0;

    for (int i = 0; i < currMetric->getNumParams(); i++) {
        const char* pname = currMetric->getParamName(i);
        if (pname != nullptr && currMetric->getParam(pname, val)) {
            fprintf(fptr, "PARAM  %d  %10s  %16.12f\n", i, pname, val);
        }
    }
    fprintf(fptr, "INIT_POS         %16.12f %16.12f %16.12f %16.12f\n", startPos[0], startPos[1], startPos[2], startPos[3]);
    fprintf(fptr, "INIT_DIR         %16.12f %16.12f %16.12f\n", startDir[0], startDir[1], startDir[2]);
    fprintf(fptr, "INIT_ANGLE_VEL   %16.12f %16.12f %16.12f\n", ksi, chi, vel);
    fprintf(fptr, "TIME_DIR          %s\n", str_time_dir[(timeDirection + 1) / 2]);
    fprintf(fptr, "AXES_ORIENT       %d\n", axes_orient);
    fprintf(fptr, "GEOD_SOLVER_TYPE  %s\n", stl_solver_names[geodSolverType]);
    //cerr << kappa << endl;
    fprintf(fptr, "GEODESIC_TYPE     %s\n", stl_geodesic_type[type]);
    if (stepsizeControlled) {
        fprintf(fptr, "STEPSIZE_CTRL     yes\n");
    } else {
        fprintf(fptr, "STEPSIZE_CTRL     no\n");
    }
    fprintf(fptr, "STEPSIZE          %16.12e\n", stepsize);
    //fprintf(fptr,"STEPSIZE_MAX_MIN  %16.12e %16.12e\n",max_stepsize,min_stepsize);
    fprintf(fptr, "STEPSIZE_MAX      %16.12e \n", max_stepsize);
    fprintf(fptr, "EPSILONS          %16.12e %16.12e\n", epsAbs, epsRel);
    fprintf(fptr, "CONSTR_EPSILON    %16.12e\n", epsConstr);
    // fprintf(fptr,"RESIZE_EPSILON    %16.12e\n",epsResize);
    fprintf(fptr, "MAX_NUM_POINTS    %d\n", maxNumPoints);
    fprintf(fptr, "TETRAD_TYPE       %s\n", stl_draw_type[tetradType]);
    for (int i = 0; i < 4; i++) {
        fprintf(fptr, "BASE_%d          %16.12f %16.12f %16.12f %16.12f\n", i, base[i][0], base[i][1], base[i][2], base[i][3]);
    }
    fprintf(fptr, "BOOST             %16.12f %16.12f %16.12f\n", boost_ksi, boost_chi, boost_beta);
    fprintf(fptr, "SPEED_OF_LIGHT    %16.12f\n", speed_of_light);
    fprintf(fptr, "GRAV_CONSTANT     %16.12e\n", grav_constant);
    fprintf(fptr, "DIELECTRIC_PERM   %16.12e\n", dielectric_perm);
    setlocale(LC_NUMERIC, oldlocale);
}

/*! Prepare a report for the current metric.
 *
 *  \param text : reference to string.
 *  \return true : success.
 *  \return false : no metric available.
 */
bool Object::makeReport(std::string  &text) {
    if (currMetric == NULL) {
        return false;
    }

    double eta = 1.0;
    double y0  = 1.0;

    if (type == enum_geodesic_timelike && fabs(vel) < 1.0) {
        y0  = 1.0 / sqrt(1.0 - vel * vel);
        eta = vel * y0;
    }
    vec4 locDir = M4D_SIGN(timeDirection) * y0 * base[0] + eta * (startDir[0] * base[1] + startDir[1] * base[2] + startDir[2] * base[3]);
    vec4 coDir;
    currMetric->localToCoord(startPos, locDir, coDir, tetradType);

  //  currMetric->report(startPos, coDir, text);
    return true;
}

void Object::printReport(FILE* fptr) {
    std::string text;
   // makeReport(text);
    fprintf(fptr, "%s\n", text.c_str());
}

} // end namespace m4d


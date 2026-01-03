// Copyright (C) 2025
// Laboratoire d'Informatique Fondamentale et Appliquée de Tours, Tours, France
//
// DIGEP, Politecnico di Torino, Corso Duca degli Abruzzi 24, Torino, Italy
// This file is part of F2-Cmax-revisited.
//
// F2-Cmax-revisited is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// F2-Cmax-revisited is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with F2-Cmax-revisited. If not, see <https://www.gnu.org/licenses/>.

//
// Created by schau on 12/2/25.
//

#ifndef CODE_INSTANCE_H
#define CODE_INSTANCE_H

#include <string>
#include <fstream>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <random>
#include <vector>

#include "F2CmaxException.h"

enum LAW{UNIFORM, NEGATIVE_BINOMIAL, GEOMETRIC,POISSON,EXPONENTIAL, GAMMA,WEIBULL,LOGNORMAL }; // law uses to generates processing times
class Instance {
public:
    typedef std::pair<double,double> Job;

private:
    double supPj = 100.0;
    std::string instanceName;
    std::filesystem::path instancePath; // the path to the instance
    unsigned int nbJobs=0; // the nb of job
    std::vector<Job> jobsSmallerOnM1;
    std::vector<Job> jobsSmallerOnM2;
    // the seed use for generate instance
    std::mt19937 numGenerator;
    LAW distribution=UNIFORM;
    double p_max_A = 0.0;
    double p_max_B = 0.0;
    double p_max = 0.0;
    double sumPA1 = 0.0;
    double sumPA2 = 0.0;
    double sumPB1 = 0.0;
    double sumPB2 = 0.0;

public:
    /**
         * Default Constructor
         * Set seed at 0 by default
         */
    explicit Instance();

    /**
     * Constructor by instance's path. It constructs an instance with setting the attribute path. Moreover,
     * the method check if the path have a corresponding file, otherwise an exception is throw
     * @param newInstancePath The path to set
     */
    explicit Instance(std::string &newInstancePath);

    void clearListJobs() {
        jobsSmallerOnM1.clear(); jobsSmallerOnM2.clear();
    }
    void addJob(double pi1, double pi2) {
        p_max = std::max(p_max,std::max(pi1,pi2));
        if (pi1<pi2) {
            p_max_A = std::max(p_max_A,pi1);
            sumPA1 += pi1;
            sumPA2 += pi2;
            jobsSmallerOnM1.emplace_back(pi1,pi2);
        }
        else {
            p_max_B = std::max(p_max_B,pi2);
            sumPB1 += pi2;
            sumPB2 += pi1;
            jobsSmallerOnM2.emplace_back(pi2,pi1); //add the job directly regarding the reverse property
        }
    }

    Job generateJob(unsigned int infPi, unsigned int supPi);

    void generateInstance(nlohmann::json &paramInstance);

    /********************/
    /*      GETTER      */
    /********************/

    [[nodiscard]] const std::filesystem::path &getInstancePath() const { return instancePath; }

    [[nodiscard]] const std::string &getInstanceName() const { return instanceName; }

    [[nodiscard]] unsigned int getNbJobs() const { return nbJobs; }
    [[nodiscard]] double getSupPj() const { return supPj; }

    [[nodiscard]] double getPMaxA() const { return p_max_A; }
    [[nodiscard]] double getPMaxB() const { return p_max_B; }
    [[nodiscard]] double getPMax() const { return p_max; }

    [[nodiscard]] double getSumPa1() const { return sumPA1; }
    [[nodiscard]] double getSumPa2() const { return sumPA2; }
    [[nodiscard]] double getSumPb1() const { return sumPB1; }
    [[nodiscard]] double getSumPb2() const { return sumPB2; }

    [[nodiscard]] std::vector<Job> & getJobsSmallerOnM1() { return jobsSmallerOnM1; }
    [[nodiscard]] std::vector<Job> & getJobsSmallerOnM2() { return jobsSmallerOnM2; }

    [[nodiscard]] std::vector<Job> getJobsSmallerOnM1() const { return jobsSmallerOnM1; }
    [[nodiscard]] std::vector<Job> getJobsSmallerOnM2() const { return jobsSmallerOnM2; }

    /********************/
    /*      SETTER      */
    /********************/

    void setInstanceName(const std::string &newInstanceName) {
        instanceName = newInstanceName;
        size_t pos1 = instanceName.find("_pmax_");
        if (pos1 != std::string::npos) {
            size_t pos2 = instanceName.find('_', pos1 + 5);
            if (pos2 != std::string::npos) {
                supPj = std::stod(instanceName.substr(pos1 + 6, pos2 - pos1 - 6));
            } else {
                throw F2CmaxException("No '*_pmax_*' was found in the name of the instance");
            }
        } else {
            throw F2CmaxException("No '*_pmax_*' was found in the name of the instance");
        }
    }

    void setInstancePath(const std::string &newInstancePath) {
        std::filesystem::path newPath = std::filesystem::path(newInstancePath);
        instancePath = newPath;
        std::filesystem::directory_entry parentDir{newPath.lexically_normal().parent_path()};
        if (std::filesystem::exists(newPath)) instancePath = newPath;
        else {
            // if the dir doesn't exist
            if (!parentDir.exists()) {
                std::filesystem::create_directories(newPath.lexically_normal().parent_path());
            }
        }
        instanceName = instancePath.stem();
    }

    void setSeed(unsigned int seed) { numGenerator = std::mt19937(seed); }

    void setNbJobs(unsigned int nbJobs) {
        this->nbJobs = nbJobs;
        jobsSmallerOnM1.reserve(nbJobs);
        jobsSmallerOnM2.reserve(nbJobs);
    }
};

inline std::ostream &operator<<(std::ostream &os, const Instance &instance) {

    os << "M1:[";
    for (auto &job : instance.getJobsSmallerOnM1()) os << job.first << " ";
    for (auto &job : std::ranges::reverse_view(instance.getJobsSmallerOnM2())) os << job.second << " ";
    os << "]" << std::endl << "M2:[" ;
    for (auto &job : instance.getJobsSmallerOnM1()) os << job.second << " ";
    for (auto &job : std::ranges::reverse_view(instance.getJobsSmallerOnM2())) os << job.first << " ";
    os << "]" << std::endl << "----------";
    return os;
}


inline std::ostream &operator<<(std::ostream &os, const std::vector<Instance::Job> &vector) {

    os << "[";
    if (vector.empty()) os << "]";
    else {
        for (size_t indexLoopVector = 0; indexLoopVector < vector.size() - 1; ++indexLoopVector) {
            auto &element = vector[indexLoopVector];
            os << element.first << ",";
        }
        os << vector.back().first << "]" << std::endl << "[";
        for (size_t indexLoopVector = 0; indexLoopVector < vector.size() - 1; ++indexLoopVector) {
            auto &element = vector[indexLoopVector];
            os << element.second << ",";
        }
        os << vector.back().second << "]" << std::endl << "----------";
    }
    return os;
}

template<typename T>
inline std::ostream &operator<<(std::ostream &os, const std::vector<T> &vector) {

    os << "[";
    if (vector.empty()) os << "]";
    else {
        for (size_t indexLoopVector = 0; indexLoopVector < vector.size() - 1; ++indexLoopVector) {
            auto &element = vector[indexLoopVector];
            os << element << ",";
        }
        os << vector.back() << "]";
    }
    return os;
}

#endif //CODE_INSTANCE_H
// Copyright (C) 2024
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

#include "Instance.h"


Instance::Instance() {
    std::random_device rd;//random number engine
    setSeed(rd());
}

Instance::Instance(std::string& newInstancePath) {
    setInstancePath(newInstancePath);
    std::random_device rd;//random number engine
    setSeed(rd());
}

Instance::Job Instance::generateJob(unsigned int infPi, unsigned int supPi) {
    double pj1,pj2;
    switch (distribution) {
        case UNIFORM: {
            std::uniform_int_distribution<> piDistribution(infPi, supPi);
            pj1= static_cast<double>(piDistribution(numGenerator));
            pj2= static_cast<double>(piDistribution(numGenerator));
            break;
        }
        case NEGATIVE_BINOMIAL: {
            constexpr unsigned int r=5u;
            const double p = static_cast<double>(r)/(static_cast<double>(r)+supPi/2.0);
            std::negative_binomial_distribution<> piDistribution(r,p);
            pj1= static_cast<double>(piDistribution(numGenerator));
            pj2= static_cast<double>(piDistribution(numGenerator));
            break;
        }case GEOMETRIC: {
            const double p = 2.0/(2.0+static_cast<double>(supPi));
            std::geometric_distribution<> piDistribution(p);
            pj1= static_cast<double>(piDistribution(numGenerator));
            pj2= static_cast<double>(piDistribution(numGenerator));
            break;
        }case POISSON: {
            const double p = supPi/2.0;
            std::poisson_distribution<> piDistribution(p);
            pj1= static_cast<double>(piDistribution(numGenerator));
            pj2= static_cast<double>(piDistribution(numGenerator));
            break;
        }case EXPONENTIAL: {
            const double p = 2.0/supPi;
            std::exponential_distribution<> piDistribution(p);
            pj1= piDistribution(numGenerator);
            pj2= piDistribution(numGenerator);
            break;
        }case GAMMA: {
            constexpr double k= 2.0;
            const double theta = supPi / 4.0;
            std::gamma_distribution<> piDistribution(k,theta);
            pj1 = piDistribution(numGenerator);
            pj2 = piDistribution(numGenerator);
            break;
        }case WEIBULL: {
            constexpr double k = 1.5;
            const double lambda = supPi / 2.0 / std::tgamma(1.0 + 1.0 / k);
            std::weibull_distribution<> piDistribution(k, lambda);
            pj1 = piDistribution(numGenerator);
            pj2 = piDistribution(numGenerator);
            break;
        }case LOGNORMAL: {
            constexpr double sigma = 0.5;
            const double mu = std::log(supPi/2.0)- sigma*sigma/2.0;
            std::lognormal_distribution<> piDistribution(mu,sigma);
            pj1 = piDistribution(numGenerator);
            pj2 = piDistribution(numGenerator);
            break;
        }
        default: throw F2CmaxException("distribution law not implemented");
    }
    return {pj1,pj2};
}

void Instance::generateInstance(nlohmann::json& paramInstance) {
    if (paramInstance["paramInstance"].contains("distribution")) {
        if (paramInstance["paramInstance"]["distribution"].is_string()) {
            auto distributionName = paramInstance["paramInstance"]["distribution"].get<std::string>();
            if (distributionName == "uniform"){
                distribution = UNIFORM;
            }else if (distributionName == "negative_binomial") {
                distribution = NEGATIVE_BINOMIAL;
            }else if (distributionName == "geometric") {
                distribution=GEOMETRIC;
            }else if (distributionName == "poisson") {
                distribution=POISSON;
            }else if (distributionName == "exponential") {
                distribution=EXPONENTIAL;
            }else if (distributionName == "gamma") {
                distribution=GAMMA;
            }else if (distributionName == "weibull") {
                distribution=WEIBULL;
            }else if (distributionName == "lognormal") {
                distribution=LOGNORMAL;
            }else throw F2CmaxException("The distribution law is not implemented");
        }
        else throw std::invalid_argument(R"(The "distribution" must be an string object)");
    }
    // set the number of N jobs
    if (paramInstance.contains("n")) {
        if (paramInstance["n"].is_number_unsigned()) setNbJobs(paramInstance["n"]);
        else throw std::invalid_argument(R"(The "N" must be an unsigned integer)");
    }

    //set pi distribution
    unsigned int infPi = 1;
    unsigned int supPi = 100;

    if (paramInstance.contains("pi")) {
        if (paramInstance["pi"].contains("inf")) {
            if (paramInstance["pi"]["inf"].is_number_unsigned()) infPi = paramInstance["pi"]["inf"];
            else throw std::invalid_argument(R"(The "inf" must be an unsigned integer in the "pi" object)");
        }
        if (paramInstance["pi"].contains("sup")) {
            if (paramInstance["pi"]["sup"].is_number_unsigned()) supPi = paramInstance["pi"]["sup"];
            else throw std::invalid_argument(R"(The "sup" must be an unsigned integer in the "pi" object)");
        }
    }

    // generate Jobs
    jobsSmallerOnM1.reserve(nbJobs);
    jobsSmallerOnM2.reserve(nbJobs);
    for (unsigned int i = 0; i < nbJobs; ++i) {
        auto newJob = generateJob(infPi, supPi);
        addJob(newJob.first,newJob.second);
    }
}

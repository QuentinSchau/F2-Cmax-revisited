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

#include "Parser.h"

Parser::Parser() {}

Instance Parser::readFromFile(std::string &filePath) const {
    Instance newInstance = Instance(filePath);
    std::fstream fileStream(newInstance.getInstancePath().lexically_normal(), std::fstream::in);
    std::string line; // new line
    // open and read the file
    if (fileStream.is_open()) {
        while (std::getline(fileStream, line)) {
            auto pos = line.find(":");
            // if we don't find the ':' then we read jobs
            if (pos == std::string::npos) {
                // create the job, the value should be separate with \t
                std::istringstream stream(line);
                double pj1, pj2;
                stream >> pj1;
                stream >> pj2;
                newInstance.addJob(pj1,pj2);
            } else {
                // we read so attributes
                std::string attribute = line.substr(0, pos);
                if (attribute == "name") newInstance.setInstanceName(line.substr(pos + 1));
                else if (attribute == "n") newInstance.setNbJobs(std::stoul(line.substr(pos + 1)));
            }
        }
    } else
        throw F2CmaxException(std::string("Can't open the file ").append(newInstance.getInstancePath().lexically_normal().string()).c_str());
    fileStream.close();

    // check if we have the right number of created job
    if (newInstance.getNbJobs() != newInstance.getJobsSmallerOnM1().size()+newInstance.getJobsSmallerOnM2().size())
        throw std::invalid_argument("The number of jobs is not equals to n");
    return newInstance;
}

void Parser::serializeInstance(Instance &instance) {
    std::fstream fileStream(instance.getInstancePath().lexically_normal().string(), std::fstream::out );
    if (fileStream.is_open()) {
        fileStream << std::setprecision(5) << "name:" << instance.getInstanceName() << std::endl
                   << "n:" << instance.getNbJobs() << std::endl << "Jobs:" << std::endl;
        for (auto &job: instance.getJobsSmallerOnM1()) {
            fileStream << job.first << "\t" << job.second << std::endl;
        }for (auto &job: instance.getJobsSmallerOnM2()) {
            fileStream << job.second << "\t" << job.first << std::endl;
        }
    } else throw F2CmaxException(std::string("Can't open the file ").append(instance.getInstancePath().lexically_normal().string()).c_str());
    fileStream.close();
}

void Parser::generateInstance(nlohmann::json &object) {

    Instance newInstance;

    // set the seed for generate
    if (object.contains("seed")) {
        if (object["seed"].is_number_unsigned()) newInstance.setSeed(object["seed"]);
        else throw std::invalid_argument(R"(The seed is not a unsigned int)");
    }

    // loop over each instances
    if (object.contains("instances")) {
        // create each kind of instance
        for (auto &paramInstance: object["instances"]) {
            // check if there is a base Path
            std::string basePath;
            if (paramInstance.contains("basePath")) {
                if (paramInstance["basePath"].is_string()) basePath = paramInstance["basePath"];
                else throw std::invalid_argument(R"(The base Path is not a string)");
            } else{
                basePath = std::filesystem::current_path().lexically_normal().string();
            }
            std::cout << "All instances will be generated at : " << basePath << std::endl;
            unsigned int nbGeneratedInstance = 0;
            unsigned int nbInstanceToGenerate = 1;
            // get the number of instance that we have to generate
            if (paramInstance.contains("numberInstance")) {
                if (paramInstance["numberInstance"].is_number_unsigned()) nbInstanceToGenerate = paramInstance["numberInstance"];
                else throw std::invalid_argument(R"(The "numberInstance" must be an unsigned integer)");
            }
            // we need to have the object "paramInstance"
            if (!paramInstance.contains("paramInstance")) throw std::invalid_argument(R"(The "paramInstance" is not defined)");
            unsigned int maxP = 100;
            if (paramInstance["paramInstance"].contains("pi") && paramInstance["paramInstance"]["pi"].contains("inf")) {
                if (paramInstance["paramInstance"]["pi"]["sup"].is_number_unsigned()) maxP = paramInstance["paramInstance"]["pi"]["sup"];
                else throw std::invalid_argument(R"(The "sup" must be an unsigned integer in the "pi" object)");
            }
            std::string distribution = "uniform"; // default distribution is uniform
            if (paramInstance["paramInstance"].contains("distribution")) {
                if (paramInstance["paramInstance"]["distribution"].is_string()) distribution = paramInstance["paramInstance"]["distribution"];
                else throw std::invalid_argument(R"(The "distribution" must be an string object)");
            }

            for (unsigned int newInstanceLoop = 0; newInstanceLoop < nbInstanceToGenerate; ++newInstanceLoop) {
                std::string path = basePath;
                path.append("instance")
                    .append(std::to_string(nbGeneratedInstance))
                    .append("_n_").append(std::to_string(paramInstance["paramInstance"]["n"].template get<unsigned int>()))
                    .append("_pmax_").append(std::to_string(maxP))
                    .append("_distribution_").append(distribution)
                    .append(".txt");
                newInstance.setInstancePath(path);
                newInstance.generateInstance(paramInstance["paramInstance"]);
                serializeInstance(newInstance);
                ++nbGeneratedInstance;
                newInstance.clearListJobs();
            }
        }
    }else throw std::invalid_argument(R"(The generate config JSON must have an instance object)");

}

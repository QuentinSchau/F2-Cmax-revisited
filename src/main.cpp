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

#include "Instance.h"
#include "Parser.h"
#include <iostream>
#include <nlohmann/json.hpp>

#include "Solver.h"


int main(int argc, char **argv) {

    if (argc < 2) {
        std::cerr << "ERROR: You need at least one argument." << std::endl;
        return -1;
    }
    char **pargv = argv + 1;
    try {
        for (; *pargv != argv[argc]; pargv++) {
            auto pathFile = std::filesystem::path(std::string(*pargv));
            if (! std::filesystem::exists(pathFile)) {
                throw F2CmaxException(std::string("The file configuration use do not exist. Path: ").append(pathFile.string()));
            }
            std::ifstream f(*pargv);
            // parse json input config file
            nlohmann::json config = nlohmann::json::parse(f);

            // create instance parser
            Parser parser = Parser();

            /********************************/
            /*      GENERATE INSTANCES      */
            /********************************/

            if (config.contains("generate")) {
                // generate them with Parser
                parser.generateInstance(config["generate"]);
            }

            /*****************************/
            /*      SOLVE INSTANCES      */
            /*****************************/

            if (config.contains("solve")) {
                // set verbose mode
                char verbose = config["solve"].contains("verbose") && config["solve"]["verbose"].is_number_unsigned()
                               ? config["solve"]["verbose"].template get<char>() : 0;

                if (config["solve"].contains("methods")) {
                    // for each method
                    for (auto &method: config["solve"]["methods"]) {
                        // set output path
                        std::string outputPath;
                        if (config["solve"].contains("output") && config["solve"]["output"].is_string())
                            outputPath = std::filesystem::path(config["solve"]["output"].template get<std::string>());
                        else {
                            outputPath = std::filesystem::current_path().parent_path().parent_path().string() +
                                         "/instances/";
                        }

                        // keep the path without the extension and add the name method;
                        outputPath.append("resultsF2Cmax.csv");
                        if (verbose >= 2) std::cout << "Save results in the path : " << outputPath << std::endl;
                        std::ofstream outputFileStream;

                        if (method.contains("instances")) {
                            // loop over each instances
                            for (auto &instance: method["instances"]) {
                                Instance newInstance;
                                if (instance.contains("path")) {
                                    if (instance["path"].is_string()) {
                                        std::string path = instance["path"];
                                        if (verbose >= 2) std::cout << "Parsing instance : " << path << std::endl;
                                        newInstance = parser.readFromFile(path);
                                        bool useRevisited = true;
                                        if (method.contains("useRevisited")) {
                                            if (method["useRevisited"].is_boolean()) useRevisited = method["useRevisited"].get<bool>();
                                            else throw std::invalid_argument(R"(The "useRevisited" must be an boolean)");
                                        }
                                        Solver solver(&newInstance,useRevisited);
                                        solver.solve();
                                        solver.printOutput(outputPath, outputFileStream);
                                    } else throw std::invalid_argument(R"(The instance path is not a string)");
                                } else throw std::invalid_argument(R"(The instance don't have attribute "path")");
                            }
                        }
                        outputFileStream.close();
                    }
                } else throw std::invalid_argument(R"(The config don't have attribute "methods")");
            }
        }
    }catch (const std::exception &e) {
        std::cerr << "Error with "<< *pargv << std::endl << "Error: " << e.what();
        return -1;
    }

    return 0;
}

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

#ifndef F2_CMAX_EXCEPTION_H
#define F2_CMAX_EXCEPTION_H

#include <exception>
#include <string>

class F2CmaxException : public std::exception {
private:
    std::string message;

public:
    explicit F2CmaxException(const char *msg) : message(msg) {}

    explicit F2CmaxException(std::string &msg) : message(msg) {}

    // Override the what() method to return our message
    virtual const char *what() const throw() {
        auto error = message.c_str();
        return error;
    }
};

#endif //F2_CMAX_EXCEPTION_H
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

#ifndef F2_CMAX_SOLVER_H
#define F2_CMAX_SOLVER_H

#include <random>

#include "Instance.h"

enum PIVOT_RULE{BFPRT};
class Solver {
    Instance * instance = nullptr;
    bool useRevisitedAlgo = true;
    std::chrono::duration<double> time_elapsed_johnson;
    std::chrono::duration<double> time_elapsed_revisited_johnson;
    enum SIDE{A,B};
    PIVOT_RULE pivotRule;
    // metrics where e have ppt1, k_a, ppt2, k_b
    std::tuple<size_t,size_t,size_t,size_t,size_t,size_t> metrics;
    double objective;
public:
    explicit Solver(Instance* instance,bool useRevisitedAlgo) : instance(instance),useRevisitedAlgo(useRevisitedAlgo), time_elapsed_johnson(0),time_elapsed_revisited_johnson(0), pivotRule(BFPRT) {}

    void solve() {
        if (useRevisitedAlgo) {
            // // shuffle list jobs to start from scratch
            std::shuffle(instance->getJobsSmallerOnM1().begin(), instance->getJobsSmallerOnM1().end(), std::mt19937(std::random_device()()));
            std::shuffle(instance->getJobsSmallerOnM2().begin(), instance->getJobsSmallerOnM2().end(), std::mt19937(std::random_device()()));
        }
        auto start = std::chrono::steady_clock::now();
        if (useRevisitedAlgo) revisitedAlgorithm();
        auto endSolve{std::chrono::steady_clock::now()};
        time_elapsed_revisited_johnson = std::chrono::duration<double>{endSolve - start};
        auto cmax1 = evaluate();

        // shuffle list jobs to start from scratch
        std::shuffle(instance->getJobsSmallerOnM1().begin(), instance->getJobsSmallerOnM1().end(), std::mt19937(std::random_device()()));
        std::shuffle(instance->getJobsSmallerOnM2().begin(), instance->getJobsSmallerOnM2().end(), std::mt19937(std::random_device()()));
        start = std::chrono::steady_clock::now();
        JohnsonAlgorithm();
        endSolve = std::chrono::steady_clock::now();
        time_elapsed_johnson = std::chrono::duration<double>{endSolve - start};
        auto cmax3 = evaluate();
        if (useRevisitedAlgo && std::fabs(cmax1-cmax3) > 1E-6) throw F2CmaxException(("Not same Cmax: revisited ->" + std::to_string(cmax1) + " johnson: " + std::to_string(cmax3)).c_str());
        bool conditionProp2 = instance->getSumPa1() <= instance->getSumPa2() - instance->getPMaxA();
        bool conditionProp3 = instance->getSumPb1() <= instance->getSumPb2() - instance->getPMaxB();
        bool conditionProp5 = instance->getSumPa1()+instance->getSumPb2() <= instance->getSumPa2() + instance->getSumPb1() - std::max(instance->getPMaxA(),instance->getPMaxB());
        bool conditionProp6 = instance->getSumPa2() + instance->getSumPb1() <= instance->getSumPa1()+instance->getSumPb2() - std::max(instance->getPMaxA(),instance->getPMaxB());
        if (conditionProp5) {
            auto [k_a,k_a_p] = compute_k_index(A);
            metrics = {5,k_a,k_a_p,0,0,0};
        }else if (conditionProp6) {
            auto [k_b,k_b_p] = compute_k_index(B);
            metrics = {0,0,0,6,k_b,k_b_p};
        }else {
            if (conditionProp2) {
                auto [k_a,k_a_p] = compute_k_index(A);
                std::get<0>(metrics) = 2;
                std::get<1>(metrics) = k_a;
                std::get<2>(metrics) = k_a_p;
            }
            if (conditionProp3) {
                auto [k_b,k_b_p] = compute_k_index(B);
                std::get<3>(metrics) = 2;
                std::get<4>(metrics) = k_b;
                std::get<5>(metrics) = k_b_p;
            }
        }
        objective = cmax3;
    }

    std::pair<size_t,size_t> compute_k_index(SIDE side) {
        auto &listJob = side == A ? instance->getJobsSmallerOnM1() : instance->getJobsSmallerOnM2();

        // identify the smallest index in johnson order
        size_t k = 1;
        for (; k < listJob.size();k++) {
            if (property_2_holds(listJob,0,k,side)) {
                break;
            }
        }
        size_t k_p = k == 0 ? k : k-1;
        // now find the smallest index with different processing time
        while (--k_p>0) {
            if (listJob[k_p].first != listJob[k_p+1].first) break;
        }
        // find the greatest index with different processing time
        while (k < listJob.size()) {
            if (listJob[k].first != listJob[k+1].first) break;
            k++;
        }
        return {k,k_p+1};
    }

    void JohnsonAlgorithm() {
        auto &jobsM1 = instance->getJobsSmallerOnM1();
        std::sort(jobsM1.begin(),jobsM1.end(),[](auto &jobLeft,auto &jobRight){return jobLeft.first < jobRight.first;});
        auto &jobsM2 = instance->getJobsSmallerOnM2();
        std::sort(jobsM2.begin(),jobsM2.end(),[](auto &jobLeft,auto &jobRight){return jobLeft.first < jobRight.first;});
    }

    void revisitedAlgorithm() {
        // Attention, on set B, we work with reverse flo shop instance, i.e. all jobs on machine M1 are in fact on machine M2 and vice versa.
        bool conditionProp2 = instance->getSumPa1() <= instance->getSumPa2() - instance->getPMaxA();
        bool conditionProp3 = instance->getSumPb1() <= instance->getSumPb2() - instance->getPMaxB();
        bool conditionProp5 = instance->getSumPa1()+instance->getSumPb2() <= instance->getSumPa2() + instance->getSumPb1() - std::max(instance->getPMaxA(),instance->getPMaxB());
        bool conditionProp6 = instance->getSumPa2() + instance->getSumPb1() <= instance->getSumPa1()+instance->getSumPb2() - std::max(instance->getPMaxA(),instance->getPMaxB());
        if (conditionProp5) {
            //with version using pivot
            // size_t k_a = find_smallest_k(instance->getJobsSmallerOnM1(),A);
            size_t k_a = find_smallest_k(instance->getJobsSmallerOnM1());
            if (static_cast<double>(k_a) * std::log(static_cast<double>(k_a)) <= static_cast<double>(instance->getNbJobs())) {
                // sort the beginning until the pivot
                std::sort(instance->getJobsSmallerOnM1().begin(), instance->getJobsSmallerOnM1().begin() + k_a);
            }else {
                //sort all the machine
                std::sort(instance->getJobsSmallerOnM1().begin(), instance->getJobsSmallerOnM1().end());
            }
            return;
        }
        if (conditionProp6) {
            //with version using pivot
            // size_t k_b = find_smallest_k(instance->getJobsSmallerOnM2(),B);
            size_t k_b = find_smallest_k(instance->getJobsSmallerOnM2());
            double k_b_p =static_cast<double>( instance->getNbJobs() - k_b + 1);
            if (k_b_p * std::log(k_b_p) <= static_cast<double>(instance->getNbJobs())) {
                // sort the beginning until the pivot
                std::sort(instance->getJobsSmallerOnM2().begin(), instance->getJobsSmallerOnM2().begin() + k_b);
            }else {
                //sort all the machine
                std::sort(instance->getJobsSmallerOnM2().begin(), instance->getJobsSmallerOnM2().end());
            }
            return;
        }
        if (not conditionProp2 && not conditionProp3) {
            //sort both
            std::sort(instance->getJobsSmallerOnM1().begin(), instance->getJobsSmallerOnM1().end());
            std::sort(instance->getJobsSmallerOnM2().begin(), instance->getJobsSmallerOnM2().end());
        }
        size_t k_a = 0;
        size_t k_b = 0;
        if (conditionProp2) {
            //with version using pivot
            // k_a = find_smallest_k(instance->getJobsSmallerOnM1(),A);
            k_a = find_smallest_k(instance->getJobsSmallerOnM1());
            if (static_cast<double>(k_a) * std::log(static_cast<double>(k_a)) <= static_cast<double>(instance->getNbJobs())) {
                // sort the beginning until the pivot
                std::sort(instance->getJobsSmallerOnM1().begin(), instance->getJobsSmallerOnM1().begin() + k_a);
            }else {
                //sort all the machine
                std::sort(instance->getJobsSmallerOnM1().begin(), instance->getJobsSmallerOnM1().end());
            }
        }
        if (conditionProp3) {
            //with version using pivot
            // k_b = find_smallest_k(instance->getJobsSmallerOnM2(),B);
            k_b = find_smallest_k(instance->getJobsSmallerOnM2());
            double k_b_p =static_cast<double>( instance->getNbJobs() - k_b + 1);
            if (k_b_p * std::log(k_b_p) <= static_cast<double>(instance->getNbJobs())) {
                // sort the beginning until the pivot
                std::sort(instance->getJobsSmallerOnM2().begin(), instance->getJobsSmallerOnM2().begin() + k_b);
            }else {
                //sort all the machine
                std::sort(instance->getJobsSmallerOnM2().begin(), instance->getJobsSmallerOnM2().end());
            }
        }
    }

    /**
     * Method that evaluate a solution
     */
    double evaluate() {
        std::vector<double> listCompletionTimeM1(instance->getNbJobs(), 0.0);
        std::vector<double> listCompletionTimeM2(instance->getNbJobs(), 0.0);
        size_t indexLoopCj = 0;
        for (auto &job : instance->getJobsSmallerOnM1()) {
            indexLoopCj == 0 ? listCompletionTimeM1[indexLoopCj] = job.first : listCompletionTimeM1[indexLoopCj] = job.first + listCompletionTimeM1[indexLoopCj-1];
            listCompletionTimeM2[indexLoopCj] = indexLoopCj == 0 ? job.second + listCompletionTimeM1[indexLoopCj] : job.second + std::max(listCompletionTimeM1[indexLoopCj],listCompletionTimeM2[indexLoopCj-1]);
            ++indexLoopCj;
        }
        // loop through job in M2 in reverse order, because we use the reverse property (inverse also the processing time use)
        for (auto &job : std::ranges::reverse_view(instance->getJobsSmallerOnM2())) {
            indexLoopCj == 0 ? listCompletionTimeM1[indexLoopCj] = job.second : listCompletionTimeM1[indexLoopCj] = job.second + listCompletionTimeM1[indexLoopCj-1];
            listCompletionTimeM2[indexLoopCj] = indexLoopCj == 0 ? job.first + listCompletionTimeM1[indexLoopCj] : job.first + std::max(listCompletionTimeM1[indexLoopCj],listCompletionTimeM2[indexLoopCj-1]);
            ++indexLoopCj;
        }
        assert(std::is_sorted(listCompletionTimeM1.begin(),listCompletionTimeM1.end()));
        assert(std::is_sorted(listCompletionTimeM2.begin(),listCompletionTimeM2.end()));
        return listCompletionTimeM2.back();
    }

    size_t find_smallest_k(std::vector<Instance::Job> & listJobs) {
        // take the 10th jobs
        double estimated_pj = std::ceil(std::max(instance->getPMaxA(),instance->getPMaxB()) / instance->getNbJobs() * 10);
        // split the list of jobs in two, those with a pj smaller than the estimate one and the other
        Instance::Job pivot_Job = {estimated_pj,estimated_pj};
        size_t loopIndexBegin = 0;
        size_t loopIndexEnd = listJobs.size()-1;
        while (loopIndexEnd > 0 && loopIndexBegin <= loopIndexEnd) {
            assert(loopIndexBegin < listJobs.size());
            assert(loopIndexEnd < listJobs.size());
            if (listJobs[loopIndexBegin].first < pivot_Job.first) loopIndexBegin++;
            else if (listJobs[loopIndexEnd].first > pivot_Job.first) loopIndexEnd--;
            else {
                std::swap(listJobs[loopIndexBegin],listJobs[loopIndexEnd]);
                loopIndexBegin++;
                loopIndexEnd--;
            }
        }
        return loopIndexEnd + 1; //return the position of the pivot
    }

    size_t find_smallest_k(std::vector<Instance::Job> & listJobs,SIDE side) {
        size_t pivot_t_minus_1 = 0;
        size_t startIndex = 0;
        size_t endIndex = listJobs.size()-1;
        size_t pivot_t = BFPRTPivot(listJobs,startIndex,endIndex);
        while (pivot_t != pivot_t_minus_1) {
            pivot_t_minus_1 = pivot_t;
            if (pivotRule == BFPRT) {
                // check if property 2 is hold
                if (property_2_holds(listJobs,startIndex,pivot_t,side)) {
                    //if pivot is 10 we stop
                    if (pivot_t <= 10) break;
                    endIndex = pivot_t;
                    pivot_t = BFPRTPivot(listJobs,startIndex,endIndex);
                }else {
                    pivot_t = BFPRTPivot(listJobs,pivot_t,endIndex);;
                    endIndex = pivot_t;
                }
            }
        }
        return pivot_t;
    }


    bool property_2_holds(std::vector<Instance::Job> & listJobs,size_t startIndex,size_t endIndex,SIDE side) {
        double sum_diff_pj = std::accumulate(listJobs.begin() + startIndex, listJobs.begin() + endIndex, 0.0, [](double sum,Instance::Job & job) {
            return sum + job.first - job.second;
        });
        return sum_diff_pj <= (side == A ? -instance->getPMaxA() : -instance->getPMaxB());
    }

    static size_t BFPRTPivot(std::vector<Instance::Job> & listJobs,size_t startIndex,size_t endIndex,size_t c=5,size_t d=2) {
        assert(startIndex <= endIndex ); // right indexes
        assert(startIndex < listJobs.size());
        assert(endIndex < listJobs.size());
        size_t pivot = 0;
        // if there is less than c elements, sort them and take the middle
        size_t n = endIndex - startIndex + 1;
        if (n <= c) {
            std::sort(listJobs.begin() + startIndex, listJobs.begin()+endIndex);
            return (startIndex + endIndex) / 2;
        }
        size_t numberJobInT = 0;
        // compute the median of median set
        while (n > c) {
            // create the column with c element, and sort them
            auto itEndColumn = listJobs.begin();
            size_t indexLoopColumn = 0;
            while (itEndColumn != listJobs.end()){
                auto itStartColumn = listJobs.begin() + startIndex + indexLoopColumn * c;
                itEndColumn= (indexLoopColumn + 1)* c < n ? listJobs.begin() + startIndex + (indexLoopColumn + 1)* c : listJobs.end();
                std::sort(itStartColumn, itEndColumn);
                // create set T by putting it at the beginning of the list
                if (indexLoopColumn * c + d < n) std::swap(listJobs[startIndex + indexLoopColumn], listJobs[startIndex + indexLoopColumn * c + d]);
                indexLoopColumn++; // try another column
            }
            // minus one for the index loop column to get the real number
            // The last index L of last added element to T: L=(indexLoopColumn - 1)  * c + d . So the number of element is L - first Element (i.e. 0) + 1 => L + 1
            // pay attention on if there is still enough job to reach the d^th element in last column
            numberJobInT = (indexLoopColumn - 1)  * c + d + 1 < n ? indexLoopColumn : indexLoopColumn >= 1 ? indexLoopColumn -1 : 0;
            endIndex = startIndex + numberJobInT - 1;
            if (numberJobInT == 1) {
                pivot = 0;
                break;
            }
            n = numberJobInT;
        }
        // return the middle of the last element
        std::sort(listJobs.begin() + startIndex, listJobs.begin()+endIndex);
        pivot = (startIndex + endIndex) / 2;

        // split in place the set of jobs
        Instance::Job pivot_Job = listJobs[pivot];
        //mark the pivot with infinity processing time in order to find him in O(n) and swap it to its right position
        listJobs[pivot].first = std::numeric_limits<double>::infinity();
        size_t loopIndexBegin = 0;
        size_t loopIndexEnd = listJobs.size()-1;
        while (loopIndexBegin <= loopIndexEnd) {
            if (listJobs[loopIndexBegin].first < pivot_Job.first) loopIndexBegin++;
            else if (listJobs[loopIndexEnd].first > pivot_Job.first) loopIndexEnd--;
            else {
                std::swap(listJobs[loopIndexBegin],listJobs[loopIndexEnd]);
                loopIndexBegin++;
                loopIndexEnd--;
            }
        }
        // find position of the pivot job
        size_t indexPivotJob = 0;
        for (; indexPivotJob < listJobs.size();indexPivotJob++) if (listJobs[indexPivotJob].first == std::numeric_limits<double>::infinity()) {listJobs[indexPivotJob].first = pivot_Job.first; break;}
        std::swap(listJobs[loopIndexEnd + 1],listJobs[indexPivotJob]);
        return loopIndexEnd + 1; //return the position of the pivot
    }

    /********************/
    /*      GETTER      */
    /********************/

    [[nodiscard]] std::string getPivotRule(){
        std::string pivotName;
        switch (pivotRule) {
        case BFPRT:
            pivotName = "BFPRT";
            break;
        }
        return pivotName;
    }

    /********************/
    /*      SETTER      */
    /********************/

    void setStrategy(std::string pivotName) {
        if (pivotName == "BFPRT") pivotRule = BFPRT;
        else throw F2CmaxException("The pivot rule is not known for the revisited Johnson's algorithm, read \"README\" file for more details on which pivot rule to use.");
    }

    void printOutput(std::string &fileOutputName, std::ofstream &outputFile) {
    bool fileExists = std::filesystem::exists(fileOutputName);
    auto filePath = std::filesystem::path(fileOutputName);
    std::filesystem::create_directories(filePath.lexically_normal().parent_path());
    outputFile.open(fileOutputName, std::ios::out | std::ios::app | std::ios::ate);
    // print header
    if (!fileExists) {
        outputFile <<
            "InstanceName"
            << "\t" << "InstancePath"
            << "\t" << "n"
            << "\t" << "pmax"
            << "\t" << "TimeJohnson";
        if (useRevisitedAlgo) outputFile << "\t" << "TimeRevisitedJohnson";
        outputFile
            << "\t" << "PptA"
            << "\t" << "K_a"
            << "\t" << "K_a_p"
            << "\t" << "PptB"
            << "\t" << "K_b"
            << "\t" << "K_b_p"
            << "\t" << "Objective" << std::endl;
    }
    auto [ppt1,k_a,k_a_p,ppt2,k_b,k_b_p] = metrics;
    // write value
    outputFile << instance->getInstanceName()
               << "\t" << instance->getInstancePath().string()
               << "\t" << instance->getNbJobs()
               << "\t" << instance->getSupPj()
               << "\t" << time_elapsed_johnson.count();
                if (useRevisitedAlgo) outputFile << "\t" << time_elapsed_revisited_johnson.count();
    outputFile
               << "\t" << ppt1
               << "\t" << k_a
               << "\t" << k_a_p
               << "\t" << ppt2
               << "\t" << k_b
               << "\t" << k_b_p
               << "\t" << objective << std::endl;
    outputFile.close();
}
};

#endif //F2_CMAX_SOLVER_H
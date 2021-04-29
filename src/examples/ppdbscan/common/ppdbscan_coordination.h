/**
 \file ppdbscan_coordination.h
 \author moellering@encrypto.cs.tu-darmstadt.de
 \brief
 */

#ifndef ABY_TRACLUS_COORDINATION_H
#define ABY_TRACLUS_COORDINATION_H

#include <iostream>
#include <string>
#include <queue>
#include <abycore/aby/abyparty.h>
#include "../../../abycore/circuit/booleancircuits.h"
#include <abycore/circuit/arithmeticcircuits.h>
#include "../../../abycore/sharing/sharing.h"
#include "../../../abycore/sharing/sharing.h"
#include <stdio.h>
#include <time.h>
#include <memory>
#include "CircuitWrapper.h"

class SharedInputRecord {
private:

public:
  SharedInputRecord();
  share_p clusterId{};
  share_p isNoise{};
  share_p notProcessed{};
  share_p neighborhoodsize{};
  share_p validNeighborhoodSize{};
  share_p simdResultDistancesCompToThresholds{};
  std::vector <share_p> resultDistancesCompToThresholds;
};

class VerificationLineSegment{
private:

public:
  VerificationLineSegment();
  uint64_t clusterId{};
  uint64_t id{};
  uint64_t isNoise{};
  uint64_t notProcessed{};
  uint64_t neighborhoodsize{};
  uint64_t validNeighborhoodSize{};
  uint64_t checkNecessary{};
  std::vector <uint64_t > resultDistancesCompToThreshold;
};

int32_t testDistancesCalculation(e_role role, const std::string& address, uint16_t port, seclvl seclvl, uint32_t bitlen, uint32_t nthreads, e_mt_gen_alg mt_alg, e_sharing distance_sharing, uint32_t noLineSegments, std::string dataset);
int32_t testClustering(e_role role, const std::string& address, uint16_t port, seclvl seclvl,
                       uint32_t bitlen, uint32_t nthreads, e_mt_gen_alg mt_alg, e_sharing distance_sharing, uint32_t noLineSegments, std::string dataset, int32_t i);

std::vector<std::vector<uint64_t>> loadData(const std::string& fileName);
std::vector<std::vector<uint64_t>> generateServer2DummyInput(int noDummys);
std::vector<std::vector<uint64_t > >transposeVector(std::vector<std::vector<uint64_t >> &b);

std::vector<share_p>
calculateDistances(uint32_t bitlen, uint32_t noLineSegments, uint64_t dim,
                   uint64_t threshold, CircuitW_p &boolcirc, CircuitW_p &arithmcirc,
                   CircuitW_p &yaocirc,
                   std::vector<std::vector<uint64_t>> &lineSegments,e_role &role, ABYParty* party);
void buildSquaredEuclideanDistCircuit(std::vector<share_p> &shrDb, std::vector<share_p> &shrCompElement, uint64_t noElements, uint64_t dim,
                                       CircuitW_p &boolcirc,CircuitW_p &arithmcirc, CircuitW_p &yaocirc, uint64_t thresh, uint64_t bitlen, e_role &role, ABYParty* party);
void buildApproxTraclusDistanceCircuit(std::vector<share_p> &shrDb, std::vector<share_p> &shrCompElement, uint64_t noElements, uint64_t dim,
                                                        CircuitW_p &boolcirc,CircuitW_p &arithmcirc, CircuitW_p &yaocirc, uint64_t thresh, uint64_t bitlen, e_role &role, ABYParty* party);
std::tuple<share_p,std::vector<SharedInputRecord>>
createLineSegments(uint32_t noLineSegments, uint64_t minLns, CircuitW_p &yaocirc,
                   std::vector<share_p> &simdDistances);
void buildClusteringCircuit(int64_t noElements, CircuitW_p &arithmcirc, CircuitW_p &yaocirc, CircuitW_p &boolcirc, e_role bitlen, ABYParty* party, int32_t i);
void buildSquaredEuclideanDistCircuit(std::vector<share_p> &shrDb, std::vector<share_p> &shrCompElement, uint64_t noElements, uint64_t dim,
                                       CircuitW_p &boolcirc,CircuitW_p &arithmcirc, CircuitW_p &yaocirc, uint64_t thresh, uint64_t bitlen, e_role &role, ABYParty* party);
void endClustering(std::vector<SharedInputRecord> &sharedLineSegments,int64_t noLineSegments, CircuitW_p &arithmcirc, CircuitW_p &yaocirc, CircuitW_p &boolcirc, ABYParty* party);
void putOutGates(std::vector<SharedInputRecord>& lineSegments, uint64_t noElements, CircuitW_p &arithmcirc, CircuitW_p &yaocirc, CircuitW_p &boolcirc);


std::vector<uint64_t > combineDbValuesForSIMD(std::vector<std::vector<uint64_t > > &lineSegments, uint64_t i, uint64_t j, uint64_t noLineSegments);
share_p combineSharesToSIMD(std::vector<share_p *> &non_simd_shares,std::size_t bitlen, CircuitW_p &circ);
share_p combineSharesToSIMD(std::vector<share_p> &non_simd_shares,std::size_t bitlen, CircuitW_p &circ);
uint32_t OrTree(std::vector<uint32_t> wires,uint32_t length, CircuitW_p &circ);


void prepareLineSegmentsForFile(std::vector<SharedInputRecord> &lineSegments, share_p &currentClusterId, uint32_t noLineSegments, CircuitW_p &boolcirc,CircuitW_p &arithmcirc);
void writeSimdDistancesToFile(share_p simdDistances, e_role role);
void writeLineSegmentsToFile(std::vector<SharedInputRecord> &lineSegments, share_p &currentClusterId, CircuitW_p &boolcirc, e_role role);
void writeLineSegmentsToFile(std::vector<SharedInputRecord> &lineSegments, share_p &currentClusterId, CircuitW_p &boolcirc, e_role role);
std::vector<share_p> getSimdDistancesFromFile(CircuitW_p &boolcirc,CircuitW_p &yaocirc,e_role role,uint32_t noElements);
std::tuple<share_p, std::vector<SharedInputRecord>> getLineSegmentsFromFile(CircuitW_p &boolcirc, CircuitW_p &arithmcirc,CircuitW_p &yaocirc,e_role role,uint32_t noElements);


void verifyMinApproxTraclusDistanceWithShares(std::vector< std::vector<uint64_t > > server1shares, std::vector< std::vector<uint64_t > > server2shares, uint64_t dbsize);
std::vector<std::vector<u_int64_t >> verifyApproxTraclusDistance(std::vector< std::vector<uint64_t > > lineSegments, uint64_t dbsize, uint64_t dim,uint64_t threshold);
void verifyClustering(std::vector<VerificationLineSegment> lineSegments, uint64_t dbsize);
void verifyClusteringAsInAby(std::vector<VerificationLineSegment> lineSegments, uint64_t dbsize);
std::vector<std::vector<u_int64_t >> verifySquaredEuclideanDistance(std::vector< std::vector<uint64_t > > lineSegments, uint64_t dbsize, uint64_t dim,uint64_t threshold);

#endif //ABY_TRACLUS_COORDINATION_H

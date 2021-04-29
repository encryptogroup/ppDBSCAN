/**
 \file ppdbscan_coordination.cpp
 \author moellering@encrypto.cs.tu-darmstadt.de
 \brief
 */

#include "ppdbscan_coordination.h"
#include <iomanip>

SharedInputRecord::SharedInputRecord() = default;

VerificationLineSegment::VerificationLineSegment() = default;

uint64_t MAX_CLUSTERS = UINT8_MAX;
uint64_t dim = 2;
uint64_t minPts = 3;
uint64_t threshold = 200000000000;


/**
 * \brief Called for calculating distances in a shared manner with ABY
 * \param role server or client
 * \param bitlen shares' bitlength
 * \param noLineSegments number of elements that shall be clustered
 * \param dataset file name that contains the data set that shall be clustered
 */
int32_t testDistancesCalculation(e_role role, const std::string& address, uint16_t port, seclvl seclvl, uint32_t bitlen, uint32_t nthreads, e_mt_gen_alg mt_alg, e_sharing distance_sharing, uint32_t noLineSegments, std::string dataset) {
  auto* party = new ABYParty(role, address, port, seclvl, bitlen, nthreads,
                             mt_alg);

  //get vector with all supported sharing in the framework
  std::vector<Sharing*>& sharings = party->GetSharings();

  //get circuit for sharing type
  CircuitW_p arithmcirc = std::make_shared<CircuitWrapper>(sharings[distance_sharing]->GetCircuitBuildRoutine());
  CircuitW_p yaocirc = std::make_shared<CircuitWrapper>(sharings[S_YAO]->GetCircuitBuildRoutine());
  CircuitW_p boolcircuit = std::make_shared<CircuitWrapper>(sharings[S_BOOL]->GetCircuitBuildRoutine());

  //Will also be created for server 1 (role 0) to verify the result
  std::vector<std::vector<uint64_t>> lineSegments =
      generateServer2DummyInput(noLineSegments);
  if (role == SERVER) {
    lineSegments = loadData("data/"+dataset+".txt");
  }

  //set shared input
  std::vector<share_p>  simdDistances =
      calculateDistances(bitlen, noLineSegments, dim, threshold, boolcircuit,arithmcirc,
                         yaocirc,  lineSegments,role,party);

  //split simd to create line segment objects
  auto[initialClusterId, sharedLineSegments] = createLineSegments(noLineSegments, minPts, yaocirc, simdDistances);

  //free memory
  simdDistances.clear();
  simdDistances.shrink_to_fit();
  lineSegments.clear();
  lineSegments.shrink_to_fit();

  //prepare data to write it to file
  prepareLineSegmentsForFile(sharedLineSegments,initialClusterId,noLineSegments,boolcircuit,arithmcirc);
  party->ExecCircuit();
  //write data to file
  writeLineSegmentsToFile(sharedLineSegments,initialClusterId,boolcircuit,role);
  //clear memory
  party->Reset();
  sharedLineSegments.clear();
  sharedLineSegments.shrink_to_fit();

  delete party;
  return 0;
}

/**
 * \brief Called for executing DBSCAN clustering on shared data stored at the two clients in files
 * \param role server or client
 * \param bitlen shares' bitlength
 * \param noLineSegments number of elements that shall be clustered
 * \param dataset file name that contains the data set that shall be clustered
 */
int32_t testClustering(e_role role, const std::string& address, uint16_t port, seclvl seclvl, uint32_t bitlen, uint32_t nthreads, e_mt_gen_alg mt_alg, e_sharing distance_sharing, uint32_t noLineSegments, std::string dataset,int32_t i) {

    auto* party = new ABYParty(role, address, port, seclvl, bitlen, nthreads,
                                   mt_alg);

    //get vector with all supported sharing in the framework
    std::vector<Sharing*>& sharings = party->GetSharings();

    //get circuit for sharing type
    CircuitW_p arithmcirc = std::make_shared<CircuitWrapper>(sharings[distance_sharing]->GetCircuitBuildRoutine());
    CircuitW_p yaocirc = std::make_shared<CircuitWrapper>(sharings[S_YAO]->GetCircuitBuildRoutine());
    CircuitW_p boolcircuit = std::make_shared<CircuitWrapper>(sharings[S_BOOL]->GetCircuitBuildRoutine());

    //Will also be created for server 1 (role 0) to verify the result
    std::vector<std::vector<uint64_t>> lineSegments =
        generateServer2DummyInput(noLineSegments);
    if (role == SERVER && i==0) {
        lineSegments = loadData("data/"+dataset+".txt");

        std::vector<std::vector<u_int64_t >> testThresholdComparison = verifyApproxTraclusDistance(lineSegments,noLineSegments,dim,threshold);//verifySquaredEuclideanDistance(lineSegments,noLineSegments,dim,threshold);

        std::vector<VerificationLineSegment> verificationLineSegments;
        for(uint32_t u=0;u<noLineSegments;u++){
          for(uint32_t t=0;t<noLineSegments;t++){
            if(u==0) {
              verificationLineSegments.emplace_back();
              verificationLineSegments[t].clusterId=MAX_CLUSTERS;
              verificationLineSegments[t].id=t;
              verificationLineSegments[t].isNoise = 0;
              verificationLineSegments[t].notProcessed = 1;
              verificationLineSegments[t].neighborhoodsize = 0;
            }
            verificationLineSegments[t].resultDistancesCompToThreshold.push_back(
                    testThresholdComparison[t][u]);
          }
        }

        for(uint32_t u=0;u<noLineSegments;u++){
          for(uint32_t t=0;t<noLineSegments;t++){//count neighbors
            verificationLineSegments[u].neighborhoodsize += verificationLineSegments[u].resultDistancesCompToThreshold[t];
          }
          if (verificationLineSegments[u].neighborhoodsize> minPts){
            verificationLineSegments[u].validNeighborhoodSize = 1;
          }else{
            verificationLineSegments[u].validNeighborhoodSize = 0;
          }
        }

        verifyClustering(verificationLineSegments, noLineSegments);
        verifyClusteringAsInAby(verificationLineSegments,noLineSegments);
    }

    //actual 2PC clustering
    buildClusteringCircuit(noLineSegments,arithmcirc,yaocirc,boolcircuit,role, party,i);
  return 0;
}

/**
 * Reads partitioned and quantized data from file
 * @param fileName
 * @return data set
 */
std::vector<std::vector<uint64_t>> loadData(const std::string& fileName){
  std::ifstream file;
  file.open(fileName);
  std::string line;

  std::vector<std::vector<uint64_t>> lineSegments;
  if (file.is_open()) {
    while (std::getline(file, line)) {
      std::string word;
      std::vector<uint64_t> substrList;
      int lineSize = line.size();
      for (int i=0;i< lineSize;i++) {
        if (line[i] != ';') {
          word += line[i];
        } else {
          if ((uint64_t) word.size() != 0) {
            substrList.push_back(stoi(word));
          }
          word = "";
        }
      }
      lineSegments.push_back(substrList);
    }
  }
  file.close();
  return lineSegments;
}

/**
 * Generate dummy input of 2nd server: all 0s
 * \param noDummys number of dummy elements that have to be created
 * @return 2-dim vector with 0s
 */
std::vector<std::vector<uint64_t>> generateServer2DummyInput(int noDummys){
  std::vector<std::vector<uint64_t>> lineSegments;
  for(int i=0;i<noDummys;i++){
    std::vector<uint64_t> segments;
    for(int j=0; j<4;j++) {
      segments.push_back(0);
    }
    lineSegments.push_back(segments);
  }
  return lineSegments;
}

/**
 * \brief Calculates pairwise distances (so far either approx. Traclus or SED) in a shared manner between input elements and determines if the distance is smaller than a threshold
 * @param bitlen shares' bitlength
 * @param noLineSegments number of elements of which the distances should be calculated
 * @param dim number of parameters per element
 * @param threshold threshold for distances
 * @param boolcirc boolean circuit
 * @param arithmcirc arithm circuit
 * @param yaocirc Garbled circuit
 * @param lineSegments input elements
 * @param role server or client
 * @param party ABY party
 * @return vector with SIMD shares indicating if the pairwise distances where smaller than the threshold
 */
std::vector<share_p>
calculateDistances(uint32_t bitlen, uint32_t noLineSegments, uint64_t dim,
                   uint64_t threshold, CircuitW_p &boolcirc, CircuitW_p &arithmcirc,
                   CircuitW_p &yaocirc,
                   std::vector<std::vector<uint64_t>> &lineSegments,e_role &role, ABYParty* party) {
  //remove old files
  remove("simdDistanceSharesStorageClient.txt");
  remove("simdDistanceSharesStorageServer.txt");

  //distance calculation
  std::vector<std::vector<uint64_t > > transposed_vector =
      transposeVector(lineSegments);
  for (uint32_t i=0;i<(noLineSegments); i++) {
    std::vector<share_p> shrCompElements;
    std::vector<share_p> newEntry;
    for (uint64_t j = 0; j < dim; j++) {
      newEntry.push_back(arithmcirc->PutSharedSIMDINGate(noLineSegments,combineDbValuesForSIMD(lineSegments,i,j,noLineSegments).data(), bitlen));
      //if(i==0) {
      shrCompElements.push_back(arithmcirc->PutSharedSIMDINGate(noLineSegments, transposed_vector[j].data(), bitlen));
      //}
    }

    //TODO: add variable for TRACLUS or SED
    buildSquaredEuclideanDistCircuit(newEntry, shrCompElements, noLineSegments,
                                     dim, boolcirc, arithmcirc, yaocirc,
                                     threshold, bitlen, role, party);
    //buildApproxTraclusDistanceCircuit(newEntry, shrCompElements, noLineSegments,
    //                                 dim, boolcirc, arithmcirc, yaocirc,
    //                                 threshold, bitlen, role, party);
    std::cout<<i<<std::endl;
  }

  //read full results from file
  return getSimdDistancesFromFile(boolcirc,yaocirc,role,noLineSegments);
}

/**
 * \brief calculates SED in a parallel manner
 * @param shrDb input elements
 * @param shrCompElement input elements
 * @param noElements number of input elements
 * @param dim number of parameters per element
 * @param boolcirc
 * @param arithmcirc
 * @param yaocirc
 * @param thresh
 * @param bitlen shares' bitlen
 * @param role server or client
 * @param party ABY party
 */
void buildSquaredEuclideanDistCircuit(std::vector<share_p> &shrDb, std::vector<share_p> &shrCompElement, uint64_t noElements, uint64_t dim,
                                          CircuitW_p &boolcirc,CircuitW_p &arithmcirc, CircuitW_p &yaocirc, uint64_t thresh, uint64_t bitlen, e_role &role, ABYParty* party){

  share_p temp;
  uint32_t j;
  uint64_t zero = 0;

  share_p compResult;
  share_p threshold = yaocirc->PutSIMDCONSGate(noElements, thresh, bitlen);
  share_p sqed = arithmcirc->PutSIMDCONSGate(noElements, zero, bitlen);
  for (j = 0; j < dim; j++) {
    temp = arithmcirc->PutSUBGate(shrDb[j], shrCompElement[j]);
    temp = arithmcirc->PutMULGate(temp, temp);
    sqed = arithmcirc->PutADDGate(temp, sqed);
  }

  // comparison to threshold
  sqed = yaocirc->PutA2YGate(sqed);
  compResult = yaocirc->PutGTGate(threshold, sqed);


  //prepare data for storing to file
  compResult = boolcirc->PutY2BGate(compResult);
  compResult = boolcirc->PutSharedOUTGate(compResult);
  party->ExecCircuit();

  //write Result to File
  writeSimdDistancesToFile(compResult,role);

  //clear memory
  party->Reset();
}

/**
 * \brief Builds a circuit for the approximated Traclus distance metric that combines 4 squared euclidean distances.
 * \param shrDb: elements to which shrCompElement will be compared by calculating euclidean distance (every element has dim parameters)
 * \param shrCompElement: element to wish all other elements will be compared, has dim parameters
 * \param dim: number of parameters per line segment, normally 4
 */
void buildApproxTraclusDistanceCircuit(std::vector<share_p> &shrDb, std::vector<share_p> &shrCompElement, uint64_t noElements, uint64_t dim,
                                                        CircuitW_p &boolcirc,CircuitW_p &arithmcirc, CircuitW_p &yaocirc, uint64_t thresh, uint64_t bitlen, e_role &role, ABYParty* party){
  share_p temp;
  uint32_t j;
  uint64_t zero = 0;

  share_p compResult;
  share_p threshold = yaocirc->PutSIMDCONSGate(noElements, thresh, bitlen);
  share_p sqed = arithmcirc->PutSIMDCONSGate(noElements, zero, bitlen);
  for (j = 0; j < dim; j++) {
    temp = arithmcirc->PutSUBGate(shrDb[j], shrCompElement[j]);
    temp = arithmcirc->PutMULGate(temp, temp);
    sqed = arithmcirc->PutADDGate(temp, sqed);
    temp = arithmcirc->PutSUBGate(shrCompElement[j],shrDb[(j+2)%dim]);//cross difference
    temp = arithmcirc->PutMULGate(temp, temp);
    sqed = arithmcirc->PutADDGate(temp,sqed);
  }

  // comparison to threshold
  sqed = yaocirc->PutA2YGate(sqed);
  compResult = yaocirc->PutGTGate(threshold, sqed);


  //prepare data for storing to file
  compResult = boolcirc->PutY2BGate(compResult);
  compResult = boolcirc->PutSharedOUTGate(compResult);
  party->ExecCircuit();

  //write Result to File
  writeSimdDistancesToFile(compResult,role);

  //clear memory
  party->Reset();
}

/**
 * \brief creates objects for each input element that later on used for clustering
 * @param noLineSegments number of objects
 * @param minPts minimum number of elements required for creating a cluster
 * @param yaocirc
 * @param simdDistances comparison result of the distances with the treshold
 * @return tuple containing the shared initialized clusterID (1) and a vector with shared lineSegements
 */
std::tuple<share_p,std::vector<SharedInputRecord>>
createLineSegments(uint32_t noLineSegments, uint64_t minPts, CircuitW_p &yaocirc,
                   std::vector<share_p> &simdDistances) {
  uint64_t zero = 0;
  uint64_t one = 1;
  std::vector<SharedInputRecord> sharedLineSegments;
  share_p minPtsShare = yaocirc->PutCONSGate(minPts,16);
  share_p sZero = yaocirc->PutCONSGate(zero,1);
  share_p sOne = yaocirc->PutCONSGate(one,1);
  share_p maxCluster = yaocirc->PutCONSGate(MAX_CLUSTERS,8);
  std::vector<std::vector<uint32_t>> wireIdMatrix;
  for(uint32_t i=0;i<noLineSegments;i++){//iterate through all pairwise distances
    for(uint32_t j=0;j<noLineSegments;j++) {//iterate through all line segments
      if(i==0){
        wireIdMatrix.push_back(std::vector<uint32_t >(noLineSegments,0));
        sharedLineSegments.emplace_back();
        sharedLineSegments[j].clusterId = maxCluster;
        sharedLineSegments[j].isNoise = sZero;
        sharedLineSegments[j].notProcessed = sOne;
      }
      uint32_t posids[1] = {j};
      share_p tmp = yaocirc->PutSubsetGate(simdDistances[i],posids,1);
      sharedLineSegments[j].resultDistancesCompToThresholds.push_back(tmp);
      wireIdMatrix[j][i] = tmp->get_wire_id(zero);

      if(i==noLineSegments-1){//last iteration done
        sharedLineSegments[j].neighborhoodsize = yaocirc->PutHammingWeightGateRec(wireIdMatrix[j].data(),noLineSegments);//TODO:remove this attribute?
        //yaocirc->PutPrintValueGate(sharedLineSegments[j].neighborhoodsize,"neighborhoodsize");
        sharedLineSegments[j].validNeighborhoodSize = yaocirc->PutGTGate(sharedLineSegments[j].neighborhoodsize, minPtsShare);
        sharedLineSegments[j].simdResultDistancesCompToThresholds = combineSharesToSIMD(sharedLineSegments[j].resultDistancesCompToThresholds,1,yaocirc);
      }
    }
  }
  share_p sCurrentClusterId = yaocirc->PutCONSGate(one,8);
  return {sCurrentClusterId,sharedLineSegments};
}

/**
 * \brief clusters inputs elements in a privacy-preserving way (DBSCAN)
 * @param noElements number of elements
 * @param arithmcirc
 * @param yaocirc
 * @param boolcirc
 * @param role server or client
 * @param party ABY party
 * @param i: number of iteration
 */
void buildClusteringCircuit(int64_t noElements, CircuitW_p &arithmcirc, CircuitW_p &yaocirc, CircuitW_p &boolcirc, e_role role, ABYParty* party, int32_t i){

  uint64_t zero = 0;
  uint64_t one = 1;

  std::cout <<i <<std::endl;
  share_p sZero = yaocirc->PutCONSGate(zero,1);
  share_p sOne = yaocirc->PutCONSGate(one,1);

  share_p simdZero = yaocirc->PutSIMDCONSGate(noElements,zero,1);
  share_p simdOne = yaocirc->PutSIMDCONSGate(noElements,one,1);
  auto[sCurrentClusterId,lineSegments] = getLineSegmentsFromFile(boolcirc,arithmcirc,yaocirc,role,noElements);

  //prepare SIMD
  std::vector<share_p *> collectionValidNeighborHoodSize;
  for(uint32_t j=0;j< noElements;j++){
    collectionValidNeighborHoodSize.push_back(&lineSegments[j].validNeighborhoodSize);
  }
  share_p simdValidNeighborhoodSize=combineSharesToSIMD(collectionValidNeighborHoodSize,1,yaocirc);

  SharedInputRecord & elem = lineSegments[i];
  share_p simdCurrentClusterId = yaocirc->PutRepeaterGate(noElements,sCurrentClusterId);

  // check if line segment is not classified & size of neigborhood is greater than minPts, set cluster index
  share_p temp= yaocirc->PutMUXGate(sCurrentClusterId, sZero, elem.validNeighborhoodSize);
  elem.clusterId = yaocirc->PutMUXGate(temp,elem.clusterId,elem.notProcessed);
  share_p sComp = yaocirc -> PutANDGate(elem.validNeighborhoodSize,elem.notProcessed);
  elem.isNoise = yaocirc->PutMUXGate(yaocirc->PutMUXGate(sZero,sOne,elem.validNeighborhoodSize),elem.isNoise,elem.notProcessed);
  elem.notProcessed = sZero;

  //update simdIsNoise, simdNotProcessed, simdClusterIds TODO:implement replacement method for values in SIMD gates
  std::vector<share_p *> newClusterIds;
  std::vector<share_p *> newIsNoise;
  std::vector<share_p *> newNotProcessed;
  for(SharedInputRecord &ls:lineSegments){
    newClusterIds.push_back(&ls.clusterId);
    newIsNoise.push_back(&ls.isNoise);
    newNotProcessed.push_back(&ls.notProcessed);
  }
  share_p simdClusterIdsOfElements = combineSharesToSIMD(newClusterIds,8,yaocirc);
  share_p simdIsNoise = combineSharesToSIMD(newIsNoise,1,yaocirc);
  share_p simdNotProcessed = combineSharesToSIMD(newNotProcessed,1,yaocirc);

  int32_t numIterations = 0;
  share_p simdComp = yaocirc->PutRepeaterGate(noElements,sComp);
  while (numIterations!=3){//iterate through all other line segments to obliviously check neighbors;
    share_p temp2 = yaocirc->PutANDGate(yaocirc->PutANDGate(yaocirc->PutORGate(simdIsNoise, simdNotProcessed),
                                                            simdComp), elem.simdResultDistancesCompToThresholds);
    simdClusterIdsOfElements =
        yaocirc->PutMUXGate(simdCurrentClusterId, simdClusterIdsOfElements, temp2);
    simdNotProcessed =
        yaocirc->PutMUXGate(simdZero, simdNotProcessed, temp2);
    simdIsNoise= yaocirc->PutMUXGate(simdZero, simdIsNoise, temp2);
    share_p simdComp2 =
        yaocirc->PutMUXGate(simdValidNeighborhoodSize, simdZero, temp2);

    for(uint32_t k=0;k<((int32_t)noElements);k++) {
      // check if other elements are neighbor of elem2 (if elem2 is valid/comp2)
      share_p simdUpdateResultDistancesCompToThresholds = yaocirc->PutMUXGate(simdOne, simdZero,yaocirc->PutANDGate(simdComp2,lineSegments[k].simdResultDistancesCompToThresholds));//TODO: update simdCollectedDistancesVector-> necessary?
      share_p splittedDistanceUpdate = yaocirc->PutSplitterGate(simdUpdateResultDistancesCompToThresholds);

      //if one of the valid neighbors had k as neighbor (simdUpdateResultDistancesCompToThresholds contains one 1), set elem.resultDistancesCompToThresholds[k] to 1
      std::vector<uint32_t> updateWires = splittedDistanceUpdate->get_wires();
      elem.resultDistancesCompToThresholds[k] = yaocirc->PutMUXGate(sOne,elem.resultDistancesCompToThresholds[k],share_p(new boolshare(std::vector<uint32_t> {OrTree(updateWires,noElements,yaocirc)},((BooleanCircuit*)yaocirc.get()))));
    }
    numIterations +=1;

    //update isNoise & notProcessed & simdResultDistancesCompToThresholds of elements
    std::vector<uint32_t > splittedIsNoiseUpdatesWires = yaocirc->PutSplitterGate(simdIsNoise)->get_wires();
    for(uint32_t l=0;l<noElements;l++)lineSegments[l].isNoise=share_p(new boolshare(std::vector<uint32_t > {splittedIsNoiseUpdatesWires[l]},((BooleanCircuit*)yaocirc.get())));

    std::vector<uint32_t > splittedNotProcessedUpdatesWires = yaocirc->PutSplitterGate(simdNotProcessed)->get_wires();;
    for(uint32_t l=0;l<noElements;l++) lineSegments[l].notProcessed=share_p(new boolshare(std::vector<uint32_t > {splittedNotProcessedUpdatesWires[l]},((BooleanCircuit*)yaocirc.get())));

    std::vector<share_p *> newsimdResultDistancesCompToThresholds;
    for(uint32_t l=0;l<noElements;l++) {
      newsimdResultDistancesCompToThresholds.push_back(&elem.resultDistancesCompToThresholds[l]);
    }
    elem.simdResultDistancesCompToThresholds=combineSharesToSIMD(newsimdResultDistancesCompToThresholds,1,yaocirc);
  }
  //update clusterIds of elements
  newClusterIds.clear();
  for(uint32_t l=0;l<noElements;l++) {
    share_p splittedClusterIdUpdates = yaocirc->PutSubsetGate(simdClusterIdsOfElements,{&l},1);
    lineSegments[l].clusterId =splittedClusterIdUpdates;
    newClusterIds.push_back(&lineSegments[l].clusterId);
  }
  if(i<noElements-1){
    simdClusterIdsOfElements=combineSharesToSIMD(newClusterIds,8,yaocirc);//TODO:replace this with more efficient way (can I replace just 1 bit in the wire?)

    sCurrentClusterId = yaocirc->PutADDGate(sCurrentClusterId,sComp);

    // prepare to intermediate results
    prepareLineSegmentsForFile(lineSegments, sCurrentClusterId, noElements,
                               boolcirc,arithmcirc);
    // execute circuit
    party->ExecCircuit();
    // write intermediate results to file
    writeLineSegmentsToFile(lineSegments, sCurrentClusterId, boolcirc, role);
    // reset
    party->Reset();
  }else{//last element checked
    endClustering(lineSegments,noElements,arithmcirc,yaocirc,boolcirc,party);
  }
}

/**
 * finalizes the clustering
 * @param sharedLineSegments secret-shared objects of the input records
 * @param noLineSegments is the data set size
 * @param arithmcirc
 * @param yaocirc
 * @param boolcirc
 * @param party is server or client
 */
void endClustering(std::vector<SharedInputRecord> &sharedLineSegments,int64_t noLineSegments, CircuitW_p &arithmcirc, CircuitW_p &yaocirc, CircuitW_p &boolcirc, ABYParty* party){
  putOutGates(sharedLineSegments,noLineSegments, arithmcirc,yaocirc,boolcirc);
  party->ExecCircuit();

  std::vector<uint64_t> results(30,0);
  for(uint64_t j = 0; j < noLineSegments; j++) {
    SharedInputRecord ls = sharedLineSegments[j];

    auto output1 = ls.clusterId->get_clear_value<uint64_t>();
    std::cout << "Cluster ID: " << output1<< std::endl;
    if(output1<noLineSegments){
      results[output1]+=1;
    }
    auto output4 = ls.notProcessed->get_clear_value<uint64_t>();
    std::cout << "Not  processed: " << output4<< std::endl;

    auto output5 = ls.isNoise->get_clear_value<uint64_t>();
    std::cout << "Marked  as noise: " << output5 << std::endl;

    auto output3 = ls.neighborhoodsize->get_clear_value<uint64_t>();
    std::cout << "Neighborhood Size: " << output3<< std::endl;

    auto output2 = ls.validNeighborhoodSize->get_clear_value<uint64_t>();
    std::cout << "Valid Neighborhood Size: " << output2<< std::endl;
    for(uint64_t i=0; i< noLineSegments;i++){
      auto output6 = ls.resultDistancesCompToThresholds[i]->get_clear_value<uint64_t>();
      std::cout << j << ". element compared with " << (i) << ". element: " << output6<< std::endl;
      //std::cout << output6<< std::endl;
    }

    uint32_t* output;
    uint32_t out_bitlen, out_nvals;
  }

  for(uint32_t j=0;j< 30;j++){
    if (j==0){
      std::cout << results[j] << " elements were marked as noise." << std::endl;
    }else{
      if (results[j]!=0) {
        std::cout << "Cluster " << j << " contains " << results[j]
                  << " line segments." << std::endl;
      }
    }
  }


  party->Reset();
  delete party;
}

/**
 * adds output gates to all secret-shared objects of input records
 * @param lineSegments are the secret-shared input recrods
 * @param noElements is the data set size
 * @param arithmcirc
 * @param yaocirc
 * @param boolcirc
 */
void putOutGates(std::vector<SharedInputRecord>& lineSegments, uint64_t noElements, CircuitW_p &arithmcirc, CircuitW_p &yaocirc, CircuitW_p &boolcirc){
  for(uint32_t i=0;i<noElements;i++){
    SharedInputRecord & elem = lineSegments[i];
    elem.clusterId = arithmcirc->PutOUTGate(arithmcirc->PutY2AGate(elem.clusterId,boolcirc->circ_),ALL);

    for(uint32_t j=0;j<noElements;j++){
      elem.resultDistancesCompToThresholds[j]=yaocirc->PutOUTGate(elem.resultDistancesCompToThresholds[j],ALL);
    }

    elem.validNeighborhoodSize = arithmcirc->PutOUTGate(arithmcirc->PutY2AGate(elem.validNeighborhoodSize,boolcirc->circ_),ALL);
    elem.neighborhoodsize = yaocirc->PutOUTGate(elem.neighborhoodsize,ALL);
    elem.isNoise = arithmcirc->PutOUTGate(arithmcirc->PutY2AGate(elem.isNoise,boolcirc->circ_),ALL);
    elem.notProcessed = arithmcirc->PutOUTGate(arithmcirc->PutY2AGate(elem.notProcessed,boolcirc->circ_),ALL);
    elem.simdResultDistancesCompToThresholds = yaocirc->PutOUTGate(elem.simdResultDistancesCompToThresholds,ALL);
  }
}

/**
 * Combines shared input elements to a SIMD share
 * @param lineSegments are the single shared input elements
 * @param i
 * @param j
 * @param noLineSegments is the data set size
 * @return combined SIMD share
 */
std::vector<uint64_t > combineDbValuesForSIMD(std::vector<std::vector<uint64_t>> &lineSegments, uint64_t i, uint64_t j, uint64_t noLineSegments){//TODO:same as transpose now?
  std::vector<uint64_t> combinedDbValues;

  for(uint32_t k=0; k<noLineSegments;k++){
    combinedDbValues.push_back(lineSegments[i][j]);
  }

  return combinedDbValues;
}

/**
 * transposes a vector
 * @param b is the vector that shall be transposed
 * @return the transposed vector
 */
std::vector<std::vector<uint64_t >> transposeVector(std::vector<std::vector<uint64_t > > &b)
{

    std::vector<std::vector<uint64_t> > trans_vec(b[0].size(), std::vector<uint64_t>());

    for (auto & i : b)
    {
        for (uint32_t j = 0; j < i.size(); j++)
        {
            trans_vec[j].push_back(i[j]);
        }
    }

    return trans_vec;
}

/**
 * puts Output gates for the objects of ShareInputRecord, s.t. the secret shares can be written into a file at each computing party. This is needed to free memory (& solve the memory problems of ABY) to cluster more data.
 * @param lineSegments are the secret-shared input data records
 * @param currentClusterId is the secret-shared current cluster id
 * @param noLineSegments is the data set size
 * @param boolcirc
 * @param arithmcirc
 */
void prepareLineSegmentsForFile(std::vector<SharedInputRecord> &lineSegments, share_p &currentClusterId, uint32_t noLineSegments, CircuitW_p &boolcirc, CircuitW_p &arithmcirc){
  //convert yao shares to GMW (needed for storing shares on drive)
  currentClusterId = boolcirc->PutSharedOUTGate(boolcirc->PutY2BGate(currentClusterId));//boolcirc->PutSharedOUTGate(boolcirc->PutA2BGate(currentClusterId,yaocirc->circ));
  for(SharedInputRecord &temp:lineSegments){
    temp.clusterId=boolcirc->PutSharedOUTGate(boolcirc->PutY2BGate(temp.clusterId));
    temp.isNoise=boolcirc->PutSharedOUTGate(boolcirc->PutY2BGate(temp.isNoise));
    temp.notProcessed=boolcirc->PutSharedOUTGate(boolcirc->PutY2BGate(temp.notProcessed));
    //boolcirc->PutPrintValueGate(temp.neighborhoodsize,"distances1");
    temp.neighborhoodsize=boolcirc->PutSharedOUTGate(boolcirc->PutY2BGate(temp.neighborhoodsize));
    //boolcirc->PutPrintValueGate(temp.neighborhoodsize,"distances2");
    temp.validNeighborhoodSize=boolcirc->PutSharedOUTGate(boolcirc->PutY2BGate(temp.validNeighborhoodSize));
    temp.simdResultDistancesCompToThresholds=boolcirc->PutSharedOUTGate(boolcirc->PutY2BGate(temp.simdResultDistancesCompToThresholds));
    for(uint32_t j=0;j<noLineSegments;j++){
      temp.resultDistancesCompToThresholds[j]=(boolcirc->PutSharedOUTGate(boolcirc->PutY2BGate(temp.resultDistancesCompToThresholds[j])));
    }
  }
}

/**
 * writes the secret-shared simd distances into a file. (also needed for freeing memory)
 * @param simdDistances are the shares that shall be written into a file
 * @param role server or client
 */
void writeSimdDistancesToFile(share_p simdDistances, e_role role){
  std::string roleString= "Client";
  if(role==SERVER) roleString="Server";
  std::ofstream serverStorage("simdDistanceSharesStorage"+roleString+".txt", std::ios::app);

  uint32_t *output, bitlen, nvals;
  simdDistances->get_clear_value_vec(&output,&bitlen,&nvals);
  for(uint32_t j=0;j<nvals;j++){
    serverStorage << output[j]<<"\n";
  }
  serverStorage.close();
}

/**
 * reads secret-shared simd distances from a file.
 * @param boolcirc
 * @param yaocirc
 * @param role server or client
 * @param noElements data set size
 * @return SIMD secret shares
 */
std::vector<share_p> getSimdDistancesFromFile(CircuitW_p &boolcirc,CircuitW_p &yaocirc,e_role role,uint32_t noElements){
  std::string roleString= "Client";
  if(role==SERVER) roleString="Server";
  std::ifstream serverStorage("simdDistanceSharesStorage"+roleString+".txt");

  std::vector<share_p> results;
  results.reserve(noElements);
  for(uint32_t i=0;i<noElements;i++){//all elements
    std::vector<uint8_t > newSimdCollection;
    for(uint32_t j=0;j<noElements;j++) {//all distances
      uint8_t temp;
      serverStorage >> temp;
      newSimdCollection.push_back(temp);
      }
    results.push_back(yaocirc->PutB2YGate(boolcirc->PutSharedSIMDINGate(noElements,newSimdCollection.data(),1)));
  }
  serverStorage.close();
  return results;
}

/**
 * writes secret-shared obeject of input data records into a file
 * @param lineSegments are the secret-shared input records
 * @param currentClusterId is the secret-shared current cluster ID
 * @param boolcirc
 * @param role is server or client
 */
void writeLineSegmentsToFile(std::vector<SharedInputRecord> &lineSegments, share_p &currentClusterId, CircuitW_p &boolcirc, e_role role){
  //convert yao shares to GMW (needed for storing shares on drive)
  std::ofstream serverStorage;
  std::string roleString= "Client";
  if(role==SERVER) roleString="Server";
  serverStorage.open ("lineSegmentsStorage"+roleString+".txt");

  serverStorage << currentClusterId->get_clear_value<uint64_t >() << "\n";

  for(SharedInputRecord &temp :lineSegments){
    serverStorage << temp.clusterId->get_clear_value<uint64_t>()<<"\n";
    serverStorage << temp.isNoise->get_clear_value<uint64_t>()<<"\n";
    serverStorage << temp.notProcessed->get_clear_value<uint64_t>()<<"\n";
    serverStorage << temp.neighborhoodsize->get_clear_value<uint64_t>()<<"\n";
    serverStorage << temp.validNeighborhoodSize->get_clear_value<uint64_t>()<<"\n";
    uint32_t *output, bitlen, nvals;
    temp.simdResultDistancesCompToThresholds->get_clear_value_vec(&output,&bitlen,&nvals);
    for(uint32_t j=0;j<nvals;j++){
      serverStorage << output[j]<<"\n";
    }
  }
  serverStorage.close();
}

/**
 * reads the secret-shares of the input data records from the file
 * @param boolcirc
 * @param arithmcirc
 * @param yaocirc
 * @param role is server or client
 * @param noElements is the dataset size
 * @return
 */
std::tuple<share_p, std::vector<SharedInputRecord>> getLineSegmentsFromFile(CircuitW_p &boolcirc, CircuitW_p &arithmcirc, CircuitW_p &yaocirc,e_role role,uint32_t noElements) {
  std::string roleString= "Client";
  if(role==SERVER) roleString="Server";
  std::ifstream serverStorage("lineSegmentsStorage"+roleString+".txt");

  //read current number of clusters from file
  uint64_t clId;
  serverStorage >> clId;
  share_p currentClusterId = yaocirc->PutB2YGate(boolcirc->PutSharedINGate(clId,8));

  std::vector<SharedInputRecord> lineSegments;
  lineSegments.reserve(noElements);
  for(uint32_t i=0;i<noElements;i++){
    lineSegments.emplace_back();
    std::vector<share_p> collectionResultDistancesComp;
    for(uint32_t j=0;j<(noElements+5);j++) {
      uint16_t temp;
      serverStorage >> temp;
      switch(j) {
        case 0:lineSegments[i].clusterId =
            yaocirc->PutB2YGate(boolcirc->PutSharedINGate(temp, 8));
            break;
        case 1:lineSegments[i].isNoise =
                 yaocirc->PutB2YGate(boolcirc->PutSharedINGate(temp, 1));
        break;
        case 2:lineSegments[i].notProcessed =yaocirc->PutB2YGate(boolcirc->PutSharedINGate(temp, 1));
        break;
        case 3:lineSegments[i].neighborhoodsize =
                 yaocirc->PutB2YGate(boolcirc->PutSharedINGate(temp, 8));//TODO:remove?
        break;
        case 4:lineSegments[i].validNeighborhoodSize =yaocirc->PutB2YGate(
                     boolcirc->PutSharedINGate(temp, 1));
        break;
        default:
          share_p dist = yaocirc->PutB2YGate(boolcirc->PutSharedINGate(temp, 1));
          lineSegments[i].resultDistancesCompToThresholds.push_back(dist);
          break;
      }
    }
    lineSegments[i].simdResultDistancesCompToThresholds=combineSharesToSIMD(lineSegments[i].resultDistancesCompToThresholds,1,yaocirc);
  }
  serverStorage.close();
  return {currentClusterId,lineSegments};
}


/**
 * combines single secret shares to SIMD secret shares
 * @param non_simd_shares is the input single secret shares
 * @param bitlen is the bit length
 * @param circ
 * @return returns the created SIMD shares
 */
share_p combineSharesToSIMD(std::vector<share_p *> &non_simd_shares,std::size_t bitlen, CircuitW_p &circ){
  using wire_t = std::uint32_t;
  // result buffer
  std::vector<wire_t> v_wires_simd;

  for(std::size_t i = 0; i < bitlen; ++i){
    std::vector<wire_t> wire_i;
    wire_i.reserve(bitlen);
    for(auto s: non_simd_shares) wire_i.emplace_back((*s)->get_wires()[i]);
    v_wires_simd.emplace_back(circ->PutCombinerGate(wire_i));
  }

  return share_p(new boolshare(v_wires_simd,((BooleanCircuit*)circ.get())));
}

/**
 * combines single secret shares to SIMD secret shares
 * @param non_simd_shares is the input single secret shares
 * @param bitlen is the bit length
 * @param circ
 * @return
 */
share_p combineSharesToSIMD(std::vector<share_p> &non_simd_shares,std::size_t bitlen, CircuitW_p &circ){
  using wire_t = std::uint32_t;
  // result buffer
  std::vector<wire_t> v_wires_simd;

  for(std::size_t i = 0; i < bitlen; ++i){
    std::vector<wire_t> wire_i;
    wire_i.reserve(bitlen);
    for(auto &s: non_simd_shares) wire_i.emplace_back(s->get_wires()[i]);
    v_wires_simd.emplace_back(circ->PutCombinerGate(wire_i));
  }

  return share_p(new boolshare(v_wires_simd,((BooleanCircuit*)circ.get())));
}

/**
 * Implements an OR-tree
 * @param wires is the input data
 * @param length is the number of objects "to OR"
 * @param circ
 * @return
 */
uint32_t OrTree(std::vector<uint32_t> wires,uint32_t length, CircuitW_p &circ) {//TODO:vector by reference?
  uint32_t result;
  float f_length = float(length);
  if (length>1) {
    for (uint32_t i = 0; i < length; i++) {
      // left
      uint32_t left = OrTree(std::vector<uint32_t> (wires.begin(), wires.begin() + (length/2)), floor(f_length/2),circ);
      // right
      uint32_t right = OrTree(std::vector<uint32_t>(wires.begin() + length/2, wires.end()), ceil(f_length/2),circ);
      return (((BooleanCircuit*)circ->circ_)->PutORGate(left,right));
    }
  }
  return wires[0];
}

/**
 * plaintext implementation of approx. TRACLUS distance measure
 * @param server1shares are the plaintext shares of party 1
 * @param server2shares are the plaintext shares of party 2
 * @param dbsize is the data base size
 */
void verifyMinApproxTraclusDistanceWithShares(std::vector< std::vector<uint64_t > > server1shares, std::vector< std::vector<uint64_t > > server2shares, uint64_t dbsize) {
    uint32_t i, j;
    uint64_t tmpdist = 0;

    std::vector<uint64_t > mindist;
    for(i=0; i < dbsize; i++) {
      mindist.push_back(MAX_CLUSTERS);

      std::vector<uint64_t > share_i {server1shares[i][0] + server2shares[i][1],server1shares[i][1] + server2shares[i][2],server1shares[i][2] + server2shares[i][2],server1shares[i][3] + server2shares[i][3]};
      std::vector<uint64_t > share_j(4);
      j = (i+1)%dbsize;
      while(j!=i){
        share_j = {server1shares[j][0] + server2shares[j][1],server1shares[j][1] + server2shares[j][2],server1shares[j][2] + server2shares[j][2],server1shares[j][3] + server2shares[j][3]};

        tmpdist = 2*(share_i[0]*share_i[0] + share_i[1]*share_i[1]+share_i[2]*share_i[2]+share_i[3]*share_i[3]+share_j[0]*share_j[0]+share_j[1]*share_j[1]+share_j[2]*share_j[2]+share_j[3]*share_j[3]);
        tmpdist -= 2*(share_i[0]*share_j[0]+share_i[2]*share_j[2]+share_i[0]*share_j[2]+share_i[2]*share_j[0]+share_i[1]*share_j[1]+share_i[3]*share_j[3]+share_i[1]*share_j[3]+share_i[3]*share_j[1]);


        if(tmpdist < mindist[i]) {
          mindist[i] = tmpdist;
        }
        j = (j+1)%dbsize;
      }
    }
}

/**
 * \brief squared euclidean distance, implementation has the same structure as clustering with ABY
 * \param lineSegments: elements
 * \param dbsize: number of line segments
 * \param dim: number of parameters per line segment
 * \param threshold: epsilon value (max. distance between neighbors)
 */
std::vector<std::vector<u_int64_t >> verifySquaredEuclideanDistance(std::vector< std::vector<uint64_t > > lineSegments, uint64_t dbsize, uint64_t dim,uint64_t threshold){
    std::vector<std::vector<uint64_t >> results;
    for(uint32_t i = 0; i<dbsize; i++){
        std::vector<uint64_t> distances;
        std::vector<uint64_t> elementA = lineSegments[i];

      for(uint32_t j = 0; j<dbsize; j++){
        std::vector<uint64_t> elementB = lineSegments[j];
        int64_t temp = 0;
        for(uint32_t k=0; k<dim;k++){
          temp += (elementA[k]- elementB[k])*(elementA[k]- elementB[k]);
        }
        distances.push_back(temp);
      }
      std::vector<u_int64_t > tmp;
      for(uint32_t l=0;l<dbsize;l++){
        if(threshold>distances[l]){
          tmp.push_back(1);
        }else{
          tmp.push_back(0);
        }
      }
      results.push_back(tmp);
    }
    return results;
}

/**
 * \brief approx. Traclus distance, implementation has the same structure as clustering with ABY
 * \param lineSegments: elements
 * \param dbsize: number of line segments
 * \param dim: number of parameters per line segment
 * \param threshold: epsilon value (max. distance between neighbors)
 */
std::vector<std::vector<u_int64_t >> verifyApproxTraclusDistance(std::vector< std::vector<uint64_t > > lineSegments, uint64_t dbsize, uint64_t dim,uint64_t threshold){
  std::vector<std::vector<uint64_t >> results;
  for(uint32_t i = 0; i<dbsize; i++){
    std::vector<uint64_t> distances;
    std::vector<uint64_t> elementA = lineSegments[i];

    for(uint32_t j = 0; j<dbsize; j++){
      std::vector<uint64_t> elementB = lineSegments[j];
      int64_t temp = 0;
      for(uint32_t k=0; k<dim;k++){
        temp += (elementA[k]- elementB[k])*(elementA[k]- elementB[k]);
        temp += (elementA[k]- elementB[((k+2)%dim)%dim])*(elementA[k]- elementB[(k+2)%dim]);
      }
      distances.push_back(temp);
    }
    std::vector<u_int64_t > tmp;
    for(uint32_t l=0;l<dbsize;l++){
      if(threshold>distances[l]){
        tmp.push_back(1);
      }else{
        tmp.push_back(0);
      }
    }
    results.push_back(tmp);
  }
  return results;
}


/**
 * \brief DBScan clustering, implementation with queue
 * \param lineSegments: elements
 * \param dbsize: number of line segments
 */
void verifyClustering(std::vector<VerificationLineSegment> lineSegments, uint64_t dbsize){
  uint64_t currentClusterId = 1;
  for(uint32_t i=0;i<dbsize;i++) {
    VerificationLineSegment &elem = lineSegments[i];
    if (elem.notProcessed) { // no cluster
      if (elem.validNeighborhoodSize) { // neighborhood big enough
        elem.clusterId = currentClusterId;
        elem.notProcessed = 0;

        std::queue<VerificationLineSegment *> queue;
        for (uint32_t j = 0; j < dbsize; j++) {
          if (elem.resultDistancesCompToThreshold[j] == 1) {
            queue.push(&lineSegments[j]);
          }
        }

        while (not queue.empty()) {
          VerificationLineSegment *elem2 = queue.front();

          if ((*elem2).isNoise) { // noise
            (*elem2).clusterId = currentClusterId;
            (*elem2).isNoise = 0;
          }

          if ((*elem2).notProcessed) { // unclustered
            if((*elem2).id>=200){
              uint32_t test = 1;
            }
            (*elem2).clusterId = currentClusterId;
            (*elem2).notProcessed = 0;

            if ((*elem2).validNeighborhoodSize) {
              for (uint32_t k = 0; k < dbsize; k++) {
                if ((*elem2).resultDistancesCompToThreshold[k] == 1) {
                  queue.push(
                      &lineSegments[k]);
                }
              }
            }
          }
          queue.pop();
        }
        currentClusterId += 1;
      } else {
        elem.clusterId = 0;
        elem.isNoise = 1;
        elem.notProcessed = 0;// mark as noise
      }
    }
  }

  for(uint32_t i=0; i<dbsize; i++){
    VerificationLineSegment ls = lineSegments[i];

    std::cout << "Cluster ID: " << ls.clusterId<< std::endl;
    std::cout << "Not  processed: " << ls.notProcessed<< std::endl;
    std::cout << "Marked  as noise: " << ls.isNoise << std::endl;
    std::cout << "Neighborhood Size: " << ls.neighborhoodsize<< std::endl;
    std::cout << "Valid Neighborhood Size: " << ls.validNeighborhoodSize<< std::endl;
    for(uint32_t j=0; j< dbsize;j++){
      std::cout << i << ". element compared with " << (j) << ". element: " << ls.resultDistancesCompToThreshold[j]<< std::endl;
    }
  }

  std::vector<uint64_t> results(currentClusterId,0);
  for (const auto& segm:lineSegments){
    results[segm.clusterId]+=1;
  }
  for(uint32_t j=0;j< currentClusterId;j++){
    if (j==0){
      std::cout << results[j] << " elements were marked as noise." << std::endl;
    }else{
      std::cout << "Cluster " << j << " contains " << results[j] << " line segments." << std::endl;
    }
  }
}

/**
 * plaintext computation of ppDBSCAN
 * @param lineSegments is the input data set
 * @param dbsize is the data set size
 */
void verifyClusteringAsInAby(std::vector<VerificationLineSegment> lineSegments, uint64_t dbsize){
  uint64_t currentClusterId = 1;
  for(int32_t i=0;i<dbsize;i++) {
    VerificationLineSegment &elem = lineSegments[i];
    if (elem.notProcessed) { // no cluster assigned
      if (elem.validNeighborhoodSize) { // neighborhood big enough
        elem.clusterId = currentClusterId;
        elem.notProcessed = 0;

        int32_t numIterations = 0;
        while(numIterations!=3){//queue replacement
          for (int32_t j=0;j<((int32_t)dbsize);j++) {
            VerificationLineSegment *elem2 = &lineSegments[j];
            if (elem.resultDistancesCompToThreshold[j] == 1) {
              if ((*elem2).isNoise) { // noise
                (*elem2).clusterId = currentClusterId;
                (*elem2).isNoise = 0;
              }
              if ((*elem2).notProcessed) { // unclustered
                (*elem2).clusterId = currentClusterId;
                (*elem2).notProcessed = 0;
                if ((*elem2).validNeighborhoodSize) {
                  for(int32_t k=0;k<dbsize;k++) {
                    if ((*elem2).resultDistancesCompToThreshold[k] == 1) {
                      elem.resultDistancesCompToThreshold
                      [k] =1;
                    }
                  }
                }
              }
            }
          }
          numIterations+=1;
        }
        currentClusterId += 1;
      } else {
        elem.clusterId = 0;// mark as noise
        elem.isNoise = 1;
        elem.notProcessed = 0;
      }
    }
  }


  for(uint32_t i=0; i<dbsize; i++){
    VerificationLineSegment ls = lineSegments[i];

    std::cout << "Cluster ID: " << ls.clusterId<< std::endl;
    std::cout << "Not  processed: " << ls.notProcessed<< std::endl;
    std::cout << "Marked  as noise: " << ls.isNoise << std::endl;
    std::cout << "Neighborhood Size: " << ls.neighborhoodsize<< std::endl;
    std::cout << "Valid Neighborhood Size: " << ls.validNeighborhoodSize<< std::endl;
    for(uint32_t j=0; j< dbsize;j++){
      std::cout << i << ". element compared with " << (j) << ". element: " << ls.resultDistancesCompToThreshold[j]<< std::endl;
    }
  }

  std::vector<uint64_t> results(currentClusterId,0);
  for (const auto& segm:lineSegments){
    if (currentClusterId>segm.clusterId) { results[segm.clusterId] += 1; }
  }
  for(uint32_t j=0;j< currentClusterId;j++){
    if (j==0){
      std::cout << results[j] << " elements were marked as noise." << std::endl;
    }else{
      std::cout << "Cluster " << j << " contains " << results[j] << " line segments." << std::endl;
    }
  }
}

/**
 * write clustering result into a new file
 * @param clusterId
 * @param lineSegment
 * @param dim
 */
void writeClusteringResultToNewFile(uint64_t clusterId, std::vector<uint64_t> lineSegment, uint64_t dim){
  std::ofstream myfile;
  time_t rawtime;
  struct tm * timeinfo;
  char buffer[80];
  time (&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(buffer,80,"%d%m%Y_%H%M%S",timeinfo);
  myfile.open ("output/clustering_result_"+std::string(buffer)+".txt",std::ios_base::app);
  myfile << clusterId << "; ";
  for(uint64_t i=0;i<dim;i++){
    myfile << lineSegment[i] <<"; ";
  }
  myfile << std::endl;

  myfile.close();

}
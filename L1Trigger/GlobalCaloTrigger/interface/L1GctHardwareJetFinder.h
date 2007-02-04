#ifndef L1GCTHARDWAREJETFINDER_H_
#define L1GCTHARDWAREJETFINDER_H_

#include "L1Trigger/GlobalCaloTrigger/interface/L1GctJetFinderBase.h"

#include <boost/cstdint.hpp> //for uint16_t
#include <vector>

/*! \class L1GctHardwareJetFinder
 * \brief Emulation of the hardware jet finder.
 *
 *  
 */
/*
 * \author Greg Heath
 * \date June 2006
 */



class L1GctHardwareJetFinder : public L1GctJetFinderBase
{
 public:

  /// id is 0-8 for -ve Eta jetfinders, 9-17 for +ve Eta, for increasing Phi.
  L1GctHardwareJetFinder(int id, std::vector<L1GctSourceCard*> sourceCards,
			 L1GctJetEtCalibrationLut* jetEtCalLut);
                 
  ~L1GctHardwareJetFinder();
   
  /// Overload << operator
  friend std::ostream& operator << (std::ostream& os, const L1GctHardwareJetFinder& algo);

  /// include additional reset functionality
  virtual void reset();

  /// get input data from sources
  virtual void fetchInput();

  /// process the data, fill output buffers
  virtual void process();

 protected:

  // Each jetFinder must define the constants as private and copy the
  // function definitions below.
  virtual unsigned maxRegionsIn() const { return MAX_REGIONS_IN; }
  virtual unsigned centralCol0() const { return CENTRAL_COL0; }
  virtual int nCols() const { return N_COLS; }

private:

  /// The real jetFinders must define these constants
  static const unsigned int MAX_REGIONS_IN; ///< Dependent on number of rows and columns.
  static const int N_COLS;
  static const unsigned int CENTRAL_COL0;

  /// Threshold value to find a local et maximum
  static const unsigned int JET_THRESHOLD;

  /// Vector of pre-clustered jets formatted as a search array for the final clustering
  RegionsVector m_protoJetRegions;

  /// Local vectors used during both stages of clustering
  RegionsVector m_localMaxima;
  /// Each local maximum becomes a cluster
  RegionsVector m_clusters;

  /// The number of local Maxima/clusters found at each stage of clustering
  unsigned m_numberOfClusters;
  
  /// The first stage of clustering, called by fetchInput()
  void findProtoJets();  
  /// The second stage of clustering, called by process()
  void findJets();  

  /// Find local maxima in the search array
  void findLocalMaxima();
  /// Convert local maxima to clusters
  void findProtoClusters();
  /// Convert protojets to final jets
  void findFinalClusters();

  // Comparison operator for sorting local maxima
  // In the case where the two et values are equal, we favour
  // the more central region
  struct etGreaterThan : public std::binary_function<L1CaloRegion, L1CaloRegion, bool> 
  {
    bool operator()(const L1CaloRegion& x, const L1CaloRegion& y) {
      return ( (x.et() > y.et()) || ((x.et() == y.et()) && (x.rctEta() < y.rctEta())) ) ;
    }
  };

  /// Fill search array for the second stage of clustering based on the pre-clustered jets
  void fillRegionsFromProtoJets();
  /// Organise the pre-clustered jets into the ones we keep and those we send to the neighbour
  void convertClustersToProtoJets();
  /// Organise the final clustered jets into L1GctJets
  void convertClustersToOutputJets();

};

std::ostream& operator << (std::ostream& os, const L1GctHardwareJetFinder& algo);

#endif /*L1GCTHARDWAREJETFINDER_H_*/

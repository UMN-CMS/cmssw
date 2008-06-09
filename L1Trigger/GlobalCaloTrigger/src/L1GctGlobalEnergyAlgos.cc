#include "L1Trigger/GlobalCaloTrigger/interface/L1GctGlobalEnergyAlgos.h"

#include "L1Trigger/GlobalCaloTrigger/interface/L1GctWheelEnergyFpga.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctWheelJetFpga.h"

#include "FWCore/Utilities/interface/Exception.h"

#include <cassert>

using std::ostream;
using std::endl;
using std::vector;
using std::max;

const unsigned int L1GctGlobalEnergyAlgos::N_JET_COUNTERS_USED=L1GctWheelJetFpga::N_JET_COUNTERS;
const unsigned int L1GctGlobalEnergyAlgos::N_JET_COUNTERS_MAX =12;

L1GctGlobalEnergyAlgos::L1GctGlobalEnergyAlgos(vector<L1GctWheelEnergyFpga*> wheelFpga,
					       vector<L1GctWheelJetFpga*> wheelJetFpga) :
  m_plusWheelFpga(wheelFpga.at(1)),
  m_minusWheelFpga(wheelFpga.at(0)),
  m_plusWheelJetFpga(wheelJetFpga.at(1)),
  m_minusWheelJetFpga(wheelJetFpga.at(0)),
  m_metComponents(),
  m_jcValPlusWheel(N_JET_COUNTERS_USED),
  m_jcVlMinusWheel(N_JET_COUNTERS_USED),
  m_outputJetCounts(N_JET_COUNTERS_MAX)
{
  if(wheelFpga.size() != 2)
  {
    throw cms::Exception("L1GctSetupError")
    << "L1GctGlobalEnergyAlgos::L1GctGlobalEnergyAlgos() : Global Energy Algos has been incorrectly constructed!\n"
    << "This class needs two wheel card pointers. "
    << "Number of wheel card pointers present is " << wheelFpga.size() << ".\n";
  }
  
  if(wheelJetFpga.size() != 2)
  {
    throw cms::Exception("L1GctSetupError")
    << "L1GctGlobalEnergyAlgos::L1GctGlobalEnergyAlgos() : Global Energy Algos has been incorrectly constructed!\n"
    << "This class needs two wheel jet fpga pointers. "
    << "Number of wheel jet fpga pointers present is " << wheelJetFpga.size() << ".\n";
  }
  
    if(m_plusWheelFpga == 0)
    {
      throw cms::Exception("L1GctSetupError")
      << "L1GctGlobalEnergyAlgos::L1GctGlobalEnergyAlgos() has been incorrectly constructed!\n"
      << "Plus Wheel Fpga pointer has not been set!\n";
    }
    if(m_minusWheelFpga == 0)
    {
      throw cms::Exception("L1GctSetupError")
      << "L1GctGlobalEnergyAlgos::L1GctGlobalEnergyAlgos() has been incorrectly constructed!\n"
      << "Minus Wheel Fpga pointer has not been set!\n";
    }
    if(m_plusWheelJetFpga == 0)
    {
      throw cms::Exception("L1GctSetupError")
      << "L1GctGlobalEnergyAlgos::L1GctGlobalEnergyAlgos() has been incorrectly constructed!\n"
      << "Plus Wheel Jet Fpga pointer has not been set!\n";
    }
    if(m_minusWheelJetFpga == 0)
    {
      throw cms::Exception("L1GctSetupError")
      << "L1GctGlobalEnergyAlgos::L1GctGlobalEnergyAlgos() has been incorrectly constructed!\n"
      << "Minus Wheel Jet Fpga pointer has not been set!\n";
    }
}

L1GctGlobalEnergyAlgos::~L1GctGlobalEnergyAlgos()
{
}

ostream& operator << (ostream& os, const L1GctGlobalEnergyAlgos& fpga)
{
  os << "===L1GctGlobalEnergyAlgos===" << endl;
  os << "WheelEnergyFpga* minus = " << fpga.m_minusWheelFpga << endl;
  os << "WheelEnergyFpga* plus  = " << fpga.m_plusWheelFpga << endl;
  os << "WheelJetFpga* minus = " << fpga.m_minusWheelJetFpga << endl;
  os << "WheelJetFpga* plus  = " << fpga.m_plusWheelJetFpga << endl;
  os << "Inputs from Plus wheel:" << endl;
  os << "  Ex " << fpga.m_exValPlusWheel << "\n  Ey " << fpga.m_eyValPlusWheel << endl;
  os << "  Et " << fpga.m_etValPlusWheel << "\n  Ht " << fpga.m_htValPlusWheel << endl; 
  os << "Inputs from Minus wheel:" << endl;
  os << "  Ex " << fpga.m_exVlMinusWheel << "\n  Ey " << fpga.m_eyVlMinusWheel << endl;
  os << "  Et " << fpga.m_etVlMinusWheel << "\n  Ht " << fpga.m_htVlMinusWheel << endl; 
  os << "Input Jet counts " << endl;
  for(unsigned i=0; i < fpga.m_jcValPlusWheel.size(); i++)
    {
      os << "  Plus wheel  " << i << ": " << fpga.m_jcValPlusWheel.at(i);
      os << "  Minus wheel " << i << ": " << fpga.m_jcVlMinusWheel.at(i) << endl;
    } 
  os << endl;
  os << "Output Etmiss " << fpga.m_outputEtMiss << endl;
  os << "Output Etmiss Phi " << fpga.m_outputEtMissPhi << endl;
  os << "Output EtSum " << fpga.m_outputEtSum << endl;
  os << "Output EtHad " << fpga.m_outputEtHad << endl;
  os << "Output Jet counts " << endl;
  for(unsigned i=0; i < fpga.m_outputJetCounts.size(); i++)
    {
      os << i << ": " << fpga.m_outputJetCounts.at(i) << endl;
    } 
  os << endl;

  return os;
}

// clear internal data
void L1GctGlobalEnergyAlgos::reset()
{
  m_exValPlusWheel.reset();
  m_exVlMinusWheel.reset();
  m_eyValPlusWheel.reset();
  m_eyVlMinusWheel.reset();
  m_etValPlusWheel.reset();
  m_etVlMinusWheel.reset();
  m_htValPlusWheel.reset();
  m_htVlMinusWheel.reset();
  for (unsigned i=0; i<N_JET_COUNTERS_USED; i++) {
    m_jcValPlusWheel.at(i).reset();
    m_jcVlMinusWheel.at(i).reset();
  }
  //
  m_outputEtMiss.reset();
  m_outputEtMissPhi.reset();
  m_outputEtSum.reset();
  m_outputEtHad.reset();
  for (unsigned i=0; i<N_JET_COUNTERS_MAX; i++) {
    m_outputJetCounts.at(i).reset();
  }
}

void L1GctGlobalEnergyAlgos::fetchInput() {
  // input from WheelEnergyFpgas
  m_exValPlusWheel = m_plusWheelFpga->getOutputEx();
  m_eyValPlusWheel = m_plusWheelFpga->getOutputEy();
  m_etValPlusWheel = m_plusWheelFpga->getOutputEt();
  m_htValPlusWheel = m_plusWheelJetFpga->getOutputHt();
  
  m_exVlMinusWheel = m_minusWheelFpga->getOutputEx();
  m_eyVlMinusWheel = m_minusWheelFpga->getOutputEy();
  m_etVlMinusWheel = m_minusWheelFpga->getOutputEt();
  m_htVlMinusWheel = m_minusWheelJetFpga->getOutputHt();

  //
  for (unsigned i=0; i<N_JET_COUNTERS_USED; i++) {
    m_jcValPlusWheel.at(i) = m_plusWheelJetFpga->getOutputJc(i);
    m_jcVlMinusWheel.at(i) = m_minusWheelJetFpga->getOutputJc(i);
  }
}


// process the event
void L1GctGlobalEnergyAlgos::process()
{
  etComponentType ExSum, EySum;
  etmiss_vec EtMissing;
  //
  //-----------------------------------------------------------------------------
  // Form the Ex and Ey sums
  ExSum = m_exValPlusWheel + m_exVlMinusWheel;
  EySum = m_eyValPlusWheel + m_eyVlMinusWheel;
  // Execute the missing Et algorithm
  m_metComponents.setComponents(-ExSum, -EySum);
  EtMissing = m_metComponents.metVector();

  //
  m_outputEtMiss    = EtMissing.mag;
  m_outputEtMissPhi = EtMissing.phi;

  //
  //-----------------------------------------------------------------------------
  // Form the Et and Ht sums
  m_outputEtSum = m_etValPlusWheel + m_etVlMinusWheel;
  m_outputEtHad = m_htValPlusWheel + m_htVlMinusWheel;

  //
  //-----------------------------------------------------------------------------
  // Add the jet counts.
  for (unsigned i=0; i<N_JET_COUNTERS_USED; i++) {
    m_outputJetCounts.at(i) =
      L1GctJetCount<5>(m_jcValPlusWheel.at(i)) +
      L1GctJetCount<5>(m_jcVlMinusWheel.at(i));
  }
  // BUT ... overwrite some of the jet counts with Hf tower sums!
  packHfTowerSumsIntoJetCountBits();
}

std::vector<unsigned> L1GctGlobalEnergyAlgos::getJetCountValues() const {
  std::vector<unsigned> jetCountValues(N_JET_COUNTERS_MAX);
  for (unsigned jc=0; jc<N_JET_COUNTERS_MAX; jc++) {
    jetCountValues.at(jc) = m_outputJetCounts.at(jc).value();
  }
  return jetCountValues;
}
  
//----------------------------------------------------------------------------------------------
// set input data per wheel: x component of missing Et
//
void L1GctGlobalEnergyAlgos::setInputWheelEx(unsigned wheel, int energy, bool overflow)
{
  if (wheel==0) {
    m_exValPlusWheel.setValue(energy);
    m_exValPlusWheel.setOverFlow(overflow);
  } else if (wheel==1) {
    m_exVlMinusWheel.setValue(energy);
    m_exVlMinusWheel.setOverFlow(overflow);
  }
}

//----------------------------------------------------------------------------------------------
// set input data per wheel: y component of missing Et
//
void L1GctGlobalEnergyAlgos::setInputWheelEy(unsigned wheel, int energy, bool overflow)
{
  if (wheel==0) {
    m_eyValPlusWheel.setValue(energy);
    m_eyValPlusWheel.setOverFlow(overflow);
  } else if (wheel==1) {
    m_eyVlMinusWheel.setValue(energy);
    m_eyVlMinusWheel.setOverFlow(overflow);
  }
}

//----------------------------------------------------------------------------------------------
// set input data per wheel: scalar sum of Et
//
void L1GctGlobalEnergyAlgos::setInputWheelEt(unsigned wheel, unsigned energy, bool overflow)
{
  if (wheel==0) {
    m_etValPlusWheel.setValue(energy);
    m_etValPlusWheel.setOverFlow(overflow);
  } else if (wheel==1) {
    m_etVlMinusWheel.setValue(energy);
    m_etVlMinusWheel.setOverFlow(overflow);
  }
}

//----------------------------------------------------------------------------------------------
// set input data per wheel: sum of transverse energy in jets (Ht)
//
void L1GctGlobalEnergyAlgos::setInputWheelHt(unsigned wheel, unsigned energy, bool overflow)
{
  if (wheel==0) {
    m_htValPlusWheel.setValue(energy);
    m_htValPlusWheel.setOverFlow(overflow);
  } else if (wheel==1) {
    m_htVlMinusWheel.setValue(energy);
    m_htVlMinusWheel.setOverFlow(overflow);
  }
}


//----------------------------------------------------------------------------------------------
// Set the jet count input values
//
void L1GctGlobalEnergyAlgos::setInputWheelJc(unsigned wheel, unsigned jcnum, unsigned count)
{
  if (jcnum<N_JET_COUNTERS_USED) {
    if (wheel==0) {
      m_jcValPlusWheel.at(jcnum).setValue(count);
    } else if (wheel==1) {
      m_jcVlMinusWheel.at(jcnum).setValue(count);
    }
  }
}

//----------------------------------------------------------------------------------
//
// The following method contains the bit-slicing
// of the Hf information for minbias triggers
// into 5-bit fields so it can be output using the 
// L1GctJetCounts object (from DataFormats/L1GlobalCaloTrigger).
//
// NOTE: the reverse operation of combining the 5-bit fields
// back into Hf information is carried out in L1GctJetCounts.
//
void L1GctGlobalEnergyAlgos::packHfTowerSumsIntoJetCountBits()
{
  assert (N_JET_COUNTERS_USED<=6 && N_JET_COUNTERS_MAX>11);
  m_outputJetCounts.at(6)  = m_plusWheelJetFpga->getOutputHfSums().nOverThreshold;
  m_outputJetCounts.at(7)  = m_minusWheelJetFpga->getOutputHfSums().nOverThreshold;

  m_outputJetCounts.at(8)  = m_plusWheelJetFpga->getOutputHfSums().etSum0;
  m_outputJetCounts.at(9)  = m_minusWheelJetFpga->getOutputHfSums().etSum0;

  m_outputJetCounts.at(10) = m_plusWheelJetFpga->getOutputHfSums().etSum1;
  m_outputJetCounts.at(11) = m_minusWheelJetFpga->getOutputHfSums().etSum1;

}



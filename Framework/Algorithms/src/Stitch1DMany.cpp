#include "MantidAlgorithms/Stitch1DMany.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidAPI/WorkspaceProperty.h"
#include "MantidKernel/ArrayProperty.h"
#include "MantidKernel/BoundedValidator.h"
#include "MantidKernel/RebinParamsValidator.h"

#include <boost/make_shared.hpp>

using namespace Mantid::Kernel;
using namespace Mantid::API;

namespace Mantid {
namespace Algorithms {
DECLARE_ALGORITHM(Stitch1DMany)

/** Initialize the algorithm's properties.
 */
void Stitch1DMany::init() {

  declareProperty(
      make_unique<ArrayProperty<std::string>>("InputWorkspaces",
                                              Direction::Input),
      "Input Workspaces. List of histogram workspaces to stitch together.");

  declareProperty(make_unique<WorkspaceProperty<Workspace>>(
                      "OutputWorkspace", "", Direction::Output),
                  "Output stitched workspace.");

  declareProperty(make_unique<ArrayProperty<double>>(
                      "Params", boost::make_shared<RebinParamsValidator>(true),
                      Direction::Input),
                  "Rebinning Parameters. See Rebin for format.");

  declareProperty(
      make_unique<ArrayProperty<double>>("StartOverlaps", Direction::Input),
      "Start overlaps for stitched workspaces.");

  declareProperty(
      make_unique<ArrayProperty<double>>("EndOverlaps", Direction::Input),
      "End overlaps for stitched workspaces.");

  declareProperty(make_unique<PropertyWithValue<bool>>("ScaleRHSWorkspace",
                                                       true, Direction::Input),
                  "Scaling either with respect to workspace 1 or workspace 2");

  declareProperty(make_unique<PropertyWithValue<bool>>("UseManualScaleFactor",
                                                       false, Direction::Input),
                  "True to use a provided value for the scale factor.");

  auto manualScaleFactorValidator =
      boost::make_shared<BoundedValidator<double>>();
  manualScaleFactorValidator->setLower(0);
  manualScaleFactorValidator->setExclusive(true);
  declareProperty(make_unique<PropertyWithValue<double>>(
                      "ManualScaleFactor", 1.0, manualScaleFactorValidator,
                      Direction::Input),
                  "Provided value for the scale factor.");

  declareProperty(
      make_unique<ArrayProperty<double>>("OutScaleFactors", Direction::Output),
      "The actual used values for the scaling factors at each stitch step.");

  auto scaleFactorFromPeriodValidator =
      boost::make_shared<BoundedValidator<size_t>>();
  scaleFactorFromPeriodValidator->setLower(0);
  declareProperty(make_unique<PropertyWithValue<size_t>>(
                      "ScaleFactorFromPeriod", 0,
                      scaleFactorFromPeriodValidator, Direction::Input),
                  "Provided index of period to obtain scale factor from.");
}

/** Load and validate the algorithm's properties.
 */
std::map<std::string, std::string> Stitch1DMany::validateInputs() {
  std::map<std::string, std::string> errors;

  const std::vector<std::string> inputWorkspacesStr =
      this->getProperty("InputWorkspaces");
  if (inputWorkspacesStr.size() < 2)
    errors["InputWorkspaces"] = "At least 2 input workspaces required.";

  for (const auto &ws : inputWorkspacesStr) {
    if (AnalysisDataService::Instance().doesExist(ws)) {
      m_inputWorkspaces.push_back(
          AnalysisDataService::Instance().retrieveWS<Workspace>(ws));
    } else {
      errors["InputWorkspaces"] = ws + " is not a valid workspace.";
      break;
    }
  }

  // Check that all the workspaces are of the same type
  if (!m_inputWorkspaces.empty()) {
    const std::string id = m_inputWorkspaces[0]->id();
    for (auto &inputWorkspace : m_inputWorkspaces) {
      if (inputWorkspace->id() != id) {
        errors["InputWorkspaces"] = "All workspaces must be the same type.";
        break;
      }
    }
  } else {
    errors["InputWorkspaces"] = "Input workspaces must be given";
  }

  m_numWorkspaces = m_inputWorkspaces.size();

  m_startOverlaps = this->getProperty("StartOverlaps");
  m_endOverlaps = this->getProperty("EndOverlaps");

  if (!m_startOverlaps.empty() && m_startOverlaps.size() != m_numWorkspaces - 1)
    errors["StartOverlaps"] = "If given, StartOverlaps must have one fewer "
                              "entries than the number of input workspaces.";

  if (m_startOverlaps.size() != m_endOverlaps.size())
    errors["EndOverlaps"] =
        "EndOverlaps must have the same number of entries as StartOverlaps.";

  m_scaleRHSWorkspace = this->getProperty("ScaleRHSWorkspace");
  m_useManualScaleFactor = this->getProperty("UseManualScaleFactor");
  m_manualScaleFactor = this->getProperty("ManualScaleFactor");
  m_params = this->getProperty("Params");

  if (m_params.empty())
    errors["Params"] = "At least one parameter must be given.";

  if (!m_scaleRHSWorkspace) {
    // Flip these around for processing
    std::reverse(m_inputWorkspaces.begin(), m_inputWorkspaces.end());
    std::reverse(m_startOverlaps.begin(), m_startOverlaps.end());
    std::reverse(m_endOverlaps.begin(), m_endOverlaps.end());
  }

  return errors;
}

/** Load and validate the algorithm's properties for workspace groups.
 */
void Stitch1DMany::validateGroupWorkspacesInputs() {
  std::string error;

  const std::vector<std::string> inputWorkspacesStr =
    this->getProperty("InputWorkspaces");
  if (inputWorkspacesStr.size() < 2)
    throw std::runtime_error("At least 2 input workspace groups required.");

  for (const auto &ws : inputWorkspacesStr) {
    if (AnalysisDataService::Instance().doesExist(ws)) {
      m_inputWorkspaceGroups.push_back(
        AnalysisDataService::Instance().retrieveWS<WorkspaceGroup>(ws));
    }
    else {
      throw std::runtime_error(ws + " is not a valid workspace group.");
    }
  }

  // Check all workspace groups are the same size
  size_t groupSize = m_inputWorkspaceGroups[0]->size();
  for (auto &inputWsGroup : m_inputWorkspaceGroups) {
    if (inputWsGroup->size() != groupSize) {
      throw std::runtime_error("All workspace groups must be the same size.");
    }
  }
}

/** Execute the algorithm.
 */
void Stitch1DMany::exec() {
  MatrixWorkspace_sptr lhsWS =
      boost::dynamic_pointer_cast<MatrixWorkspace>(m_inputWorkspaces[0]);

  for (size_t i = 1; i < m_numWorkspaces; ++i) {
    MatrixWorkspace_sptr rhsWS =
        boost::dynamic_pointer_cast<MatrixWorkspace>(m_inputWorkspaces[i]);
    double outScaleFactor;

    doStitch1D(lhsWS, rhsWS, i, m_startOverlaps, m_endOverlaps, m_params,
        m_scaleRHSWorkspace, m_useManualScaleFactor, m_manualScaleFactor, lhsWS,
        outScaleFactor);

    m_scaleFactors.push_back(outScaleFactor);
  }

  if (!isChild()) {
    // Copy each input workspace's history into our output workspace's history
    for (auto &inputWorkspace : m_inputWorkspaces)
      lhsWS->history().addHistory(inputWorkspace->getHistory());
  }
  // We're a child algorithm, but we're recording history anyway
  else if (isRecordingHistoryForChild() && m_parentHistory) {
    m_parentHistory->addChildHistory(m_history);
  }

  m_outputWorkspace = lhsWS;

  // Save output
  this->setProperty("OutputWorkspace", m_outputWorkspace);
  this->setProperty("OutScaleFactors", m_scaleFactors);
}

/** Performs the Stitch1D algorithm at a specific workspace index
*/
void Stitch1DMany::doStitch1D(MatrixWorkspace_sptr lhsWS,
    MatrixWorkspace_sptr rhsWS, size_t wsIndex,
    std::vector<double> startOverlaps, std::vector<double> endOverlaps,
    std::vector<double> params, bool scaleRhsWS, bool useManualScaleFactor,
    double manualScaleFactor, MatrixWorkspace_sptr &outWS,
    double &outScaleFactor) {

  IAlgorithm_sptr alg = createChildAlgorithm("Stitch1D");
  alg->initialize();
  alg->setProperty("LHSWorkspace", lhsWS);
  alg->setProperty("RHSWorkspace", rhsWS);
  if (startOverlaps.size() > wsIndex - 1) {
    alg->setProperty("StartOverlap", startOverlaps[wsIndex - 1]);
    alg->setProperty("EndOverlap", endOverlaps[wsIndex - 1]);
  }
  alg->setProperty("Params", params);
  alg->setProperty("ScaleRHSWorkspace", scaleRhsWS);
  alg->setProperty("UseManualScaleFactor", useManualScaleFactor);
  if (useManualScaleFactor)
    alg->setProperty("ManualScaleFactor", manualScaleFactor);
  alg->execute();

  outWS = alg->getProperty("OutputWorkspace");
  outScaleFactor = alg->getProperty("OutScaleFactor");
}

bool Stitch1DMany::checkGroups() {
  std::vector<std::string> wsNames = getProperty("InputWorkspaces");

  try {
    if (AnalysisDataService::Instance().retrieveWS<WorkspaceGroup>(wsNames[0]))
      return true;
  }
  catch (...) { }
  return false;
}

bool Stitch1DMany::processGroups() {
  validateGroupWorkspacesInputs();

  // List of workspaces to be grouped
  std::vector<std::string> toGroup;

  const std::string groupName = this->getProperty("OutputWorkspace");

  size_t numWSPerGroup = m_inputWorkspaceGroups[0]->size();

  for (size_t i = 0; i < numWSPerGroup; ++i) {
    // List of workspaces to stitch
    std::vector<std::string> toProcess;
    // The name of the resulting workspace
    std::string outName = groupName;

    for (auto &groupWs : m_inputWorkspaceGroups) {
      const std::string wsName = groupWs->getItem(i)->name();
      toProcess.push_back(wsName);
      outName += "_" + wsName;
    }

    IAlgorithm_sptr stitchAlg = createChildAlgorithm("Stitch1DMany");
    stitchAlg->initialize();
    stitchAlg->setAlwaysStoreInADS(true);
    for (auto &prop : this->getProperties()) {
      stitchAlg->setProperty(prop->name(), prop->value());
    }
    stitchAlg->setProperty("InputWorkspaces", toProcess);
    stitchAlg->setProperty("OutputWorkspace", outName);
    stitchAlg->execute();

    // Add the resulting workspace to the list to be grouped together
    toGroup.push_back(outName);

    // Add the scalefactors to the list so far
    const std::vector<double> scaleFactors =
      stitchAlg->getProperty("OutScaleFactors");
    m_scaleFactors.insert(m_scaleFactors.end(), scaleFactors.begin(),
      scaleFactors.end());
  }

  IAlgorithm_sptr groupAlg = createChildAlgorithm("GroupWorkspaces");
  groupAlg->initialize();
  groupAlg->setAlwaysStoreInADS(true);
  groupAlg->setProperty("InputWorkspaces", toGroup);
  groupAlg->setProperty("OutputWorkspace", groupName);
  groupAlg->execute();

  m_outputWorkspace =
    AnalysisDataService::Instance().retrieveWS<Workspace>(groupName);

  this->setProperty("OutputWorkspace", m_outputWorkspace);
  this->setProperty("OutScaleFactors", m_scaleFactors);
  return true;
}
} // namespace Algorithms
} // namespace Mantid

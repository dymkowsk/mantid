#include "MantidAPI/MatrixWorkspace.h"
#include "MantidDataHandling/LoadInstrument.h"
#include "MantidDataObjects/EventWorkspace.h"
#include "MantidGeometry/IInstrument.h"
#include "MantidKernel/cow_ptr.h"
#include "MantidKernel/DateAndTime.h"
#include "MantidKernel/Utils.h"
#include "MantidMDEvents/BoxController.h"
#include "MantidMDEvents/MDEventWorkspace.h"
#include "MantidTestHelpers/DLLExport.h"
#include "MantidTestHelpers/MDEventsTestHelper.h"
#include "MantidTestHelpers/WorkspaceCreationHelper.h"
#include "MantidGeometry/MDGeometry/MDTypes.h"


using Mantid::DataObjects::EventWorkspace_sptr;
using Mantid::Kernel::DateAndTime;
using Mantid::DataHandling::LoadInstrument;
using Mantid::DataObjects::EventWorkspace;

namespace Mantid
{
namespace MDEvents
{


/** Set of helper methods for testing MDEventWorkspace things
 *
 * @author Janik Zikovsky
 * @date March 29, 2011
 * */
namespace MDEventsTestHelper
{

  //-------------------------------------------------------------------------------------
  /** Create an EventWorkspace containing fake data
   * of single-crystal diffraction.
   * Instrument is MINITOPAZ
   *
   * @return EventWorkspace_sptr
   */
  EventWorkspace_sptr createDiffractionEventWorkspace(int numEvents)
  {
    Mantid::Kernel::ConfigService::Instance().setString("default.facility", "TEST");
    int numPixels = 10000;
    int numBins = 1600;
    double binDelta = 10.0;

    EventWorkspace_sptr retVal(new EventWorkspace);
    retVal->initialize(numPixels,1,1);

    // --------- Load the instrument -----------
    LoadInstrument * loadInst = new LoadInstrument();
    loadInst->initialize();
    loadInst->setPropertyValue("Filename", "IDFs_for_UNIT_TESTING/MINITOPAZ_Definition.xml");
    loadInst->setProperty<Mantid::API::MatrixWorkspace_sptr> ("Workspace", retVal);
    loadInst->execute();
    delete loadInst;
    // Populate the instrument parameters in this workspace - this works around a bug
    retVal->populateInstrumentParameters();

    DateAndTime run_start("2010-01-01");

    for (int pix = 0; pix < numPixels; pix++)
    {
      for (int i=0; i<numEvents; i++)
      {
        retVal->getEventListAtPixelID(pix) += Mantid::DataObjects::TofEvent((i+0.5)*binDelta, run_start+double(i));
      }

    }
    retVal->doneLoadingData();

    //Create the x-axis for histogramming.
    Mantid::MantidVecPtr x1;
    Mantid::MantidVec& xRef = x1.access();
    xRef.resize(numBins);
    for (int i = 0; i < numBins; ++i)
    {
      xRef[i] = i*binDelta;
    }

    //Set all the histograms at once.
    retVal->setAllX(x1);

    // Give it a crystal and goniometer
    WorkspaceCreationHelper::SetGoniometer(retVal, 0., 0., 0.);
    WorkspaceCreationHelper::SetOrientedLattice(retVal, 1., 1., 1.);

    // Some sanity checks
    if (retVal->getInstrument()->getName() != "MINITOPAZ")
      throw std::runtime_error("MDEventsTestHelper::createDiffractionEventWorkspace(): Wrong instrument loaded.");
    Mantid::detid2det_map dets;
    retVal->getInstrument()->getDetectors(dets);
    if ( dets.size() != 100*100)
      throw std::runtime_error("MDEventsTestHelper::createDiffractionEventWorkspace(): Wrong instrument size.");


    return retVal;
  }



  //-------------------------------------------------------------------------------------
  /** Generate an empty MDBox */
  MDBox<MDEvent<1>,1> * makeMDBox1(size_t splitInto)
  {
    // Split at 5 events
    BoxController_sptr splitter(new BoxController(1));
    splitter->setSplitThreshold(5);
    // Splits into 10 boxes
    splitter->setSplitInto(splitInto);
    // Set the size
    MDBox<MDEvent<1>,1> * out = new MDBox<MDEvent<1>,1>(splitter);
    out->setExtents(0, 0.0, 10.0);
    out->calcVolume();
    return out;
  }

  //-------------------------------------------------------------------------------------
  /** Generate an empty MDBox with 3 dimensions, split 10x5x2 */
  MDBox<MDEvent<3>,3> * makeMDBox3()
  {
    // Split at 5 events
    BoxController_sptr splitter(new BoxController(3));
    splitter->setSplitThreshold(5);
    // Splits into 10x5x2 boxes
    splitter->setSplitInto(10);
    splitter->setSplitInto(1,5);
    splitter->setSplitInto(2,2);
    // Set the size to 10.0 in all directions
    MDBox<MDEvent<3>,3> * out = new MDBox<MDEvent<3>,3>(splitter);
    for (size_t d=0; d<3; d++)
      out->setExtents(d, 0.0, 10.0);
    return out;
  }


  //-------------------------------------------------------------------------------------
  /** Return a vector with this many MDEvents, spaced evenly from 0.5, 1.5, etc. */
  std::vector<MDEvent<1> > makeMDEvents1(size_t num)
  {
    std::vector<MDEvent<1> > out;
    for (double i=0; i<num; i++)
    {
      coord_t coords[1] = {i*1.0+0.5};
      out.push_back( MDEvent<1>(1.0, 1.0, coords) );
    }
    return out;
  }



} // namespace
}
}

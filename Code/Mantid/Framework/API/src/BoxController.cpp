#include "MantidKernel/Strings.h"
#include "MantidKernel/System.h"
#include "MantidKernel/VectorHelper.h"
#include "MantidAPI/BoxController.h"
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <Poco/File.h>
#include <Poco/DOM/Attr.h>
#include <Poco/DOM/AutoPtr.h>
#include <Poco/DOM/Document.h>
#include <Poco/DOM/DOMParser.h>
#include <Poco/DOM/DOMWriter.h>
#include <Poco/DOM/Element.h>
#include <Poco/DOM/Text.h>
#include <Poco/XML/XMLWriter.h>
#include <sstream>

using namespace Mantid::Kernel;
using Mantid::Kernel::Strings::convert;
using Mantid::Kernel::VectorHelper::splitStringIntoVector;

namespace Mantid
{
namespace API
{

  //-----------------------------------------------------------------------------------
  /** Copy constructor
   */
  BoxController::BoxController(const BoxController & other)
  : nd(other.nd), m_maxId(other.m_maxId),
    m_SplitThreshold(other.m_SplitThreshold),
    m_maxDepth(other.m_maxDepth), m_splitInto(other.m_splitInto),
    m_numSplit(other.m_numSplit),
    m_addingEvents_eventsPerTask(other.m_addingEvents_eventsPerTask),
    m_addingEvents_numTasksPerBlock(other.m_addingEvents_numTasksPerBlock),
    m_numMDBoxes(other.m_numMDBoxes),
    m_numMDGridBoxes(other.m_numMDGridBoxes),
    m_maxNumMDBoxes(other.m_maxNumMDBoxes),
    m_filename(other.m_filename),
    m_file(other.m_file),
    m_bytesPerEvent(other.m_bytesPerEvent)
  {
  }

  /// Destructor
  BoxController::~BoxController()
  {
    this->closeFile();
  }
   /**reserve range of id-s for use on set of adjacent boxes. 
    * Needed to be thread safe as adjacent boxes have to have subsequent ID-s
    * @param range  --range number of box-id-s to lock
    * @returns initial ID to use in the range
    */
   size_t BoxController::claimIDRange(size_t range)
   {
     m_idMutex.lock();
     size_t tmp=m_maxId;
     m_maxId+=range; 
     m_idMutex.unlock();
     return tmp;
   }
  /** Serialize to an XML string
   * @return XML string
   */
  std::string BoxController::toXMLString() const
  {
    using namespace Poco::XML;

    //Create the root element for this fragment.
    AutoPtr<Document> pDoc = new Document;
    AutoPtr<Element> pBoxElement = pDoc->createElement("BoxController");
    pDoc->appendChild(pBoxElement);

    AutoPtr<Element> element;
    AutoPtr<Text> text;
    std::string vecStr;

    element = pDoc->createElement("NumDims");
    element->appendChild( pDoc->createTextNode(boost::str(boost::format("%d") % this->getNDims())) );
    pBoxElement->appendChild(element);

    element = pDoc->createElement("MaxId");
    element->appendChild( pDoc->createTextNode(boost::str(boost::format("%d") % this->getMaxId())) );
    pBoxElement->appendChild(element);

    element = pDoc->createElement("SplitThreshold");
    element->appendChild( pDoc->createTextNode(boost::str(boost::format("%d") % this->getSplitThreshold())) );
    pBoxElement->appendChild(element);

    element = pDoc->createElement("MaxDepth");
    element->appendChild( pDoc->createTextNode(boost::str(boost::format("%d") % this->getMaxDepth())) );
    pBoxElement->appendChild(element);

    element = pDoc->createElement("SplitInto");
    vecStr = Kernel::Strings::join( this->m_splitInto.begin(), this->m_splitInto.end(), ",");
    element->appendChild( pDoc->createTextNode( vecStr ) );
    pBoxElement->appendChild(element);

    element = pDoc->createElement("NumMDBoxes");
    vecStr = Kernel::Strings::join( this->m_numMDBoxes.begin(), this->m_numMDBoxes.end(), ",");
    element->appendChild( pDoc->createTextNode( vecStr ) );
    pBoxElement->appendChild(element);

    element = pDoc->createElement("NumMDGridBoxes");
    vecStr = Kernel::Strings::join( this->m_numMDGridBoxes.begin(), this->m_numMDGridBoxes.end(), ",");
    element->appendChild( pDoc->createTextNode( vecStr ) );
    pBoxElement->appendChild(element);

    //Create a string representation of the DOM tree.
    std::stringstream xmlstream;
    DOMWriter writer;
    writer.writeNode(xmlstream, pDoc);

    return xmlstream.str().c_str();
  }


  //------------------------------------------------------------------------------------------------------
  /** Static method that sets the data inside this BoxController from an XML string
   *
   * @param xml :: string generated by BoxController::toXMLString()
   */
  void BoxController::fromXMLString(const std::string & xml)
  {
    using namespace Poco::XML;
    Poco::XML::DOMParser pParser;
    Poco::XML::Document* pDoc = pParser.parseString(xml);
    Poco::XML::Element* pBoxElement = pDoc->documentElement();

    std::string s;
    s = pBoxElement->getChildElement("NumDims")->innerText();
    Strings::convert(s, nd);
    if (nd <= 0 || nd > 20) throw std::runtime_error("BoxController::fromXMLString(): Bad number of dimensions found.");

    size_t ival;
    Strings::convert(pBoxElement->getChildElement("MaxId")->innerText(), ival);
    this->setMaxId(ival);
    Strings::convert(pBoxElement->getChildElement("SplitThreshold")->innerText(), ival);
    this->setSplitThreshold(ival);
    Strings::convert(pBoxElement->getChildElement("MaxDepth")->innerText(), ival);
    this->setMaxDepth(ival);

    s = pBoxElement->getChildElement("SplitInto")->innerText();
    this->m_splitInto = splitStringIntoVector<size_t>(s);

    s = pBoxElement->getChildElement("NumMDBoxes")->innerText();
    this->m_numMDBoxes = splitStringIntoVector<size_t>(s);

    s = pBoxElement->getChildElement("NumMDGridBoxes")->innerText();
    this->m_numMDGridBoxes = splitStringIntoVector<size_t>(s);

    this->calcNumSplit();
  }


  //------------------------------------------------------------------------------------------------------
  /** Close the open file for the back-end, if any
   * Note: this does not save any data that might be, e.g., in the MRU.
   * @param deleteFile :: if true, will delete the file. Default false.
   */
  void BoxController::closeFile(bool deleteFile)
  {
    if (m_file)
    {
      m_file->close();
      m_file = NULL;
    }
    if (deleteFile && !m_filename.empty())
    {
      Poco::File file(m_filename);
      if (file.exists()) file.remove();
      m_filename = "";
    }
  }

} // namespace Mantid

} // namespace API




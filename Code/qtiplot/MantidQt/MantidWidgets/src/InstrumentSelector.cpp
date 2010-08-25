//------------------------------------------------------
// Includes
//------------------------------------------------------
#include "MantidQtMantidWidgets/InstrumentSelector.h"
#include "MantidKernel/ConfigService.h"
#include "MantidKernel/FacilityInfo.h"
#include "MantidKernel/InstrumentInfo.h"
#include "MantidKernel/Exception.h"

namespace MantidQt
{
  namespace MantidWidgets
  {
     using namespace Mantid::Kernel;

    //------------------------------------------------------
    // Public member functions
    //------------------------------------------------------

    /**
     * Default constructor
     * @param parent A widget to act as this widget's parent (default = NULL)
     * @param init If true then the widget will be populated with the instrument list (default = true)
     */
    InstrumentSelector::InstrumentSelector(QWidget *parent, bool init) : QComboBox(parent), m_techniques(), m_currentFacility(NULL)
    {
      setEditable(false);
      if( init )
      {
        fillWithInstrumentsFromFacility();
        connect(this, SIGNAL(currentIndexChanged(const QString &)), this, SLOT(updateDefaultInstrument(const QString &)));
        connect(this, SIGNAL(currentIndexChanged(const QString &)), this, SIGNAL(instrumentSelectionChanged(const QString &)));
      }
    }

    /**
     * Return the list of techniques that are supported by the instruments in the widget
     * @returns A list of supported techniques
     */
    QStringList InstrumentSelector::getTechniques() const
    {
      return m_techniques;
    }

    /**
     * Set the list of techniques
     * @param techniques Only those instruments that support these techniques will be shown
     */
    void InstrumentSelector::setTechniques(const QStringList & techniques)
    {
      m_techniques = techniques;
      if( count() > 0 && m_currentFacility ) 
      {
        filterByTechniquesAtFacility(techniques, *m_currentFacility);
      }      
    }

    //------------------------------------------------------
    // Public slot member functions
    //------------------------------------------------------

     /**
      * Populate list with instruments from the named facility. Note the current list is cleared.
      * @param name The name of the facility whose instruments should be placed in the list. An empty string uses the default
      * facility defined in Mantid.
      */
    void InstrumentSelector::fillWithInstrumentsFromFacility(const QString & name)
    {
      ConfigServiceImpl & mantidSettings = ConfigService::Instance(); 
      
      clear();
      if( name.isEmpty() )
      {
        m_currentFacility = &(mantidSettings.Facility());
      }
      else
      {
        m_currentFacility = &(mantidSettings.Facility(name.toStdString()));
      }

      const std::vector<InstrumentInfo> & instruments = m_currentFacility->Instruments();
      std::vector<InstrumentInfo>::const_iterator iend = instruments.end();
      std::set<std::string> alphabetizedNames;
      for( std::vector<InstrumentInfo>::const_iterator itr = instruments.begin(); itr != iend; ++itr )
      {
        alphabetizedNames.insert(itr->name());
      }
      std::set<std::string>::const_iterator namesEnd = alphabetizedNames.end();
      for( std::set<std::string>::const_iterator itr = alphabetizedNames.begin(); itr != namesEnd; ++itr )
      {
        QString name = QString::fromStdString(*itr);
        std::string prefix = m_currentFacility->Instrument(*itr).shortName();
        QString shortName = QString::fromStdString(prefix);
        this->addItem(name, QVariant(shortName));
      }
      filterByTechniquesAtFacility(m_techniques, *m_currentFacility);

      QString defaultName;
      try
      {
        defaultName = QString::fromStdString(m_currentFacility->Instrument().name());
      }
      catch( Exception::NotFoundError &)
      {
        defaultName = "";
      }
      int index = this->findText(defaultName);
      if( index < 0 )
      {
        index = 0;
      }
      // Don't affect the default instrument
      this->blockSignals(true);
      this->setCurrentIndex(index);
      this->blockSignals(false);
    }

    //------------------------------------------------------
    // Private slot member functions
    //------------------------------------------------------
    /**
     * Set the named instrument as the default for Mantid  
     * @param name A string containing the new instrument to set as the default 
     */
    void InstrumentSelector::updateDefaultInstrument(const QString & name) const
    {
      if( !name.isEmpty() )
      {
        ConfigService::Instance().setString("default.instrument", name.toStdString());
      }
    }

    //------------------------------------------------------
    // Private non-slot member functions
    //------------------------------------------------------

    /**
     * Filter the list to only show those supporting the given technique
     * @param techniques A string list containing the names of a techniques to filter the instrument list by
     * @param facility A FacilityInfo object
     */
    void InstrumentSelector::filterByTechniquesAtFacility(const QStringList & techniques, const Mantid::Kernel::FacilityInfo & facility)
    {
      if( techniques.isEmpty() ) return;

      QStringList supportedInstruments;
      QStringListIterator techItr(techniques);
      while( techItr.hasNext() )
      {
        const std::vector<InstrumentInfo> instruments = facility.Instruments(techItr.next().toStdString());
        const size_t nInstrs = instruments.size();
        for( size_t i = 0; i < nInstrs; ++i )
        {
          supportedInstruments.append(QString::fromStdString(instruments[i].name()));
        }
      }

      // Remove those not supported
      for( int i = 0 ; i < this->count(); )
      {
        if( !supportedInstruments.contains(itemText(i)) )
        {
          removeItem(i);
        }
        else
        {
          ++i;
        }
      }
    }

  }
}
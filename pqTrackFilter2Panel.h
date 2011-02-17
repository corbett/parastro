#ifndef _pqTrackFilter2Panel_h
#define _pqTrackFilter2Panel_h

#include "pqLoadedFormObjectPanel.h"
#include "pqObjectPanelInterface.h"

#include <QObject>

class pqDisplayPanel;
class pqDoubleRangeWidget;

class pqTrackFilter2Panel : public pqLoadedFormObjectPanel
{
  Q_OBJECT

  //typedef pqLoadedFormObjectPanel Superclass;
  //typedef QObject Superclass;
  
public:
  pqTrackFilter2Panel(pqProxy* proxy, QWidget* p);
  ~pqTrackFilter2Panel();

private:
	//Ui::TrackFilter2Panel Widgets;
};

#endif
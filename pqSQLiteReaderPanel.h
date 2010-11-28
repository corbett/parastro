
#ifndef _pqSQLiteReaderPanel_h
#define _pqSQLiteReaderPanel_h

#include "pqLoadedFormObjectPanel.h"
#include "pqObjectPanelInterface.h"

class pqSQLiteReaderPanel : public pqLoadedFormObjectPanel
{
  Q_OBJECT
public:
  pqSQLiteReaderPanel(pqProxy* proxy, QWidget* p);
  ~pqSQLiteReaderPanel();
};

#endif


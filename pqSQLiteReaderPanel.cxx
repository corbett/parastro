
#include "pqSQLiteReaderPanel.h"

#include <QLayout>
#include <QLabel>

#include "vtkSMProxy.h"
#include "pqProxy.h"

pqSQLiteReaderPanel::pqSQLiteReaderPanel(pqProxy* pxy, QWidget* p)
  : pqLoadedFormObjectPanel(":/pqSQLiteReaderPanel.ui", pxy, p)
{
  //this->layout()->addWidget(new QLabel("This is from a plugin", this));
}

pqSQLiteReaderPanel::~pqSQLiteReaderPanel()
{
}



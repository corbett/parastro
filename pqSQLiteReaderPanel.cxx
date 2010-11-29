
#include "pqSQLiteReaderPanel.h"

#include <QLayout>
#include <QLabel>

#include "vtkSMProxy.h"
#include "pqProxy.h"

#include "QRadioButton.h"
#include "QCheckBox.h"


pqSQLiteReaderPanel::pqSQLiteReaderPanel(pqProxy* pxy, QWidget* p)
  : pqLoadedFormObjectPanel(":/pqSQLiteReaderPanel.ui", pxy, p)
{
  //this->layout()->addWidget(new QLabel("This is from a plugin", this));
	this->DisplayMode = this->findChild<QRadioButton*>("DisplayMode_0");
	this->DisplayAllData1 = this->findChild<QCheckBox*>("DisplayAllData1");
	this->linkServerManagerProperties();
}

pqSQLiteReaderPanel::~pqSQLiteReaderPanel()
{
}



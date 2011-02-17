#include "pqTrackFilter2Panel.h"

#include "pqProxy.h"
#include "vtkSMProxy.h"

#include <QLabel>
#include <QLayout>

#include "pqDoubleRangeWidget.h"

/*
#include "vtkSMProperty.h"



//#include <pqPropertyHelper.h>
#include "pqSMAdaptor.h"



#include <QComboBox>
#include <QMessageBox>

#include <iostream>
#include <set>
*/

 
pqTrackFilter2Panel::pqTrackFilter2Panel(pqProxy* proxy, QWidget* p) :
  pqLoadedFormObjectPanel(":/ParaViewResources/pqTrackFilter2Panel.ui", proxy, p)
{

	this->layout()->addWidget(new QLabel("This is from a plugin", this));
	this->layout()->addWidget(new pqDoubleRangeWidget(this));

	this->linkServerManagerProperties();

}

pqTrackFilter2Panel::~pqTrackFilter2Panel()
{
}
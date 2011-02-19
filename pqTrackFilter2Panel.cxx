#include "pqTrackFilter2Panel.h"

#include "pqProxy.h"
#include "vtkSMProxy.h"

#include <QLayout>

#include "pqDoubleRangeWidget.h"
#include "pqSMAdaptor.h"

/*
#include "vtkSMProperty.h"



//#include <pqPropertyHelper.h>



#include <QMessageBox>

#include <iostream>
#include <set>
*/

 
pqTrackFilter2Panel::pqTrackFilter2Panel(pqProxy* proxy, QWidget* p) :
  pqLoadedFormObjectPanel(":/ParaViewResources/pqTrackFilter2Panel.ui", proxy, p)
{
	//pqDoubleRangeWidget * newRW;
	QWidget * widget;
	QComboBox * comboB;
	QGridLayout *gridlayoutF = this->findChild<QGridLayout*>("filterlayout");
	QGridLayout *gridlayoutR = this->findChild<QGridLayout*>("restrictionlayout");

	this->RestrictionArray = this->findChild<QComboBox*>("RestrictionArray");;
	this->label4 = this->findChild<QLabel*>("label_4");
	this->label5 = this->findChild<QLabel*>("label_5");

// adjust GUI

	widget = this->findChild<QWidget*>("Filter_0");
	widget->hide();
	this->Filter_0 = new pqDoubleRangeWidget();
	this->Filter_0->setObjectName(QString::fromUtf8("Filter_0"));
	QObject::connect(
		this->Filter_0,
		SIGNAL(valueEdited(double)),
		this,
		SLOT(lowerFChanged(double))
		);
	gridlayoutF->addWidget(this->Filter_0,0,1,1,1);


	widget = this->findChild<QWidget*>("Filter_1");
	widget->hide();
	this->Filter_1 = new pqDoubleRangeWidget();
	this->Filter_1->setObjectName(QString::fromUtf8("Filter_1"));
	QObject::connect(
		this->Filter_1,
		SIGNAL(valueEdited(double)),
		this,
		SLOT(upperFChanged(double))
		);
	gridlayoutF->addWidget(this->Filter_1,2,1,1,1);


	widget = this->findChild<QWidget*>("Restriction_0");
	widget->hide();
	this->Restriction_0 = new pqDoubleRangeWidget();
	this->Restriction_0->setObjectName(QString::fromUtf8("Restriction_0"));
	QObject::connect(
		this->Restriction_0,
		SIGNAL(valueEdited(double)),
		this,
		SLOT(lowerRChanged(double))
		);
	gridlayoutR->addWidget(this->Restriction_0,0,1,1,1);

	widget = this->findChild<QWidget*>("Restriction_1");
	widget->hide();
	this->Restriction_1 = new pqDoubleRangeWidget();
	this->Restriction_1->setObjectName(QString::fromUtf8("Restriction_1"));
	QObject::connect(
		this->Restriction_1,
		SIGNAL(valueEdited(double)),
		this,
		SLOT(upperRChanged(double))
		);
	gridlayoutR->addWidget(this->Restriction_1,2,1,1,1);


	comboB = this->findChild<QComboBox*>("ModeSelection");
	comboB->addItem("some point on the track");
	comboB->addItem("every point on the track");
	comboB->addItem("the points on the track, where");
	QObject::connect(
		comboB,
		SIGNAL(currentIndexChanged(int)),
		this,
		SLOT(selectionModeChanged(int))
		);


	comboB = this->findChild<QComboBox*>("FilterArray");
	comboB->addItem("f1");
	comboB->addItem("f2");
	comboB->addItem("f3");


	comboB = this->findChild<QComboBox*>("RestrictionArray");
	comboB->addItem("r1");
	comboB->addItem("r2");
	comboB->addItem("r3");

	vtkSMProperty* prop = this->proxy()->GetProperty("FilterArray");
	double val = pqSMAdaptor::getElementProperty(prop).toDouble();
	pqSMAdaptor::setElementProperty(prop,5);




	//this->layout()->addWidget(new QLabel("This is from a plugin", this));
	//this->layout()->addWidget(new pqDoubleRangeWidget(this));

	this->linkServerManagerProperties();

}

pqTrackFilter2Panel::~pqTrackFilter2Panel()
{
}

/*
void pqTrackFilter2Panel::accept()
{
	Superclass::accept();    
}

void pqTrackFilter2Panel::reset()
{
	Superclass::reset();
}
*/
void pqTrackFilter2Panel::lowerFChanged(double val)
{
  // clamp the lower threshold if we need to
	if(this->Filter_1->value() < val)
	{
		this->Filter_1->setValue(val);
	}
}

void pqTrackFilter2Panel::upperFChanged(double val)
{
  // clamp the lower threshold if we need to
	if(this->Filter_0->value() > val)
    {
		this->Filter_0->setValue(val);
    }
}

void pqTrackFilter2Panel::lowerRChanged(double val)
{
	// clamp the lower threshold if we need to
	if(this->Restriction_1->value() < val)
	{
		this->Restriction_1->setValue(val);
	}
}

void pqTrackFilter2Panel::upperRChanged(double val)
{
	// clamp the lower threshold if we need to
	if(this->Restriction_0->value() > val)
	{
		this->Restriction_0->setValue(val);
	}
}
void pqTrackFilter2Panel::variableChanged()
{
	// when the user changes the variable, adjust the ranges on the ThresholdBetween
	vtkSMProperty* prop = this->proxy()->GetProperty("ThresholdBetween");
	QList<QVariant> range = pqSMAdaptor::getElementPropertyDomain(prop);
	if(range.size() == 2 && range[0].isValid() && range[1].isValid())
	{
		this->Filter_0->setValue(range[0].toDouble());
		this->Filter_1->setValue(range[1].toDouble());
	}
}

void pqTrackFilter2Panel::selectionModeChanged(int selection)
{
	if (selection==2)
	{
		this->Restriction_0->setEnabled(true);
		this->Restriction_1->setEnabled(true);
		this->RestrictionArray->setEnabled(true);
		this->label4->setEnabled(true);
		this->label5->setEnabled(true);
	}
	else
	{
		this->Restriction_0->setEnabled(false);
		this->Restriction_1->setEnabled(false);
		this->RestrictionArray->setEnabled(false);
		this->label4->setEnabled(false);
		this->label5->setEnabled(false);
	}
}
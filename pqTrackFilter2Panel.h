#ifndef _pqTrackFilter2Panel_h
#define _pqTrackFilter2Panel_h

#include "pqLoadedFormObjectPanel.h"
#include "pqObjectPanelInterface.h"

#include "pqPropertyManager.h"
#include "pqNamedWidgets.h"

#include <QDockWidget>
#include <QObject>
#include <QLabel>
#include <QComboBox>



class pqDisplayPanel;
class pqDoubleRangeWidget;

class pqTrackFilter2Panel : public pqLoadedFormObjectPanel
{
	Q_OBJECT

	typedef pqLoadedFormObjectPanel Superclass;
	//typedef QObject Superclass;
  
public:
	pqTrackFilter2Panel(pqProxy* proxy, QWidget* p);
	~pqTrackFilter2Panel();

protected slots:
	void lowerFChanged(double);
	void upperFChanged(double);
	void lowerRChanged(double);
	void upperRChanged(double);
	void variableChanged();
	void selectionModeChanged(int);

private slots:
	//virtual void accept();
	//virtual void reset();

protected:
	pqDoubleRangeWidget* Filter_0;
	pqDoubleRangeWidget* Filter_1;
	pqDoubleRangeWidget* Restriction_0;
	pqDoubleRangeWidget* Restriction_1;
	QComboBox *RestrictionArray;
	QLabel *label4;
	QLabel *label5;


private:
	//Ui::pqTrackFilter2Panel Widgets;


};

#endif
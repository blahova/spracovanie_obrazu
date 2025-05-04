#pragma once

#include <QtWidgets/QMainWindow>
#include <QtWidgets>
#include "ui_ImageViewer.h"
#include "ViewerWidget.h"
#include "ImageProcessing.h"

class ImageViewer : public QMainWindow
{
	Q_OBJECT

public:
	ImageViewer(QWidget* parent = Q_NULLPTR);

private:
	Ui::ImageViewerClass* ui;
	ViewerWidget* vW;

	QSettings settings;
	QMessageBox msgBox;

	QImage original;
	ImageProcessing img_proc;

	//ImageViewer Events
	void closeEvent(QCloseEvent* event);

	//Image functions
	bool openImage(QString filename);
	bool saveImage(QString filename);
	bool invertColors();
	bool showProcImg(int time=-1);

private slots:
	void on_actionOpen_triggered();
	void on_actionSave_as_triggered();
	void on_actionExit_triggered();
	void on_actionInvert_triggered();
	void on_actionEH_triggered();
	void on_actionFSHS_triggered();
	void on_actionOriginal_triggered();
	void on_actionKonvolucia_triggered();
	void on_actionHranovy_detektor_triggered();
	void on_actionVzdialenostna_funkcia_triggered();

	void on_pushButton_RVT_clicked();
	void on_spinBox_T_lin_valueChanged(int i);
	void on_horizontalSlider_lin_valueChanged(int i);
	void on_spinBox_stepIMG_lin_valueChanged(int i);

	void on_pushButton_nelin_clicked();
	void on_horizontalSlider_nelin_valueChanged(int i);
	void on_spinBox_stepIMG_nelin_valueChanged(int i);

	void on_pushButton_MCF_clicked();
	void on_pushButton_GMCF_clicked();

};

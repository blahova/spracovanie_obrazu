#include "ImageViewer.h"


ImageViewer::ImageViewer(QWidget* parent)
	: QMainWindow(parent), ui(new Ui::ImageViewerClass)
{
	ui->setupUi(this);
	vW = new ViewerWidget(QSize(500, 500));
	ui->scrollArea->setWidget(vW);

	ui->scrollArea->setBackgroundRole(QPalette::Dark);
	ui->scrollArea->setWidgetResizable(true);
	ui->scrollArea->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
	ui->scrollArea->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);

	vW->setObjectName("ViewerWidget");
}

//ImageViewer Events
void ImageViewer::closeEvent(QCloseEvent* event)
{
	if (QMessageBox::Yes == QMessageBox::question(this, "Close Confirmation", "Are you sure you want to exit?", QMessageBox::Yes | QMessageBox::No))
	{
		event->accept();
	}
	else {
		event->ignore();
	}
}

//Image functions
bool ImageViewer::openImage(QString filename)
{
	QImage loadedImg(filename);
	if (!loadedImg.isNull()) {
		original = loadedImg.copy();
		img_proc = ImageProcessing(original.bits(), original.width(), original.height());

		return vW->setImage(loadedImg);
	}
	return false;
}
bool ImageViewer::saveImage(QString filename)
{
	QFileInfo fi(filename);
	QString extension = fi.completeSuffix();

	QImage* img = vW->getImage();
	return img->save(filename, extension.toStdString().c_str());
}

bool ImageViewer::invertColors()
{
	if (vW->isEmpty()) {
		return false;
	}

	uchar* data = vW->getData();

	int row = vW->getImage()->bytesPerLine();
	int depth = vW->getImage()->depth();

	for (int i = 0; i < vW->getImgHeight(); i++)
	{
		for (int j = 0; j < vW->getImgWidth(); j++)
		{
			//Grayscale
			if (depth == 8) {
				vW->setPixel(j, i, static_cast<uchar>(255 - data[i * row + j]));
			}
			//RGBA
			else {
				uchar r = static_cast<uchar>(255 - data[i * row + j * 4]);
				uchar g = static_cast<uchar>(255 - data[i * row + j * 4 + 1]);
				uchar b = static_cast<uchar>(255 - data[i * row + j * 4 + 2]);
				vW->setPixel(j, i, r, g, b);
			}
		}
	}
	vW->update();
	return true;
}

bool ImageViewer::showProcImg(int time)
{
	//ak je time=true, tak potom sa zobrazuje casovy krok z pola krokov, ak false, tak sa zobrazuje procImage


	if (vW->isEmpty()) {
		return false;
	}

	double** procData = nullptr;

	if (time==-1) procData = img_proc.getProcData();
	else procData = img_proc.getTimeData(time);

	if (!procData) return false;

	int width = img_proc.getwidth();
	int height = img_proc.getheight();

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			uchar pixelValue = static_cast<uchar>(procData[i][j] * 255.0+0.5);
			vW->setPixel(j, i, pixelValue);
		}
	}

	vW->update();
	return true;
}


//Slots
void ImageViewer::on_actionOpen_triggered()
{
	QString folder = settings.value("folder_img_load_path", "").toString();

	QString fileFilter = "Image data (*.bmp *.gif *.jpg *.jpeg *.png *.pbm *.pgm *.ppm .*xbm .* xpm);;All files (*)";
	QString fileName = QFileDialog::getOpenFileName(this, "Load image", folder, fileFilter);
	if (fileName.isEmpty()) { return; }

	QFileInfo fi(fileName);
	settings.setValue("folder_img_load_path", fi.absoluteDir().absolutePath());

	if (!openImage(fileName)) {
		msgBox.setText("Unable to open image.");
		msgBox.setIcon(QMessageBox::Warning);
		msgBox.exec();
	}
	ui->horizontalSlider_lin->blockSignals(true);
	ui->horizontalSlider_lin->setValue(ui->horizontalSlider_lin->maximum());
	ui->horizontalSlider_lin->blockSignals(false);
	ui->horizontalSlider_lin->setEnabled(false);

	ui->spinBox_stepIMG_lin->blockSignals(true);
	ui->spinBox_stepIMG_lin->setValue(ui->spinBox_stepIMG_lin->maximum());
	ui->spinBox_stepIMG_lin->blockSignals(false);
	ui->spinBox_stepIMG_lin->setEnabled(false);

	ui->horizontalSlider_nelin->blockSignals(true);
	ui->horizontalSlider_nelin->setValue(ui->horizontalSlider_nelin->maximum());
	ui->horizontalSlider_nelin->blockSignals(false);
	ui->horizontalSlider_nelin->setEnabled(false);

	ui->spinBox_stepIMG_nelin->blockSignals(true);
	ui->spinBox_stepIMG_nelin->setValue(ui->spinBox_stepIMG_nelin->maximum());
	ui->spinBox_stepIMG_nelin->blockSignals(false);
	ui->spinBox_stepIMG_nelin->setEnabled(false);

}
void ImageViewer::on_actionSave_as_triggered()
{
	QString folder = settings.value("folder_img_save_path", "").toString();

	QString fileFilter = "Image data (*.bmp *.gif *.jpg *.jpeg *.png *.pbm *.pgm *.ppm .*xbm .* xpm);;All files (*)";
	QString fileName = QFileDialog::getSaveFileName(this, "Save image", folder, fileFilter);
	if (!fileName.isEmpty()) {
		QFileInfo fi(fileName);
		settings.setValue("folder_img_save_path", fi.absoluteDir().absolutePath());

		if (!saveImage(fileName)) {
			msgBox.setText("Unable to save image.");
			msgBox.setIcon(QMessageBox::Warning);
		}
		else {
			msgBox.setText(QString("File %1 saved.").arg(fileName));
			msgBox.setIcon(QMessageBox::Information);
		}
		msgBox.exec();
	}
}
void ImageViewer::on_actionExit_triggered()
{
	this->close();
}

void ImageViewer::on_actionInvert_triggered()
{
	invertColors();
}

void ImageViewer::on_actionEH_triggered()
{
	img_proc.EH();
	showProcImg();
}

void ImageViewer::on_actionFSHS_triggered()
{
	img_proc.FSHS();
	showProcImg();
}

void ImageViewer::on_actionOriginal_triggered()
{
	vW->setImage(original);
}

void ImageViewer::on_actionKonvolucia_triggered()
{
	if (img_proc.getProcData() == nullptr) return;

	img_proc.konvolucia();
	showProcImg();
}

void ImageViewer::on_actionHranovy_detektor_triggered()
{
	img_proc.hranovy();
	showProcImg();
}

void ImageViewer::on_actionVzdialenostna_funkcia_triggered()
{
	img_proc.vzdialenostna();
	showProcImg();
}

void ImageViewer::on_pushButton_RVT_clicked()
{
	if (!img_proc.getProcData()) return;
	int T = ui->spinBox_T_lin->value();
	double tau = ui->doubleSpinBox_tau_lin->value();

	if (tau <= 0.2)
	{
		img_proc.explicitna(T, tau);
	}
	else
	{
		img_proc.implicitna(T, tau);
	}

	ui->horizontalSlider_lin->setEnabled(true);
	ui->spinBox_stepIMG_lin->setEnabled(true);

	ui->horizontalSlider_lin->setRange(0, ui->spinBox_T_lin->value());
	ui->horizontalSlider_lin->setValue(ui->spinBox_T_lin->value());
	ui->spinBox_stepIMG_lin->setMaximum(ui->spinBox_T_lin->value());
	ui->spinBox_stepIMG_lin->setValue(ui->spinBox_T_lin->value());

	ui->horizontalSlider_nelin->blockSignals(true);
	ui->horizontalSlider_nelin->setValue(ui->horizontalSlider_nelin->maximum());
	ui->horizontalSlider_nelin->blockSignals(false);
	ui->horizontalSlider_nelin->setEnabled(false);

	ui->spinBox_stepIMG_nelin->blockSignals(true);
	ui->spinBox_stepIMG_nelin->setValue(ui->spinBox_stepIMG_nelin->maximum());
	ui->spinBox_stepIMG_nelin->blockSignals(false);
	ui->spinBox_stepIMG_nelin->setEnabled(false);

	showProcImg();
}

void ImageViewer::on_spinBox_T_lin_valueChanged(int i)
{
	
}

void ImageViewer::on_horizontalSlider_lin_valueChanged(int i)
{
	if (ui->spinBox_stepIMG_lin->value() != i)
	{
		ui->spinBox_stepIMG_lin->blockSignals(true);
		ui->spinBox_stepIMG_lin->setValue(i);
		ui->spinBox_stepIMG_lin->blockSignals(false);
	}
	showProcImg(i);
}

void ImageViewer::on_spinBox_stepIMG_lin_valueChanged(int i)
{
	if (ui->horizontalSlider_lin->value() != i) 
	{
		ui->horizontalSlider_lin->blockSignals(true);
		ui->horizontalSlider_lin->setValue(i);
		ui->horizontalSlider_lin->blockSignals(false);
	}
	showProcImg(i);
}

void ImageViewer::on_pushButton_nelin_clicked()
{
	if (!img_proc.getProcData()) return;

	int T = ui->spinBox_T_nelin->value();
	double tau = ui->doubleSpinBox_tau_nelin->value();
	double sigma = ui->doubleSpinBox_sigma_nelin->value();
	double K = ui->doubleSpinBox_K_nelin->value();

	img_proc.semi_implicitna(T, tau, sigma, K);

	ui->horizontalSlider_nelin->setEnabled(true);
	ui->spinBox_stepIMG_nelin->setEnabled(true);

	ui->horizontalSlider_nelin->setRange(0, ui->spinBox_T_nelin->value());
	ui->horizontalSlider_nelin->setValue(ui->spinBox_T_nelin->value());
	ui->spinBox_stepIMG_nelin->setMaximum(ui->spinBox_T_nelin->value());
	ui->spinBox_stepIMG_nelin->setValue(ui->spinBox_T_nelin->value());


	ui->horizontalSlider_lin->blockSignals(true);
	ui->horizontalSlider_lin->setValue(ui->horizontalSlider_lin->maximum());
	ui->horizontalSlider_lin->blockSignals(false);
	ui->horizontalSlider_lin->setEnabled(false);

	ui->spinBox_stepIMG_lin->blockSignals(true);
	ui->spinBox_stepIMG_lin->setValue(ui->spinBox_stepIMG_lin->maximum());
	ui->spinBox_stepIMG_lin->blockSignals(false);
	ui->spinBox_stepIMG_lin->setEnabled(false);

	showProcImg();


}

void ImageViewer::on_horizontalSlider_nelin_valueChanged(int i)
{
	if (ui->spinBox_stepIMG_nelin->value() != i)
	{
		ui->spinBox_stepIMG_nelin->blockSignals(true);
		ui->spinBox_stepIMG_nelin->setValue(i);
		ui->spinBox_stepIMG_nelin->blockSignals(false);
	}
	showProcImg(i);
}

void ImageViewer::on_spinBox_stepIMG_nelin_valueChanged(int i)
{
	if (ui->horizontalSlider_nelin->value() != i)
	{
		ui->horizontalSlider_nelin->blockSignals(true);
		ui->horizontalSlider_nelin->setValue(i);
		ui->horizontalSlider_nelin->blockSignals(false);
	}
	showProcImg(i);
}

void ImageViewer::on_pushButton_MCF_clicked()
{
	if (!img_proc.getProcData()) return;
	int T = ui->spinBox_T_lin->value();
	double tau = ui->doubleSpinBox_tau_lin->value();

	img_proc.MCF(T,tau);

	ui->horizontalSlider_lin->setEnabled(true);
	ui->spinBox_stepIMG_lin->setEnabled(true);

	ui->horizontalSlider_lin->setRange(0, ui->spinBox_T_lin->value());
	ui->horizontalSlider_lin->setValue(ui->spinBox_T_lin->value());
	ui->spinBox_stepIMG_lin->setMaximum(ui->spinBox_T_lin->value());
	ui->spinBox_stepIMG_lin->setValue(ui->spinBox_T_lin->value());

	ui->horizontalSlider_nelin->blockSignals(true);
	ui->horizontalSlider_nelin->setValue(ui->horizontalSlider_nelin->maximum());
	ui->horizontalSlider_nelin->blockSignals(false);
	ui->horizontalSlider_nelin->setEnabled(false);

	ui->spinBox_stepIMG_nelin->blockSignals(true);
	ui->spinBox_stepIMG_nelin->setValue(ui->spinBox_stepIMG_nelin->maximum());
	ui->spinBox_stepIMG_nelin->blockSignals(false);
	ui->spinBox_stepIMG_nelin->setEnabled(false);

	showProcImg();
}

void ImageViewer::on_pushButton_GMCF_clicked()
{
	if (!img_proc.getProcData()) return;

	int T = ui->spinBox_T_nelin->value();
	double tau = ui->doubleSpinBox_tau_nelin->value();
	double sigma = ui->doubleSpinBox_sigma_nelin->value();
	double K = ui->doubleSpinBox_K_nelin->value();

	img_proc.GMCF(T, tau, sigma, K);

	ui->horizontalSlider_nelin->setEnabled(true);
	ui->spinBox_stepIMG_nelin->setEnabled(true);

	ui->horizontalSlider_nelin->setRange(0, ui->spinBox_T_nelin->value());
	ui->horizontalSlider_nelin->setValue(ui->spinBox_T_nelin->value());
	ui->spinBox_stepIMG_nelin->setMaximum(ui->spinBox_T_nelin->value());
	ui->spinBox_stepIMG_nelin->setValue(ui->spinBox_T_nelin->value());


	ui->horizontalSlider_lin->blockSignals(true);
	ui->horizontalSlider_lin->setValue(ui->horizontalSlider_lin->maximum());
	ui->horizontalSlider_lin->blockSignals(false);
	ui->horizontalSlider_lin->setEnabled(false);

	ui->spinBox_stepIMG_lin->blockSignals(true);
	ui->spinBox_stepIMG_lin->setValue(ui->spinBox_stepIMG_lin->maximum());
	ui->spinBox_stepIMG_lin->blockSignals(false);
	ui->spinBox_stepIMG_lin->setEnabled(false);

	showProcImg();
}




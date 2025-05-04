#pragma once
#include <QVector>

//vytvorit triedu ImageProcessing kde sa bude ukladat obrazok ako kopia z imagevieweru ale nie ako unsigned char ale ako double v rozsahu 0-1

class ImageProcessing {
private:
	double** imageData = nullptr;	//originalne data v 0-1
	double** procData = nullptr;	//data ktore boli neakym sopsobom upravene

	int width, height;

	std::vector<int>histogram;

	std::vector<double**> steps;

public:
	ImageProcessing() {};
	ImageProcessing(uchar* data, int w, int h);

	double** getProcData() const { return procData; }
	double** getTimeData(int i)const { return steps[i]; }
	void mirror(int d, double** img_data, double** mirrored);
	void FSHS();
	void EH();
	void konvolucia();
	void explicitna(int T,double tau);
	void implicitna(int T, double tau);
	void hranovy();
	void implicitna_sigma(double sigma, double** vstup,double** vysledok);
	void explicitna_sigma(double sigma, double** vstup, double** vysledok);
	void semi_implicitna(int T, double tau, double sigma, double K);
	void MCF(int T, double tau);
	void GMCF(int T, double tau, double sigma, double K);
	void vzdialenostna();

	double find_min();
	double find_max();
	double mean(double** img);


	int getwidth() { return width; };
	int getheight() { return height; };
};
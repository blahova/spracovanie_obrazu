#include "ImageProcessing.h"
#include <iostream>
#include <fstream>
#include <iomanip>




ImageProcessing::ImageProcessing(uchar* data, int w, int h)
{
    width = w;
    height = h;
    
    imageData = new double* [height];
    procData = new double* [height];
    histogram.assign(256,0);

    for (int i = 0; i < height; ++i) {
        imageData[i] = new double[width];
        procData[i] = new double[width];
        for (int j = 0; j < width; ++j) {
            imageData[i][j] = data[i * width + j] / 255.0;
            histogram[data[i * width + j]]++;
        }
    }

    //std::ofstream outFile("test.txt");

    //if (!outFile) {
    //    std::cerr << "Error opening file for writing.\n";
    //    return;
    //}

    //for (int i = 0; i < height; ++i) {
    //    for (int j = 0; j < width; ++j) {
    //        outFile << imageData[i][j] << " ";
    //    }
    //    outFile << "\n";
    //}

    //outFile.close();

}



void ImageProcessing::FSHS()
{
    if (!imageData)return;
    double min = find_min();
    double max = find_max();

    for (int i = 0; i < height; ++i) 
    {
        for (int j = 0; j < width; ++j) {
            procData[i][j] = (imageData[i][j] - min) / (max - min);
        }
    }

}

void ImageProcessing::EH()
{
    if (!imageData)return;
    int numPixels = width * height;
    std::vector<double> cdf(256, 0); //kumulativna distribucna funkcia z histogramu

    cdf[0] = histogram[0];
    for (int i = 1; i < 256; i++) {
        cdf[i] = cdf[i - 1] + histogram[i];
    }

    for (int i = 0; i < 256; i++) {
        cdf[i] /= numPixels;
    }

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            int oldIntensity = (int)(imageData[i][j] * 255.0 +0.5);//index
            procData[i][j] = cdf[oldIntensity]; 
        }
    }
}

void ImageProcessing::konvolucia()
{
    int d = 3;

    int newHeight = height + 2 * d;
    int newWidth = width + 2 * d;

    double** mirrored=new double* [newHeight];
    for (int i = 0; i < newHeight; i++) {
        mirrored[i] = new double[newWidth];
        for (int j = 0; j < newWidth; ++j) {
            mirrored[i][j] = 0;
        }
    }

    mirror(d, imageData,mirrored); //tu urobim zrkadlenie, teraz uz mozem robit konvoluciu

    double mask[7][7] = {
        {0.000035724962295, 0.000362193668795, 0.001444830988919, 0.002288756141850, 0.001444830988919, 0.000362193668795, 0.000035724962295},
        {0.000362193668795, 0.003672061362370, 0.014648262812588, 0.023204306757593, 0.014648262812588, 0.003672061362370, 0.000362193668795},
        {0.001444830988919, 0.014648262812588, 0.058433556047159, 0.092564570748284, 0.058433556047159, 0.014648262812588, 0.001444830988919},
        {0.002288756141850, 0.023204306757593, 0.092564570748284, 0.147561796159379, 0.092564570748284, 0.023204306757593, 0.002288756141850},
        {0.001444830988919, 0.014648262812588, 0.058433556047159, 0.092564570748284, 0.058433556047159, 0.014648262812588, 0.001444830988919},
        {0.000362193668795, 0.003672061362370, 0.014648262812588, 0.023204306757593, 0.014648262812588, 0.003672061362370, 0.000362193668795},
        {0.000035724962295, 0.000362193668795, 0.001444830988919, 0.002288756141850, 0.001444830988919, 0.000362193668795, 0.000035724962295}
    };

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            double sum = 0.0;

            for (int mi = 0; mi < 2 * d + 1; mi++) {
                for (int mj = 0; mj < 2 * d + 1; mj++) {
                    if ((i + mi) >= 0 && (i + mi) < newHeight && (j + mj) >= 0 && (j + mj) < newWidth && mi <7 && mj<7) {
                        sum += mirrored[i + mi][j + mj] * mask[mi][mj];
                    }
                }
            }
            procData[i][j] = sum;

        }
    }

}

void ImageProcessing::implicitna(int T, double tau)
{
    int d = 1;
    int newHeight = height + 2 * d;
    int newWidth = width + 2 * d;

    //najprv zmazem ak uz boli nejake data predtym ulozene, ulozene su v std::vector<double**> steps;
    for (double** step : steps) {
        for (int i = 0; i < height; i++) {
            delete[] step[i];
        }
        delete[] step;
    }
    steps.clear();

    steps.resize(T + 1);

    for (int t = 0; t <= T; t++)
    {
        steps[t] = new double* [height];
        for (int i = 0; i < height; i++)
        {
            steps[t][i] = new double[width]();
        }
    }
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            steps[0][i][j] = imageData[i][j];
        }
    }

    double** img_old = new double* [newHeight];
    double** img_new = new double* [height];

    for (int i = 0; i < newHeight; i++) {
        img_old[i] = new double[newWidth]();
    }

    for (int i = 0; i < height; i++) {
        img_new[i] = new double[width]();
    }

    double sum = 0.0;
    mirror(d, imageData, img_old);
 
    double omega = 1.25;
    double sigma = 0.0;
    double tolerancia = 1.0E-6;
    int max_iter = 10000;
    double rezid;
    double* b=new double[newHeight*newWidth]; //prava strana
    double* phi= new double[newHeight * newWidth]();


    for (int i = 0; i < newHeight; i++)
    {
        for (int j = 0; j < newWidth; j++)
        {
            b[i * newWidth + j] = img_old[i][j];
            phi[i * newWidth + j] = b[i * newWidth + j];
            if(i>=d and i<newHeight-d and j>=d and j<newWidth-d) sum += img_old[i][j];
        }
    }

    std::cout << "Povodna stredna hodnota:" <<std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << sum / (width * height) << std::endl;
    
    for (int cas = 1; cas <= T; cas++)
    {
        int iter=0;
        sum = 0.0;
        for(iter=1;iter<max_iter;iter++)
        {
            rezid = 0.0;
            for (int i = d; i < height + d; i++) {
                for (int j = d; j < width + d; j++) {
                    sigma = 0.0;
                    sigma += -tau * (phi[i*newWidth+j+1] + phi[(i-1)*newWidth+j] + phi[(i+1)*newWidth+j] + phi[i*newWidth+j-1]);

                    int idx = i * newWidth + j; 
                    phi[idx] = (1 - omega) * phi[idx] + omega * (b[idx] - sigma) / (1 + 4 * tau);
                }
            }

            //zapis vysledku iteracie do 2d
            for (int i = d; i < height+d; i++)
            {
                for (int j = d; j < width+d; j++)
                {
                    img_new[i-d][j-d] = phi[i * newWidth + j];
                }
            }
            mirror(d, img_new, img_old);//mirror novo vypocitaneho obrazku
            for (int i = 0; i < newHeight; i++)//prepis na 1d aby sa dalo pocitat v 1d
            {
                for (int j = 0; j < newWidth; j++)
                {
                    phi[i * newWidth + j] = img_old[i][j];
                }
            }

            for (int i = d; i < height + d; i++) {//pocitanie rezidua
                for (int j = d; j < width + d; j++) {
                    int idx = i * newWidth + j;
                    double Ax = (1 + 4 * tau) * phi[idx] - tau * (phi[i * newWidth + j + 1] + phi[(i - 1) * newWidth + j] + phi[(i + 1) * newWidth + j] + phi[i * newWidth + j - 1]) - b[idx];
                    rezid += Ax * Ax;
                }
            }
            //std::cout << sqrt(rezid) << std::endl;
            if (sqrt(rezid) < tolerancia) break;  
        }
        
        std::cout << "casovy krok: " << cas << " iteracia: " << iter << " reziduum: " << sqrt(rezid) << std::endl;

        for (int i = 0; i < newHeight; i++) {//update pravej strany
            for (int j = 0; j < newWidth; j++) {
                b[i * newWidth + j] = img_old[i][j]; 
                if (i >= d and i < newHeight - d and j >= d and j < newWidth-d) sum += img_old[i][j];
            }
        }
        std::cout << "Stredna hodnota:" << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << sum / (width * height) << std::endl;

        for (int i = 0; i < height; i++) {//zapis jednotlivych krokov aby sa dali po jednom vykreslit
            for (int j = 0; j < width; j++) {
                steps[cas][i][j] = img_new[i][j];
            }
        }

    }
    //ulozenie posledneho na vykreslenie hned
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            procData[i][j] = img_new[i][j];

        }
    }

    for (int i = 0; i < newHeight; i++) {
        delete[] img_old[i];
    }
    delete[] img_old;

    for (int i = 0; i < height; i++) {
        delete[] img_new[i];
    }
    delete[] img_new;
    delete[] b;
    delete[] phi;

}

void ImageProcessing::explicitna(int T, double tau)
{
    int d = 1;
    int newHeight = height + 2 * d;
    int newWidth = width + 2 * d;

    //najprv zmazem ak uz boli nejake data predtym ulozene, ulozene su v std::vector<double**> steps;
    for (double** step : steps) {
        for (int i = 0; i < height; i++) {
            delete[] step[i];
        }
        delete[] step;
    }
    steps.clear();


    steps.resize(T + 1);

    for (int t = 0; t <= T; t++) 
    {
        steps[t] = new double* [height];
        for (int i = 0; i < height; i++)
        {
            steps[t][i] = new double[width]();  
        }
    }


    double** img_old = new double* [newHeight];
    double** img_new = new double* [height];

    for (int i = 0; i < newHeight; i++) {
        img_old[i] = new double[newWidth]();
    }

    for (int i = 0; i < height; i++) {
        img_new[i] = new double[width]();
    }

    double sum = 0.0;
    mirror(d, imageData, img_old);
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            steps[0][i][j] = imageData[i][j];

            sum += imageData[i][j];
        }
    }
    std::cout << "Step 0" << " average: " << sum / (height * width) << std::endl;

    for (int cas = 1; cas <= T; cas++)
    {
        sum = 0.0;
        for (int i = d; i < height + d; i++)
        {
            for (int j = d; j < width + d; j++)
            {
                img_new[i - d][j - d] = (1 - 4 * tau) * img_old[i][j] + tau*(img_old[i - 1][j] + img_old[i + 1][j] + img_old[i][j - 1] + img_old[i][j + 1]);
                sum += img_new[i - d][j - d];
            }
        }
        std::cout << "Step " << cas << " average: " << sum / (height * width) << std::endl;

        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                steps[cas][i][j] = img_new[i][j];
            }
        }
        mirror(d, img_new, img_old);  
    }

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            procData[i][j] = img_new[i][j];
        }
    }

    for (int i = 0; i < newHeight; i++) {
        delete[] img_old[i];
    }
    delete[] img_old;

    for (int i = 0; i < height; i++) {
        delete[] img_new[i]; 
    }
    delete[] img_new;
}

void ImageProcessing::implicitna_sigma(double sigma, double** vstup, double** vysledok)
{
    //vstup je velkosti height x width, vystup bude ten vacsi (ako zmirrorovany) 
    int d = 1;
    int newHeight = height + 2 * d;
    int newWidth = width + 2 * d;

    double** img_old = new double* [newHeight];
    double** img_new = new double* [height];

    for (int i = 0; i < newHeight; i++) {
        img_old[i] = new double[newWidth]();
    }

    for (int i = 0; i < height; i++) {
        img_new[i] = new double[width]();
    }

    double sum = 0.0;
    mirror(d, vstup, img_old);

    double omega = 1.25;
    double sigmaSOR = 0.0;
    double tolerancia = 1.0E-6;
    int max_iter = 10000;
    double rezid;
    double* b = new double[newHeight * newWidth]; //prava strana
    double* phi = new double[newHeight * newWidth]();


    for (int i = 0; i < newHeight; i++)
    {
        for (int j = 0; j < newWidth; j++)
        {
            b[i * newWidth + j] = img_old[i][j];
            phi[i * newWidth + j] = b[i * newWidth + j];
            if (i >= d and i < newHeight - d and j >= d and j < newWidth - d) sum += img_old[i][j];
        }
    }

    //std::cout << " IMPLICITNA\n Povodna stredna hodnota:" << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << sum / (width * height) << std::endl;

    int iter = 0;
    sum = 0.0;
    for (iter = 1; iter < max_iter; iter++)
    {
        rezid = 0.0;
        for (int i = d; i < height + d; i++) {
            for (int j = d; j < width + d; j++) {
                sigmaSOR = 0.0;
                sigmaSOR += -sigma * (phi[i * newWidth + j + 1] + phi[(i - 1) * newWidth + j] + phi[(i + 1) * newWidth + j] + phi[i * newWidth + j - 1]);

                int idx = i * newWidth + j;
                phi[idx] = (1 - omega) * phi[idx] + omega * (b[idx] - sigmaSOR) / (1 + 4 * sigma);
            }
        }


        for (int i = d; i < height + d; i++)
        {
            for (int j = d; j < width + d; j++)
            {
                img_new[i - d][j - d] = phi[i * newWidth + j];
            }
        }
        mirror(d, img_new, img_old);
        for (int i = 0; i < newHeight; i++)
        {
            for (int j = 0; j < newWidth; j++)
            {
                phi[i * newWidth + j] = img_old[i][j];
            }
        }

        for (int i = d; i < height + d; i++) {
            for (int j = d; j < width + d; j++) {
                int idx = i * newWidth + j;
                double Ax = (1 + 4 * sigma) * phi[idx] - sigma * (phi[i * newWidth + j + 1] + phi[(i - 1) * newWidth + j] + phi[(i + 1) * newWidth + j] + phi[i * newWidth + j - 1]) - b[idx];
                rezid += Ax * Ax;
            }
        }
        //std::cout << sqrt(rezid) << std::endl;
        if (sqrt(rezid) < tolerancia) break;
    }

    //std::cout <<" reziduum: " << sqrt(rezid) << std::endl;

    for (int i = 0; i < newHeight; i++) {
        for (int j = 0; j < newWidth; j++) {
            b[i * newWidth + j] = img_old[i][j];
            if (i >= d and i < newHeight - d and j >= d and j < newWidth - d) sum += img_old[i][j];
        }
    }
    //std::cout << "IMPLICITNA\nStredna hodnota:" << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << sum / (width * height) << std::endl;

    for (int i = 0; i < newHeight; i++) {
        for (int j = 0; j < newWidth; j++) {
            vysledok[i][j] = img_old[i][j];
        }
    }

    for (int i = 0; i < newHeight; i++) {
        delete[] img_old[i];
    }
    delete[] img_old;

    for (int i = 0; i < height; i++) {
        delete[] img_new[i];
    }
    delete[] img_new;
    delete[] b;
    delete[] phi;

}

void ImageProcessing::semi_implicitna(int T, double tau, double sigma,double K) //perona-malik
{
    int d = 1;
    int newHeight = height + 2 * d;
    int newWidth = width + 2 * d;

    double sum = 0.0;

    for (double** step : steps) {
        for (int i = 0; i < height; i++) {
            delete[] step[i];
        }
        delete[] step;
    }
    steps.clear();

    steps.resize(T + 1);

    for (int t = 0; t <= T; t++)
    {
        steps[t] = new double* [height];
        for (int i = 0; i < height; i++)
        {
            steps[t][i] = new double[width]();
        }
    }

    double** img_before = new double* [height]; //onrazok pred pripravou implicitnou/explicitnou


    for (int i = 0; i < height; i++) {
        img_before[i] = new double[width]();
    }
    
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            steps[0][i][j] = imageData[i][j];
            img_before[i][j] = imageData[i][j];
            sum += imageData[i][j];
        }
    }

    double** prepared_img = new double* [newHeight];    //obrazok prejdeny implicitnou/explicitnou
    for (int i = 0; i < newHeight; i++) {
        prepared_img[i] = new double[newWidth]();
    }


    double omega = 1.25;
    double sigmaSOR = 0.0;
    double tolerancia = 1.0E-6;
    int max_iter = 10000;
    double rezid;
    double* b = new double[newHeight * newWidth](); //prava strana
    double* phi = new double[newHeight * newWidth]();

    double uE, uW, uS, uN, uNE, uNW, uSW, uSE, uxE, uxN, uxW, uxS, uyE, uyN, uyW, uyS, uC, uE_grad, uN_grad, uW_grad, uS_grad;
    double gE, gN, gW, gS;

    if (sigma <= 0.2)
    {
        explicitna_sigma(sigma, img_before, prepared_img);
    }
    else
    {
        implicitna_sigma(sigma, img_before, prepared_img);
    }

    for (int i = 0; i < newHeight; i++)
    {
        for (int j = 0; j < newWidth; j++)
        {
            b[i * newWidth + j] = prepared_img[i][j];
            phi[i * newWidth + j] = b[i * newWidth + j];
        }
    }

    std::cout << "Povodna stredna hodnota:" << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << sum / (width * height) << std::endl;

    for (int cas = 1; cas <= T; cas++)
    {
        int iter = 0;
        sum = 0.0;
        for (iter = 1; iter < max_iter; iter++)
        {
            rezid = 0.0;
            for (int i = d; i < height + d; i++) {
                for (int j = d; j < width + d; j++) {
                    sigmaSOR = 0.0;
                    //najst g-cka pre vsetky strany E,N,W,S
                    uE = prepared_img[i][j + 1];
                    uN = prepared_img[i - 1][j];
                    uW = prepared_img[i][j - 1];
                    uS = prepared_img[i + 1][j];

                    uNE = prepared_img[i - 1][j + 1];
                    uNW = prepared_img[i - 1][j - 1];
                    uSW = prepared_img[i + 1][j - 1];
                    uSE = prepared_img[i + 1][j + 1];

                    uC = prepared_img[i][j];//u current index (stred)

                    uxE = uE - uC;
                    uyE = (uNE - uSE - uS + uN) / 4;

                    uxW = uW - uC;
                    uyW = (uNW + uN - uS - uSW) / 4;

                    uyN = uN - uC;
                    uxN = (uNE + uE - uNW - uW) / 4;

                    uyS = uS - uC;
                    uxS = (uE + uSE - uSW - uW) / 4;

                    uE_grad = uxE * uxE + uyE * uyE;
                    uW_grad = uxW * uxW + uyW * uyW;
                    uN_grad = uxN * uxN + uyN * uyN;
                    uS_grad = uxS * uxS + uyS * uyS;

                    gE = 1.0 / (1.0 + K * uE_grad);
                    gN = 1.0 / (1.0 + K * uN_grad);
                    gW = 1.0 / (1.0 + K * uW_grad);
                    gS = 1.0 / (1.0 + K * uS_grad);

                    sigmaSOR += -tau * (gE*phi[i * newWidth + j + 1] + gN*phi[(i - 1) * newWidth + j] + gS*phi[(i + 1) * newWidth + j] + gW*phi[i * newWidth + j - 1]);

                    int idx = i * newWidth + j;
                    phi[idx] = (1 - omega) * phi[idx] + omega * (b[idx] - sigmaSOR) / (1 + tau*(gE+gN+gW+gS));
                }
            }

            //zapis vysledku iteracie do 2d
            for (int i = d; i < height + d; i++)
            {
                for (int j = d; j < width + d; j++)
                {
                    img_before[i - d][j - d] = phi[i * newWidth + j];
                }
            }
            mirror(d, img_before, prepared_img);
            for (int i = 0; i < newHeight; i++)
            {
                for (int j = 0; j < newWidth; j++)
                {
                    phi[i * newWidth + j] = prepared_img[i][j];
                }
            }

            //rezidua
            for (int i = d; i < height + d; i++) {
                for (int j = d; j < width + d; j++) {
                    uE = prepared_img[i][j + 1];
                    uN = prepared_img[i - 1][j];
                    uW = prepared_img[i][j - 1];
                    uS = prepared_img[i + 1][j];

                    uNE = prepared_img[i - 1][j + 1];
                    uNW = prepared_img[i - 1][j - 1];
                    uSW = prepared_img[i + 1][j - 1];
                    uSE = prepared_img[i + 1][j + 1];

                    uC = prepared_img[i][j];//u current index (stred)

                    uxE = uE - uC;
                    uyE = (uNE - uSE - uS + uN) / 4;

                    uxW = uW - uC;
                    uyW = (uNW + uN - uS - uSW) / 4;

                    uyN = uN - uC;
                    uxN = (uNE + uE - uNW - uW) / 4;

                    uyS = uS - uC;
                    uxS = (uE + uSE - uSW - uW) / 4;

                    uE_grad = uxE * uxE + uyE * uyE;
                    uW_grad = uxW * uxW + uyW * uyW;
                    uN_grad = uxN * uxN + uyN * uyN;
                    uS_grad = uxS * uxS + uyS * uyS;

                    gE = 1 / (1 + K * uE_grad);
                    gN = 1 / (1 + K * uN_grad);
                    gW = 1 / (1 + K * uW_grad);
                    gS = 1 / (1 + K * uS_grad);

                    int idx = i * newWidth + j;
                    double Ax = (1 + tau * (gE + gN + gW + gS)) * phi[idx] - tau * (gE * phi[i * newWidth + j + 1] + gN * phi[(i - 1) * newWidth + j] + gS * phi[(i + 1) * newWidth + j] + gW * phi[i * newWidth + j - 1]) - b[idx];
                    rezid += Ax * Ax;
                }
            }
            //std::cout << sqrt(rezid) << std::endl;
            if (sqrt(rezid) < tolerancia) break;
        }

        std::cout << "casovy krok: " << cas << " iteracia: " << iter << " reziduum: " << sqrt(rezid) << std::endl;

        for (int i = 0; i < newHeight; i++) {//update pravej strany
            for (int j = 0; j < newWidth; j++) {
                b[i * newWidth + j] = prepared_img[i][j];
                if (i >= d and i < newHeight - d and j >= d and j < newWidth - d) sum += prepared_img[i][j];
            }
        }
        std::cout << "Stredna hodnota:" << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << sum / (width * height) << std::endl;

        for (int i = d; i < height+d; i++) {//zapis jednotlivych krokov aby sa dali po jednom vykreslit
            for (int j = d; j < width+d; j++) {
                img_before[i - d][j - d] = prepared_img[i][j];
            }
        }

        for (int i = 0; i < height; i++) {//zapis jednotlivych krokov aby sa dali po jednom vykreslit
            for (int j = 0; j < width; j++) {
                steps[cas][i][j] = img_before[i][j];
            }
        }

        if (sigma <= 0.2)
        {
            explicitna_sigma(sigma, img_before, prepared_img);
        }
        else
        {
            implicitna_sigma(sigma, img_before, prepared_img);
        }
    }

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            procData[i][j] = img_before[i][j];

        }
    }

    for (int i = 0; i < newHeight; i++) {
        delete[] prepared_img[i];
    }
    delete[] prepared_img;

    for (int i = 0; i < height; i++) {
        delete[] img_before[i];
    }
    delete[] img_before;
    delete[] b;
    delete[] phi;
}

void ImageProcessing::MCF(int T, double tau)
{
    int d = 1;
    int newHeight = height + 2 * d;
    int newWidth = width + 2 * d;
    double eps = 0.01 ;

    for (double** step : steps) {
        for (int i = 0; i < height; i++) {
            delete[] step[i];
        }
        delete[] step;
    }
    steps.clear();

    steps.resize(T + 1);

    for (int t = 0; t <= T; t++)
    {
        steps[t] = new double* [height];
        for (int i = 0; i < height; i++)
        {
            steps[t][i] = new double[width]();
        }
    }

    double** img_old = new double* [newHeight];
    double** img_new = new double* [height];

    for (int i = 0; i < newHeight; i++) {
        img_old[i] = new double[newWidth]();
    }

    for (int i = 0; i < height; i++) {
        img_new[i] = new double[width]();
    }

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            steps[0][i][j] = imageData[i][j];
        }
    }

    mirror(d, imageData, img_old);

    double omega = 1.25;
    double sigmaSOR = 0.0;
    double tolerancia = 1.0E-6;
    int max_iter = 10000;
    double rezid;
    double* b = new double[newHeight * newWidth](); //prava strana
    double* phi = new double[newHeight * newWidth]();

    double uE, uW, uS, uN, uNE, uNW, uSW, uSE, uxE, uxN, uxW, uxS, uyE, uyN, uyW, uyS, uC, uE_grad, uN_grad, uW_grad, uS_grad,u_avg_grad;
    double gE, gN, gW, gS;


    for (int i = 0; i < newHeight; i++)
    {
        for (int j = 0; j < newWidth; j++)
        {
            b[i * newWidth + j] = img_old[i][j];
            phi[i * newWidth + j] = b[i * newWidth + j];
        }
    }


    for (int cas = 1; cas <= T; cas++)
    {
        int iter = 0;
        for (iter = 1; iter < max_iter; iter++)
        {
            rezid = 0.0;
            for (int i = d; i < height + d; i++) {
                for (int j = d; j < width + d; j++) {
                    sigmaSOR = 0.0;

                    uE = img_old[i][j + 1];
                    uN = img_old[i - 1][j];
                    uW = img_old[i][j - 1];
                    uS = img_old[i + 1][j];

                    uNE = img_old[i - 1][j + 1];
                    uNW = img_old[i - 1][j - 1];
                    uSW = img_old[i + 1][j - 1];
                    uSE = img_old[i + 1][j + 1];

                    uC = img_old[i][j];//u current index (stred)

                    uxE = uE - uC;
                    uyE = (uNE - uSE - uS + uN) / 4;

                    uxW = uW - uC;
                    uyW = (uNW + uN - uS - uSW) / 4;

                    uyN = uN - uC;
                    uxN = (uNE + uE - uNW - uW) / 4;

                    uyS = uS - uC;
                    uxS = (uE + uSE - uSW - uW) / 4;

                    uE_grad = uxE * uxE + uyE * uyE;
                    uW_grad = uxW * uxW + uyW * uyW;
                    uN_grad = uxN * uxN + uyN * uyN;
                    uS_grad = uxS * uxS + uyS * uyS;

                    gE = 1.0;
                    gN = 1.0;
                    gW = 1.0;
                    gS = 1.0;

                    //chcem abz to bol grad_epsilon
                    uE_grad = sqrt(uE_grad * uE_grad + eps * eps);
                    uW_grad = sqrt(uW_grad * uW_grad + eps * eps);
                    uS_grad = sqrt(uS_grad * uS_grad + eps * eps);
                    uN_grad = sqrt(uN_grad * uN_grad + eps * eps);

                    u_avg_grad = (uE_grad + uW_grad + uS_grad + uN_grad) / 4;
                    
                    sigmaSOR += -tau *u_avg_grad* ((1/ uE_grad)*gE * phi[i * newWidth + j + 1] + (1 / uN_grad)* gN * phi[(i - 1) * newWidth + j] + 
                        (1 / uS_grad)*gS * phi[(i + 1) * newWidth + j] + (1 / uW_grad)*gW * phi[i * newWidth + j - 1]);

                    int idx = i * newWidth + j;
                    phi[idx] = (1 - omega) * phi[idx] + omega * (b[idx] - sigmaSOR) / (1 + tau* u_avg_grad * ((1 / uE_grad) * gE + (1 / uN_grad) * gN + (1 / uW_grad) * gW + (1 / uS_grad) * gS));
                }
            }

            //zapis vysledku iteracie do 2d
            for (int i = d; i < height + d; i++)
            {
                for (int j = d; j < width + d; j++)
                {
                    img_new[i - d][j - d] = phi[i * newWidth + j];
                }
            }
            mirror(d, img_new, img_old);
            for (int i = 0; i < newHeight; i++)
            {
                for (int j = 0; j < newWidth; j++)
                {
                    phi[i * newWidth + j] = img_old[i][j];
                }
            }

            //rezidua
            for (int i = d; i < height + d; i++) {
                for (int j = d; j < width + d; j++) {
                    uE = img_old[i][j + 1];
                    uN = img_old[i - 1][j];
                    uW = img_old[i][j - 1];
                    uS = img_old[i + 1][j];

                    uNE = img_old[i - 1][j + 1];
                    uNW = img_old[i - 1][j - 1];
                    uSW = img_old[i + 1][j - 1];
                    uSE = img_old[i + 1][j + 1];

                    uC = img_old[i][j];//u current index (stred)

                    uxE = uE - uC;
                    uyE = (uNE - uSE - uS + uN) / 4;

                    uxW = uW - uC;
                    uyW = (uNW + uN - uS - uSW) / 4;

                    uyN = uN - uC;
                    uxN = (uNE + uE - uNW - uW) / 4;

                    uyS = uS - uC;
                    uxS = (uE + uSE - uSW - uW) / 4;

                    uE_grad = uxE * uxE + uyE * uyE;
                    uW_grad = uxW * uxW + uyW * uyW;
                    uN_grad = uxN * uxN + uyN * uyN;
                    uS_grad = uxS * uxS + uyS * uyS;

                    gE = 1.0;
                    gN = 1.0;
                    gW = 1.0;
                    gS = 1.0;

                    //chcem abz to bol grad_epsilon
                    uE_grad = sqrt(uE_grad * uE_grad + eps * eps);
                    uW_grad = sqrt(uW_grad * uW_grad + eps * eps);
                    uS_grad = sqrt(uS_grad * uS_grad + eps * eps);
                    uN_grad = sqrt(uN_grad * uN_grad + eps * eps);

                    u_avg_grad = (uE_grad + uW_grad + uS_grad + uN_grad) / 4;

                    int idx = i * newWidth + j;
                    double Ax = (1 + tau * u_avg_grad * ((1 / uE_grad) * gE + (1 / uN_grad) * gN + (1 / uW_grad) * gW + (1 / uS_grad) * gS)) * phi[idx] - 
                        tau * u_avg_grad*((1 / uE_grad) * gE * phi[i * newWidth + j + 1] + (1 / uN_grad) * gN * phi[(i - 1) * newWidth + j] +
                            (1 / uS_grad) * gS * phi[(i + 1) * newWidth + j] + (1 / uW_grad) * gW * phi[i * newWidth + j - 1]) - b[idx];
                    rezid += Ax * Ax;
                }
            }
            //std::cout << sqrt(rezid) << std::endl;
            if (sqrt(rezid) < tolerancia) break;
        }

        std::cout << "casovy krok: " << cas << " iteracia: " << iter << " reziduum: " << sqrt(rezid) << std::endl;

        for (int i = 0; i < newHeight; i++) {//update pravej strany
            for (int j = 0; j < newWidth; j++) {
                b[i * newWidth + j] = img_old[i][j];
            }
        }

        for (int i = d; i < height + d; i++) {
            for (int j = d; j < width + d; j++) {
                img_new[i - d][j - d] = img_old[i][j];
            }
        }

        for (int i = 0; i < height; i++) {//zapis jednotlivych krokov aby sa dali po jednom vykreslit
            for (int j = 0; j < width; j++) {
                steps[cas][i][j] = img_new[i][j];
            }
        }
    }

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            procData[i][j] = img_new[i][j];

        }
    }

    for (int i = 0; i < newHeight; i++) {
        delete[] img_old[i];
    }
    delete[] img_old;

    for (int i = 0; i < height; i++) {
        delete[] img_new[i];
    }
    delete[] img_new;
    delete[] b;
    delete[] phi;
}

void ImageProcessing::GMCF(int T, double tau, double sigma, double K)
{
    int d = 1;
    int newHeight = height + 2 * d;
    int newWidth = width + 2 * d;
    double eps = 0.01;

    for (double** step : steps) {
        for (int i = 0; i < height; i++) {
            delete[] step[i];
        }
        delete[] step;
    }
    steps.clear();

    steps.resize(T + 1);

    for (int t = 0; t <= T; t++)
    {
        steps[t] = new double* [height];
        for (int i = 0; i < height; i++)
        {
            steps[t][i] = new double[width]();
        }
    }

    double** img_old = new double* [newHeight];
    double** img_new = new double* [height];

    for (int i = 0; i < newHeight; i++) {
        img_old[i] = new double[newWidth]();
    }

    for (int i = 0; i < height; i++) {
        img_new[i] = new double[width]();
    }

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            steps[0][i][j] = imageData[i][j];
        }
    }

    mirror(d, imageData, img_old);

    double omega = 1.25;
    double sigmaSOR = 0.0;
    double tolerancia = 1.0E-6;
    int max_iter = 10000;
    double rezid;
    double* b = new double[newHeight * newWidth](); //prava strana
    double* phi = new double[newHeight * newWidth]();

    double uE, uW, uS, uN, uNE, uNW, uSW, uSE, uxE, uxN, uxW, uxS, uyE, uyN, uyW, uyS, uC, uE_grad, uN_grad, uW_grad, uS_grad, u_avg_grad;
    double gE, gN, gW, gS;


    for (int i = 0; i < newHeight; i++)
    {
        for (int j = 0; j < newWidth; j++)
        {
            b[i * newWidth + j] = img_old[i][j];
            phi[i * newWidth + j] = b[i * newWidth + j];
        }
    }


    for (int cas = 1; cas <= T; cas++)
    {
        int iter = 0;
        for (iter = 1; iter < max_iter; iter++)
        {
            rezid = 0.0;
            for (int i = d; i < height + d; i++) {
                for (int j = d; j < width + d; j++) {
                    sigmaSOR = 0.0;

                    uE = img_old[i][j + 1];
                    uN = img_old[i - 1][j];
                    uW = img_old[i][j - 1];
                    uS = img_old[i + 1][j];

                    uNE = img_old[i - 1][j + 1];
                    uNW = img_old[i - 1][j - 1];
                    uSW = img_old[i + 1][j - 1];
                    uSE = img_old[i + 1][j + 1];

                    uC = img_old[i][j];//u current index (stred)

                    uxE = uE - uC;
                    uyE = (uNE - uSE - uS + uN) / 4;

                    uxW = uW - uC;
                    uyW = (uNW + uN - uS - uSW) / 4;

                    uyN = uN - uC;
                    uxN = (uNE + uE - uNW - uW) / 4;

                    uyS = uS - uC;
                    uxS = (uE + uSE - uSW - uW) / 4;

                    uE_grad = uxE * uxE + uyE * uyE;
                    uW_grad = uxW * uxW + uyW * uyW;
                    uN_grad = uxN * uxN + uyN * uyN;
                    uS_grad = uxS * uxS + uyS * uyS;

                    gE = 1.0 / (1 + K * uE_grad * uE_grad);
                    gN = 1.0 / (1 + K * uN_grad * uN_grad);
                    gW = 1.0 / (1 + K * uW_grad * uW_grad);
                    gS = 1.0 / (1 + K * uS_grad * uS_grad);

                    //chcem abz to bol grad_epsilon
                    uE_grad = sqrt(uE_grad * uE_grad + eps * eps);
                    uW_grad = sqrt(uW_grad * uW_grad + eps * eps);
                    uS_grad = sqrt(uS_grad * uS_grad + eps * eps);
                    uN_grad = sqrt(uN_grad * uN_grad + eps * eps);

                    u_avg_grad = (uE_grad + uW_grad + uS_grad + uN_grad) / 4;

                    sigmaSOR += -tau * u_avg_grad * ((1 / uE_grad) * gE * phi[i * newWidth + j + 1] + (1 / uN_grad) * gN * phi[(i - 1) * newWidth + j] +
                        (1 / uS_grad) * gS * phi[(i + 1) * newWidth + j] + (1 / uW_grad) * gW * phi[i * newWidth + j - 1]);

                    int idx = i * newWidth + j;
                    phi[idx] = (1 - omega) * phi[idx] + omega * (b[idx] - sigmaSOR) / (1 + tau * u_avg_grad * ((1 / uE_grad) * gE + (1 / uN_grad) * gN + (1 / uW_grad) * gW + (1 / uS_grad) * gS));
                }
            }

            //zapis vysledku iteracie do 2d
            for (int i = d; i < height + d; i++)
            {
                for (int j = d; j < width + d; j++)
                {
                    img_new[i - d][j - d] = phi[i * newWidth + j];
                }
            }
            mirror(d, img_new, img_old);
            for (int i = 0; i < newHeight; i++)
            {
                for (int j = 0; j < newWidth; j++)
                {
                    phi[i * newWidth + j] = img_old[i][j];
                }
            }

            //rezidua
            for (int i = d; i < height + d; i++) {
                for (int j = d; j < width + d; j++) {
                    uE = img_old[i][j + 1];
                    uN = img_old[i - 1][j];
                    uW = img_old[i][j - 1];
                    uS = img_old[i + 1][j];

                    uNE = img_old[i - 1][j + 1];
                    uNW = img_old[i - 1][j - 1];
                    uSW = img_old[i + 1][j - 1];
                    uSE = img_old[i + 1][j + 1];

                    uC = img_old[i][j];//u current index (stred)

                    uxE = uE - uC;
                    uyE = (uNE - uSE - uS + uN) / 4;

                    uxW = uW - uC;
                    uyW = (uNW + uN - uS - uSW) / 4;

                    uyN = uN - uC;
                    uxN = (uNE + uE - uNW - uW) / 4;

                    uyS = uS - uC;
                    uxS = (uE + uSE - uSW - uW) / 4;

                    uE_grad = uxE * uxE + uyE * uyE;
                    uW_grad = uxW * uxW + uyW * uyW;
                    uN_grad = uxN * uxN + uyN * uyN;
                    uS_grad = uxS * uxS + uyS * uyS;

                    gE = 1.0/(1+K*uE_grad*uE_grad);
                    gN = 1.0 / (1 + K * uN_grad * uN_grad);
                    gW = 1.0 / (1 + K * uW_grad * uW_grad);
                    gS = 1.0 / (1 + K * uS_grad * uS_grad);

                    //chcem abz to bol grad_epsilon
                    uE_grad = sqrt(uE_grad * uE_grad + eps * eps);
                    uW_grad = sqrt(uW_grad * uW_grad + eps * eps);
                    uS_grad = sqrt(uS_grad * uS_grad + eps * eps);
                    uN_grad = sqrt(uN_grad * uN_grad + eps * eps);

                    u_avg_grad = (uE_grad + uW_grad + uS_grad + uN_grad) / 4;

                    int idx = i * newWidth + j;
                    double Ax = (1 + tau * u_avg_grad * ((1 / uE_grad) * gE + (1 / uN_grad) * gN + (1 / uW_grad) * gW + (1 / uS_grad) * gS)) * phi[idx] -
                        tau * u_avg_grad * ((1 / uE_grad) * gE * phi[i * newWidth + j + 1] + (1 / uN_grad) * gN * phi[(i - 1) * newWidth + j] +
                            (1 / uS_grad) * gS * phi[(i + 1) * newWidth + j] + (1 / uW_grad) * gW * phi[i * newWidth + j - 1]) - b[idx];
                    rezid += Ax * Ax;
                }
            }
            //std::cout << sqrt(rezid) << std::endl;
            if (sqrt(rezid) < tolerancia) break;
        }

        std::cout << "casovy krok: " << cas << " iteracia: " << iter << " reziduum: " << sqrt(rezid) << std::endl;

        for (int i = 0; i < newHeight; i++) {//update pravej strany
            for (int j = 0; j < newWidth; j++) {
                b[i * newWidth + j] = img_old[i][j];
            }
        }

        for (int i = d; i < height + d; i++) {
            for (int j = d; j < width + d; j++) {
                img_new[i - d][j - d] = img_old[i][j];
            }
        }

        for (int i = 0; i < height; i++) {//zapis jednotlivych krokov aby sa dali po jednom vykreslit
            for (int j = 0; j < width; j++) {
                steps[cas][i][j] = img_new[i][j];
            }
        }
    }

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            procData[i][j] = img_new[i][j];

        }
    }

    for (int i = 0; i < newHeight; i++) {
        delete[] img_old[i];
    }
    delete[] img_old;

    for (int i = 0; i < height; i++) {
        delete[] img_new[i];
    }
    delete[] img_new;
    delete[] b;
    delete[] phi;
}

void ImageProcessing::vzdialenostna()
{
    int d = 1;
    int newHeight = height + 2 * d;
    int newWidth = width + 2 * d;

    double** mirrored = new double* [newHeight];
    
    bool** F = new bool* [newHeight];   //true bude ze je tam uz vypocitene, false ze este nie
    double** vzdialenost = new double*[newHeight];    //toto budu vzdialenostne hodnoty
    double** vzdialenost_tmp = new double* [height];

    double hore, dole, vpravo, vlavo,aktual;
    double stvorlistok = 50.0 / 255.0;  //farba, ktorou je zafarbeny stvorlistok

    for (int i = 0; i < newHeight; i++) {
        mirrored[i] = new double[newWidth]();
        vzdialenost[i] = new double[newWidth]();
        F[i] = new bool[newWidth];
        for (int j = 0; j < newWidth; j++)
        {
            F[i][j] = false;    //defaultna hodnota je ze este neviem ci je to hrana alebo nie takze vsetko beriem ako nie
        }
    }
    for (int i = 0; i < newHeight; i++) {
        vzdialenost_tmp[i] = new double[width]();
    }


    mirror(d, imageData, mirrored);

    //inicializacia-> najdenie hrany a zapisanie hodnoty vzdialenosti na hrane (ulozenie ze to in fact je hrana a ze je tam vzdialenost nulova)
    for (int i = d; i < height+d; i++)
    {
        for (int j = d; j < width+d; j++)
        {
            if (mirrored[i][j] == stvorlistok)
            {
                hore = mirrored[i - 1][j];
                dole = mirrored[i + 1][j];
                vpravo = mirrored[i][j + 1];
                vlavo = mirrored[i][j - 1];
                aktual = mirrored[i][j];

                if (aktual != hore || aktual != dole || aktual != vpravo || aktual != vlavo)
                {
                    F[i][j] = true; //zapisujem ze ta hodnota je uz najdena vypocitana viem vzdialenost, lebo na hrane viem vzdialenost = je nulova
                    vzdialenost[i][j] = 0;  //a aj si tu vzdialenost zapisem, je nula
                }
            }            
        }
    }

    //tu budem pokracovat tym realnym hladanim vzdialenosti

    bool end_switch;
    double M_vlavo, M_vpravo, M_hore, M_dole, d_nove,d_stare,tolerancia=0.001;
    double tau = 0.4;
    double minimum=DBL_MAX, maximum=-DBL_MAX;
    int maxIter = 1000;
    int iter = 1;

    do
    {
        iter++;
        end_switch = true;  //defaultne zapnem, ak najdem nejaký ktorý ešte trepa vypoèíta, tak vypnem

        for (int i = d; i < height + d; i++)
        {
            for (int j = d; j < width + d; j++)
            {
                if (!F[i][j]) //iba ak ešte nie je vypoèitané
                {
                    end_switch = false; //vypnem konèovaciu podmienku
                    M_vlavo = std::min(vzdialenost[i][j-1] - vzdialenost[i][j], 0.0);
                    M_vlavo = M_vlavo * M_vlavo;

                    M_vpravo = std::min(vzdialenost[i][j+1] - vzdialenost[i][j], 0.0);
                    M_vpravo = M_vpravo * M_vpravo;

                    M_hore = std::min(vzdialenost[i-1][j] - vzdialenost[i][j], 0.0);
                    M_hore = M_hore * M_hore;

                    M_dole = std::min(vzdialenost[i + 1][j] - vzdialenost[i][j], 0.0);
                    M_dole = M_dole * M_dole;

                    d_stare = vzdialenost[i][j];
                    d_nove = d_stare + tau - tau * sqrt(std::max(M_vlavo, M_vpravo) + std::max(M_dole, M_hore));
                    vzdialenost_tmp[i-d][j-d] = d_nove;

                    if (fabs(d_stare - d_nove) < tolerancia)
                    {
                        F[i][j] = true;
                        
                        if (d_nove < minimum)minimum = d_nove;
                        if (d_nove > maximum)maximum = d_nove;
                    }
                } 
            }
        }
        mirror(d, vzdialenost_tmp, vzdialenost);

    } while (!end_switch && iter<maxIter);

    double rozdiel = maximum - minimum;


    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            procData[i][j] = (vzdialenost_tmp[i][j]-minimum) / rozdiel;
        }
    }


}

void ImageProcessing::explicitna_sigma(double sigma, double** vstup, double** vysledok)
{
    //funkcia na jeden krok explicitnej schemy s krokom sigma, ako predpriprava do perona malika semi implicitnej
    int d = 1;
    int newHeight = height + 2 * d;
    int newWidth = width + 2 * d;

    double** img_old = new double* [newHeight];
    double** img_new = new double* [height];

    for (int i = 0; i < newHeight; i++) {
        img_old[i] = new double[newWidth]();
    }

    for (int i = 0; i < height; i++) {
        img_new[i] = new double[width]();
    }

    double sum = 0.0;
    mirror(d, vstup, img_old);

    //std::cout << "Step 0" << " average: " << sum / (height * width) << std::endl;

    sum = 0.0;
    for (int i = d; i < height + d; i++)
    {
        for (int j = d; j < width + d; j++)
        {
            img_new[i - d][j - d] = (1 - 4 * sigma) * img_old[i][j] + sigma * (img_old[i - 1][j] + img_old[i + 1][j] + img_old[i][j - 1] + img_old[i][j + 1]);
            sum += img_new[i - d][j - d];
        }
    }
   // std::cout << "Step " << sigma << " average: " << sum / (height * width) << std::endl;

    mirror(d, img_new, img_old);

    for (int i = 0; i < newHeight; i++) {
        for (int j = 0; j < newWidth; j++) {
            vysledok[i][j] = img_old[i][j];
        }
    }

    for (int i = 0; i < newHeight; i++) {
        delete[] img_old[i];
    }
    delete[] img_old;

    for (int i = 0; i < height; i++) {
        delete[] img_new[i];
    }
    delete[] img_new;
}

void ImageProcessing::hranovy()
{
    int d = 1;
    double K = 200;
    int newWidth = width + 2 * d;
    int newHeight = height + 2 * d;
    double** u = new double* [newHeight];
    double uE, uW, uS, uN, uNE, uNW, uSW, uSE, uxE, uxN, uxW, uxS, uyE, uyN, uyW, uyS, uC,uE_grad,uN_grad,uW_grad,uS_grad,grad_priemer,g;

    for (int i = 0; i < newHeight; i++) {
        u[i] = new double[newWidth]();
    }
    
    implicitna_sigma(1, imageData, u);    //linearna filtracia pred hranovou detekciou
    //u uz je vacsej velkosti, teda ma hrany zmirrorovane a mozem pristupovat k susedom

    for (int i = d; i < height + d; i++) 
    {
        for (int j = d; j < width + d; j++) 
        {
            uE = u[i][j + 1];
            uN = u[i - 1][j];
            uW = u[i][j - 1];
            uS = u[i + 1][j];

            uNE = u[i - 1][j + 1];
            uNW = u[i - 1][j - 1];
            uSW = u[i + 1][j - 1];
            uSE = u[i + 1][j + 1];

            uC = u[i][j];//u current index (stred)

            uxE = uE - uC;
            uyE = (uNE - uSE - uS + uN) / 4;

            uxW = uW - uC;
            uyW = (uNW + uN - uS - uSW) / 4;

            uyN = uN - uC;
            uxN = (uNE + uE - uNW - uW) / 4;

            uyS = uS - uC;
            uxS = (uE + uSE - uSW - uW) / 4;

            uE_grad = uxE * uxE + uyE * uyE;
            uW_grad = uxW * uxW + uyW * uyW;
            uN_grad = uxN * uxN + uyN * uyN;
            uS_grad = uxS * uxS + uyS * uyS;

            grad_priemer = (uE_grad + uN_grad + uW_grad + uS_grad) / 4;

            g = 1 / (1 + K * grad_priemer);

            procData[i - d][j - d] = g;
        }
    }

}

double ImageProcessing::find_min()
{
    for (int i = 0; i < 255; i++)
    {
        if (histogram[i] != 0)
        {
            return i / 255.0;
        }
    }
}

double ImageProcessing::find_max()
{
    for (int i = 255; i >=0; i--)
    {
        if (histogram[i] != 0)
        {
            return i / 255.0;
        }
    }
}

double ImageProcessing::mean(double** img)
{
    double sum = 0.0;
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            sum += img[i][j];
        }
    }
    return sum / (height * width);
}


void ImageProcessing::mirror(int d, double** img_data, double** mirrored)
{
    int newHeight = height + 2 * d;
    int newWidth = width + 2 * d;


    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            mirrored[i + d][j + d] = img_data[i][j];
        }
    }

    for (int i = d; i < height + d; i++)
    {
        for (int j = 0; j < d; j++)
        {
            mirrored[i][j] = mirrored[i][2 * d - j -1];
        }

        for (int j = 0; j < d; j++)
        {
            mirrored[i][newWidth - 1 - j] = mirrored[i][newWidth  - 2 * d + j];
        }
    }

    for (int j = 0; j < newWidth; j++)
    {
        for (int i = 0; i < d; i++)
        {
            mirrored[i][j] = mirrored[2 * d - i -1][j];
        }

        for (int i = 0; i < d; i++)
        {
            mirrored[newHeight - 1 - i][j] = mirrored[newHeight  - 2 * d + i][j];
        }
    }

    //Zapis do suboru na kontrolu
//std::ofstream outFile("test.pgm", std::ios::binary);
//
//if (!outFile) {
//    std::cerr << "Error opening file for writing.\n";
//    return;
//}
//
//
//outFile << "P2\n" << newWidth << " " << newHeight << "\n255\n";
//
//for (int i = 0; i < newHeight; ++i) {
//    for (int j = 0; j < newWidth; ++j) {
//        outFile << (int)(mirrored[i][j] * 255 +0.5) << " ";
//    }
//    outFile << "\n";
//}
//
//outFile.close();

}


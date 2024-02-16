#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QStandardPaths>
#include <QFileDialog>
#include <QLabel>

#include <QApplication>
#include <QImage>
#include <QLabel>
#include <QPainter>
#include <QTextStream>
#include <QString>

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

const int INTENS_MIN = 0;
const int INTENS_MAX = 255;

const double EPS = 1.0e-14;

bool circle = false;
bool line = false;

void calcHisto(QImage &image, double histo[])
{
    for (int i = 0; i <= INTENS_MAX; i++)
    {
        histo[i] = 0;
    }
    for (int indx_row = 0; indx_row < image.height(); indx_row++)
    {
        quint8* ptr_row = (quint8*)(image.bits()
                                      + indx_row * image.bytesPerLine());
        for (int indx_col = 0; indx_col < image.width(); indx_col++)
        {
            histo[ptr_row[indx_col]]++;
        }
    }
    int numb_pix = image.height() * image.width();
    for (int i = 0; i <= INTENS_MAX; i++)
    {
        histo[i] /= numb_pix;
    }
}//calcHisto

void thresh(QImage &image, int thr)
{
    for (int indx_row = 0; indx_row < image.height(); indx_row++)
    {
        quint8* ptr_row = (quint8*)(image.bits()
                                      + indx_row * image.bytesPerLine());
        for (int indx_col = 0; indx_col < image.width(); indx_col++)
        {
            ptr_row[indx_col] =
                (ptr_row[indx_col] < thr) ? INTENS_MIN : INTENS_MAX;
        }
    }
}//thresh

int otsu(const double histo[])
{
    // compute cumulative sums
    double p_1[INTENS_MAX + 1] = {0};
    p_1[0] = histo[0];
    for (int i = 1; i <= INTENS_MAX; i++)
    {
        p_1[i] = p_1[i - 1] + histo[i];
    }

    // cumulative mean
    double m[INTENS_MAX + 1] = {0};
    for (int i = 1; i <= INTENS_MAX; i++)
    {
        m[i] = m[i - 1] + i * histo[i];
    }

    // global mean
    double m_g = m[INTENS_MAX];

    // between-class
    double b_c[INTENS_MAX + 1] = {0};
    for (int i = 1; i <= INTENS_MAX; i++)
    {
        double div = (p_1[i] * (1 - p_1[i]));
        b_c[i] =
            fabs(div < EPS) ? 0 :
                ((m_g * p_1[i] - m[i]) * (m_g * p_1[i] - m[i])) / div;
    }

    // find max
    double max = 0;
    int max_i = 0;
    for (int i = 0; i <= INTENS_MAX; i++)
    {
        if (b_c[i] > max)
        {
            max = b_c[i];
            max_i = i;
        }
    }
    QTextStream(stdout) << "Th: " << max_i << Qt::endl;

    return max_i;
}

inline double degToRad(double deg)
{
    return (M_PI / 180.0) * deg;
}//degToRad

void houghTransform(QImage &image, int numb_rho, int numb_theta)
{
    const double MIN_THETA = -90;
    const double MAX_THETA = 90;

    int diag = sqrt(image.height() * image.height() + image.width() * image.width());

    double delta_rho = 1.0 * (2 * (diag)) / numb_rho;
    double delta_theta = 180.0 / numb_theta;

    // construct accumulator space
    vector<vector<int>> accum(numb_rho, vector<int>(numb_theta));
    for (int indx_row = 0; indx_row < numb_rho; indx_row++)
    {
        for (int indx_col = 0; indx_col < numb_theta; indx_col++)
        {
            accum[indx_row][indx_col] = 0;
        }
    }

    // for each edge (object) pixel in the image
    for (int indx_row = 0; indx_row < image.height(); indx_row++)
    {
        quint8* ptr_row = (quint8*)(image.bits()
                                      + indx_row * image.bytesPerLine());
        for (int indx_col = 0; indx_col < image.width(); indx_col++)
        {
            if (ptr_row[indx_col] == INTENS_MIN)
            {
                // voting procedure
                for (int indx_theta = 0; indx_theta < numb_theta; indx_theta++)
                {
                    double theta = MIN_THETA + (indx_theta * delta_theta);
                    double rho = indx_row * cos(degToRad(theta)) + indx_col * sin(degToRad(theta));
                    int indx_rho = (int)((diag) + rho) / delta_rho;
                    accum[indx_rho][indx_theta]++;
                }
            }
        }
    }
}// Hough Transfrom for Straight Lines


void houghTransformCircles(QImage &image, int minRadius, int maxRadius) {
    int width = image.width();
    int height = image.height();
    vector<vector<vector<int>>> accum(width, vector<vector<int>>(height, vector<int>(maxRadius - minRadius + 1, 0)));
    int maxVotes = 0;

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            quint8 pixel = image.pixelColor(x, y).red(); // Assuming image is pre-processed for edges
            if (pixel == INTENS_MAX)
            {
                for (int r = minRadius; r <= maxRadius; ++r)
                {
                    for (int theta = 0; theta < 360; ++theta)
                    {
                        int a = int(x - r * cos(degToRad(theta)));
                        int b = int(y - r * sin(degToRad(theta)));
                        if (a >= 0 && a < width && b >= 0 && b < height)
                        {
                            accum[a][b][r - minRadius]++;
                            maxVotes = max(maxVotes, accum[a][b][r - minRadius]);
                        }
                    }
                }
            }
        }
    }
} // Hough Transform for Circles

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_clicked()
{
    const QString downloadsFolder = QStandardPaths::writableLocation(QStandardPaths::DownloadLocation);
    QString file_name = QFileDialog::getOpenFileName(this, "Upload an image.", downloadsFolder);

    QImage image;

    if (image.load(file_name))
    {
        QTextStream(stdout) << "Image loaded: " << file_name << Qt::endl;
        QTextStream(stdout) << "Format: " << image.format() << Qt::endl;

        if (image.format() != QImage::Format_Grayscale8)
        {
            image = image.convertToFormat(QImage::Format_Grayscale8);
            QTextStream(stdout) << "Format not grayscale, converted to grayscale: " << image.format() << Qt::endl;
        }

        double histo[INTENS_MAX + 1];
        calcHisto(image, histo);

        // Otsu's method to determine global threshold
        int th = otsu(histo);

        // transform image to binary
        thresh(image, th);

        if(line){
            // Hough Transform for Straight Lines
            houghTransform(image, 180, 180);
            line = false;
        }
        else if (circle)
        {
            // Hough Transform for Circles
            houghTransformCircles(image, 20, 100);
            circle = false;
        }

        QLabel *label = new QLabel(this);
        label->setWindowFlags(Qt::Window);
        label ->setPixmap(QPixmap::fromImage(image));
        label->show();
    }
    else
    {
        QTextStream(stdout) << "Cannot load image: " << file_name << Qt::endl;
    }

}


void MainWindow::on_circle_toggled(bool checked)
{
    circle = checked;
}


void MainWindow::on_line_toggled(bool checked)
{
    line = checked;
}


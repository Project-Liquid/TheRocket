#include "transducer.h"

Transducer::Transducer(int pin, float P_MIN, float P_MAX) {
  this->pin = pin;
  this->V_REF = V_REF;
  this->P_MIN = P_MIN;
  this->P_MAX = P_MAX;
  I_MAX = 20;
  _computeSGCoeffs();
}

static float Transducer::barToPSI(float bar) {
  return bar * 14.504;
}

static float Transducer::PSIToBar(float psi) {
  return psi / 14.504;
}

// Call _computeSGCoeffs() in your constructor or begin() method.

void Transducer::_computeSGCoeffs() {
    const int m = SG_WINDOW;
    const int half = m / 2;
    const int poly = SG_POLY;

    // Build Vandermonde matrix A (m x (poly+1)), compute (A^T A)^-1 A^T
    // For the center point, SG coefficients are the middle row of the hat matrix.
    // We compute this numerically using the normal equations.
    const int p = poly + 1;
    double ATA[4][4] = {};     // max poly+1 = 4
    double ATe[4] = {};        // e = unit vector at center row

    for (int i = 0; i < m; i++) {
        double x = i - half;
        double row[4];
        row[0] = 1.0;
        for (int j = 1; j < p; j++) row[j] = row[j-1] * x;

        for (int r = 0; r < p; r++)
            for (int c = 0; c < p; c++)
                ATA[r][c] += row[r] * row[c];

        // Center point: i == half
        if (i == half)
            for (int r = 0; r < p; r++)
                ATe[r] = row[r];
    }

    // Solve ATA * beta = ATe  (Gaussian elimination)
    double mat[4][5] = {};
    for (int r = 0; r < p; r++) {
        for (int c = 0; c < p; c++) mat[r][c] = ATA[r][c];
        mat[r][p] = ATe[r];
    }
    for (int col = 0; col < p; col++) {
        // Pivot
        int pivot = col;
        for (int r = col+1; r < p; r++)
            if (fabs(mat[r][col]) > fabs(mat[pivot][col])) pivot = r;
        for (int c = 0; c <= p; c++) { double t = mat[col][c]; mat[col][c] = mat[pivot][c]; mat[pivot][c] = t; }
        for (int r = 0; r < p; r++) {
            if (r == col) continue;
            double f = mat[r][col] / mat[col][col];
            for (int c = col; c <= p; c++) mat[r][c] -= f * mat[col][c];
        }
    }
    double beta[4];
    for (int r = 0; r < p; r++) beta[r] = mat[r][p] / mat[r][r];

    // _sgCoeffs[i] = dot(beta, row_i)  — these are the convolution weights
    for (int i = 0; i < m; i++) {
        double x = i - half;
        double val = 0.0, xi = 1.0;
        for (int j = 0; j < p; j++, xi *= x) val += beta[j] * xi;
        _sgCoeffs[i] = (float)val;
    }
}

float Transducer::readPressure() {
    raw = analogRead(pin);
    v = raw * V_REF / 1023.0;
    i_mA = (v / R_SHUNT) * 1000.0;
    if (i_mA < 0) {
        i_mA = 0;
    } else if (i_mA > 25) {
        i_mA = 25;
    }
    pressure = (v - 1.0) * (P_MAX - P_MIN) / (V_MAX - V_MIN) + P_MIN;

    // --- Savitzky-Golay filter ---
    // _pressureBuffer[_bufferIndex] = pressure;
    // _bufferIndex = (_bufferIndex + 1) % SG_WINDOW;
    // if (_bufferCount < SG_WINDOW) _bufferCount++;

    // if (_bufferCount < SG_WINDOW) return pressure; // buffer not full yet

    // // Apply SG convolution; oldest sample is at _bufferIndex (circular buffer)
    // float filtered = 0.0f;
    // for (int i = 0; i < SG_WINDOW; i++) {
    //     int idx = (_bufferIndex + i) % SG_WINDOW;
    //     filtered += _sgCoeffs[i] * _pressureBuffer[idx];
    // }

    // Return the estimate at window center (window_length / 2 samples ago)
    //pressure = filtered;
    return pressure;
}

String Transducer::status() {
  readPressure();
  String output = "";
  output += "raw = " + String(raw) + " V = " + String(v, 3) + "V  I = " + String(i_mA, 2) + " mA  P = " + String(pressure, 3) + " PSI";
  return output;
}

String Transducer::value() {
  return String(readPressure(), 3);
}


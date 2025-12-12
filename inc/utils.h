#ifndef UTILS_H
#define UTILS_H

#include "core.h"
#include <vector>
#include <string>

// Гамма-коррекция и тональное отображение
std::vector<unsigned char> tonemapAndGammaCorrect(
    const std::vector<Color> &hdrImage,
    float exposure = 1.0f,
    float gamma = 2.2f,
    bool autoExposure = true);

// Сохранение в PPM формат
void savePPM(const std::string &filename, const std::vector<unsigned char> &image,
             int width, int height);

#endif // UTILS_H
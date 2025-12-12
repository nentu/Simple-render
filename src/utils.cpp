#include "inc/utils.h"

// Гамма-коррекция и тональное отображение
std::vector<unsigned char> tonemapAndGammaCorrect(
    const std::vector<Color> &hdrImage,
    float exposure = 1.0f,
    float gamma = 2.2f,
    bool autoExposure = true)
{

    std::vector<unsigned char> ldrImage(hdrImage.size() * 3);

    // Вычисление средней яркости для автоэкспозиции
    float avgLuminance = 0;
    if (autoExposure)
    {
        for (const auto &c : hdrImage)
        {
            float lum = 0.299f * c.x + 0.587f * c.y + 0.114f * c.z;
            avgLuminance += std::log(0.0001f + lum);
        }
        avgLuminance = std::exp(avgLuminance / hdrImage.size());
        exposure = 0.18f / avgLuminance; // Нормализация к средне-серому 0.18
        std::cout << "Автоэкспозиция: " << exposure << std::endl;
    }

    // Применение тонального отображения и гамма-коррекции
    for (size_t i = 0; i < hdrImage.size(); ++i)
    {
        Color c = hdrImage[i] * exposure;

        // Простое отсечение (clipping)
        c.x = std::min(c.x, 1.0f);
        c.y = std::min(c.y, 1.0f);
        c.z = std::min(c.z, 1.0f);

        // Гамма-коррекция
        float invGamma = 1.0f / gamma;
        c.x = std::pow(std::max(0.0f, c.x), invGamma);
        c.y = std::pow(std::max(0.0f, c.y), invGamma);
        c.z = std::pow(std::max(0.0f, c.z), invGamma);

        // Преобразование в 8-бит
        ldrImage[i * 3 + 0] = static_cast<unsigned char>(std::min(255.0f, c.x * 255.0f));
        ldrImage[i * 3 + 1] = static_cast<unsigned char>(std::min(255.0f, c.y * 255.0f));
        ldrImage[i * 3 + 2] = static_cast<unsigned char>(std::min(255.0f, c.z * 255.0f));
    }

    return ldrImage;
}

// Сохранение в PPM формат
void savePPM(const std::string &filename, const std::vector<unsigned char> &image,
             int width, int height)
{
    std::ofstream file(filename, std::ios::binary);
    if (!file)
    {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << std::endl;
        return;
    }

    file << "P6\n"
         << width << " " << height << "\n255\n";
    file.write(reinterpret_cast<const char *>(image.data()), image.size());

    std::cout << "Изображение сохранено: " << filename << std::endl;
}

#include "inc/pathtracer.h"

class PathTracer
{
private:
    Scene &scene;
    Camera &camera;
    int maxDepth;
    int samplesPerPixel;

public:
    PathTracer(Scene &s, Camera &c, int depth = 10, int spp = 64)
        : scene(s), camera(c), maxDepth(depth), samplesPerPixel(spp) {}

    void setSamplesPerPixel(int spp) { samplesPerPixel = spp; }
    void setMaxDepth(int depth) { maxDepth = depth; }

    // Трассировка одного пути
    Color tracePath(float ox, float oy, float oz,
                    float dx, float dy, float dz,
                    RandomGenerator &rng) const
    {
        Color throughput(1, 1, 1); // Вес пути
        Color radiance(0, 0, 0);   // Накопленное излучение

        float rayOx = ox, rayOy = oy, rayOz = oz;
        float rayDx = dx, rayDy = dy, rayDz = dz;

        for (int depth = 0; depth < maxDepth; ++depth)
        {
            float t;
            int triId;

            if (!scene.intersect(rayOx, rayOy, rayOz, rayDx, rayDy, rayDz, t, triId))
            {
                // Луч ушел в пустоту - фоновый цвет
                break;
            }

            const Triangle &tri = scene.getTriangle(triId);
            const Material &mat = scene.getMaterial(tri.materialId);

            // Точка пересечения
            Vec3 hitPoint(rayOx + rayDx * t, rayOy + rayDy * t, rayOz + rayDz * t);
            Vec3 normal = tri.normal;
            Vec3 rayDir(rayDx, rayDy, rayDz);

            // Если нормаль направлена от луча, переворачиваем
            if (normal.dot(rayDir) > 0)
            {
                normal = -normal;
            }

            // Добавляем эмиссию, если это источник света
            if (mat.type == MaterialType::EMISSIVE)
            {
                radiance += throughput * mat.emission;
                break; // Попали в источник света
            }

            // Прямое освещение (NEE - Next Event Estimation)
            Color directLight = sampleDirectLight(hitPoint, normal, mat, rng);
            radiance += throughput * directLight;

            // Выбор события: диффузное или зеркальное отражение
            bool isSpecular = false;
            Color brdfColor;

            if (mat.type == MaterialType::SPECULAR)
            {
                isSpecular = true;
                brdfColor = mat.specularColor;
            }
            else if (mat.type == MaterialType::MIXED)
            {
                // Русская рулетка для выбора типа отражения
                if (rng.uniform() < mat.specularProb)
                {
                    isSpecular = true;
                    brdfColor = mat.specularColor / mat.specularProb;
                }
                else
                {
                    brdfColor = mat.diffuseColor / (1.0f - mat.specularProb);
                }
            }
            else
            {
                brdfColor = mat.diffuseColor;
            }

            // Генерация нового направления
            Vec3 newDir;
            float pdf;

            if (isSpecular)
            {
                // Идеальное зеркальное отражение
                newDir = rayDir - normal * (2.0f * rayDir.dot(normal));
                pdf = 1.0f; // Дельта-функция
            }
            else
            {
                // Косинусно-взвешенная выборка для диффузного отражения
                newDir = rng.cosineHemisphere(normal);
                float cosTheta = newDir.dot(normal);
                pdf = cosTheta / M_PI; // PDF косинусной выборки

                // BRDF Ламберта: f = albedo / pi
                // Для косинусной выборки: throughput *= f * cos / pdf = f * cos / (cos/pi) = albedo
                // Уже учтено в brdfColor
            }

            // Обновляем throughput
            if (!isSpecular)
            {
                // Для диффузного: BRDF = color / pi, pdf = cos / pi
                // throughput *= (color / pi) * cos / (cos / pi) = color
                throughput = throughput * brdfColor;
            }
            else
            {
                throughput = throughput * brdfColor;
            }

            // Русская рулетка для завершения пути
            if (depth > 3)
            {
                float survivalProb = std::min(throughput.max_component(), 0.95f);
                if (rng.uniform() > survivalProb)
                {
                    break;
                }
                throughput = throughput / survivalProb;
            }

            // Устанавливаем новый луч
            rayOx = hitPoint.x;
            rayOy = hitPoint.y;
            rayOz = hitPoint.z;
            rayDx = newDir.x;
            rayDy = newDir.y;
            rayDz = newDir.z;
        }

        return radiance;
    }

    // Выборка прямого освещения (Next Event Estimation)
    Color sampleDirectLight(const Vec3 &point, const Vec3 &normal,
                            const Material &mat, RandomGenerator &rng) const
    {
        const auto &lights = scene.getLights();
        if (lights.empty() || mat.type == MaterialType::SPECULAR)
        {
            return Color(0, 0, 0);
        }

        // Выбираем случайный источник с вероятностью пропорциональной мощности
        float lightSample = rng.uniform() * scene.getTotalLightPower();
        float cumPower = 0;
        const AreaLight *selectedLight = nullptr;
        float lightPdf = 0;

        for (const auto &light : lights)
        {
            cumPower += light.power;
            if (cumPower >= lightSample)
            {
                selectedLight = &light;
                lightPdf = light.power / scene.getTotalLightPower();
                break;
            }
        }

        if (!selectedLight)
        {
            selectedLight = &lights.back();
            lightPdf = selectedLight->power / scene.getTotalLightPower();
        }

        // Случайная точка на выбранном источнике
        Vec3 lightPoint = selectedLight->triangle.samplePoint(rng);
        Vec3 toLight = lightPoint - point;
        float distSq = toLight.lengthSquared();
        float dist = std::sqrt(distSq);
        toLight = toLight / dist;

        // Проверка геометрии
        float cosAtSurface = normal.dot(toLight);
        if (cosAtSurface <= 0)
            return Color(0, 0, 0);

        float cosAtLight = -selectedLight->triangle.normal.dot(toLight);
        if (cosAtLight <= 0)
            return Color(0, 0, 0);

        // Проверка видимости
        if (!scene.isVisible(point, lightPoint))
        {
            return Color(0, 0, 0);
        }

        // Вычисление освещенности
        // L_direct = Le * BRDF * G / pdf
        // G = cos(theta_surface) * cos(theta_light) / r^2
        // pdf = 1 / area * (light selection pdf)

        float G = cosAtSurface * cosAtLight / distSq;
        float areaPdf = 1.0f / selectedLight->triangle.area;
        float pdf = lightPdf * areaPdf;

        Color brdf;
        if (mat.type == MaterialType::MIXED)
        {
            brdf = mat.diffuseColor / M_PI; // Используем только диффузную часть для NEE
        }
        else
        {
            brdf = mat.diffuseColor / M_PI;
        }

        return selectedLight->emission * brdf * G / pdf;
    }

    // Рендеринг изображения
    std::vector<Color> render()
    {
        std::vector<Color> image(camera.width * camera.height);

        auto startTime = std::chrono::high_resolution_clock::now();
        int totalPixels = camera.width * camera.height;
        int completed = 0;

        std::cout << "Рендеринг: " << camera.width << "x" << camera.height
                  << ", " << samplesPerPixel << " лучей на пиксель" << std::endl;

#pragma omp parallel for schedule(dynamic, 16)
        for (int y = 0; y < camera.height; ++y)
        {
            RandomGenerator rng(y * 1234567);

            for (int x = 0; x < camera.width; ++x)
            {
                Color pixelColor(0, 0, 0);

                for (int s = 0; s < samplesPerPixel; ++s)
                {
                    float ox, oy, oz, dx, dy, dz;
                    camera.generateRay(x, y, ox, oy, oz, dx, dy, dz, rng);
                    pixelColor += tracePath(ox, oy, oz, dx, dy, dz, rng);
                }

                image[y * camera.width + x] = pixelColor / samplesPerPixel;
            }

#pragma omp critical
            {
                completed += camera.width;
                if (completed % (totalPixels / 10) < camera.width)
                {
                    int progress = (completed * 100) / totalPixels;
                    std::cout << "  Прогресс: " << progress << "%" << std::endl;
                }
            }
        }

        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
        std::cout << "Время рендеринга: " << duration.count() / 1000.0 << " сек" << std::endl;

        return image;
    }
};

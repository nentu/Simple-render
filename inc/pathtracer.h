#ifndef PATHTRACER_H
#define PATHTRACER_H

#include "core.h"
#include <vector>

class PathTracer
{
private:
    Scene &scene;
    Camera &camera;
    int maxDepth;
    int samplesPerPixel;

    // Выборка прямого освещения (Next Event Estimation)
    Color sampleDirectLight(const Vec3 &point, const Vec3 &normal,
                            const Material &mat, RandomGenerator &rng) const;

public:
    PathTracer(Scene &s, Camera &c, int depth = 10, int spp = 64);

    void setSamplesPerPixel(int spp);
    void setMaxDepth(int depth);

    // Трассировка одного пути
    Color tracePath(float ox, float oy, float oz,
                    float dx, float dy, float dz,
                    RandomGenerator &rng) const;

    // Рендеринг изображения
    std::vector<Color> render();
};

#endif // PATHTRACER_H
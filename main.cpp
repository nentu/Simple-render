/**
 * Лабораторная работа №5: Трассировка путей (Path Tracing)
 *
 * Цель: Освоить простейшие методы синтеза изображений трехмерных сцен
 * с учетом глобального освещения методом трассировки путей.
 *
 * Реализовано:
 * - Трассировка путей с глобальным освещением
 * - Геометрия сцены на основе треугольных сеток (Embree API)
 * - Материалы: диффузное (Ламберт) и зеркальное отражение
 * - Выборка по значимости и русская рулетка
 * - Протяженные источники света с диаграммой Ламберта
 * - Точечная камера с антиалиасингом
 * - HDR -> LDR преобразование с гамма-коррекцией
 * - Вывод в формате PPM
 */

#include <embree4/rtcore.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <algorithm>
#include <limits>
#include <chrono>
#include <string>
#include <sstream>

// ==================== Математические структуры ====================

struct Vec3
{
    float x, y, z;

    Vec3() : x(0), y(0), z(0) {}
    Vec3(float x, float y, float z) : x(x), y(y), z(z) {}

    Vec3 operator+(const Vec3 &v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3 operator-(const Vec3 &v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    Vec3 operator*(float s) const { return Vec3(x * s, y * s, z * s); }
    Vec3 operator*(const Vec3 &v) const { return Vec3(x * v.x, y * v.y, z * v.z); }
    Vec3 operator/(float s) const { return Vec3(x / s, y / s, z / s); }
    Vec3 operator-() const { return Vec3(-x, -y, -z); }
    Vec3 &operator+=(const Vec3 &v)
    {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }
    Vec3 &operator*=(float s)
    {
        x *= s;
        y *= s;
        z *= s;
        return *this;
    }

    float dot(const Vec3 &v) const { return x * v.x + y * v.y + z * v.z; }

    Vec3 cross(const Vec3 &v) const
    {
        return Vec3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }

    float length() const { return std::sqrt(x * x + y * y + z * z); }
    float lengthSquared() const { return x * x + y * y + z * z; }

    Vec3 normalized() const
    {
        float len = length();
        if (len > 0)
            return *this / len;
        return *this;
    }

    float max_component() const { return std::max({x, y, z}); }

    static Vec3 min(const Vec3 &a, const Vec3 &b)
    {
        return Vec3(std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z));
    }

    static Vec3 max(const Vec3 &a, const Vec3 &b)
    {
        return Vec3(std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z));
    }
};

// RGB цвет (используется как для яркости, так и для коэффициентов)
using Color = Vec3;

// ==================== Генератор случайных чисел ====================

class RandomGenerator
{
private:
    std::mt19937 gen;
    std::uniform_real_distribution<float> dist;

public:
    RandomGenerator() : gen(std::random_device{}()), dist(0.0f, 1.0f) {}
    RandomGenerator(unsigned seed) : gen(seed), dist(0.0f, 1.0f) {}

    float uniform() { return dist(gen); }

    // Равномерная точка на полусфере
    Vec3 uniformHemisphere(const Vec3 &normal)
    {
        float z = uniform();
        float r = std::sqrt(std::max(0.0f, 1.0f - z * z));
        float phi = 2.0f * M_PI * uniform();

        Vec3 localDir(r * std::cos(phi), r * std::sin(phi), z);
        return localToWorld(localDir, normal);
    }

    // Косинусно-взвешенное направление (для выборки по значимости Ламберта)
    Vec3 cosineHemisphere(const Vec3 &normal)
    {
        float r1 = uniform();
        float r2 = uniform();

        float sinTheta = std::sqrt(r1);
        float cosTheta = std::sqrt(1.0f - r1);
        float phi = 2.0f * M_PI * r2;

        Vec3 localDir(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
        return localToWorld(localDir, normal);
    }

    // Равномерная точка внутри треугольника
    Vec3 uniformTriangle(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2)
    {
        float r1 = uniform();
        float r2 = uniform();

        if (r1 + r2 > 1.0f)
        {
            r1 = 1.0f - r1;
            r2 = 1.0f - r2;
        }

        return v0 * (1.0f - r1 - r2) + v1 * r1 + v2 * r2;
    }

private:
    // Преобразование из локальной системы координат (z = нормаль) в мировую
    Vec3 localToWorld(const Vec3 &localDir, const Vec3 &normal)
    {
        Vec3 n = normal.normalized();
        Vec3 tangent = (std::abs(n.x) > 0.9f) ? Vec3(0, 1, 0) : Vec3(1, 0, 0);
        Vec3 bitangent = n.cross(tangent).normalized();
        tangent = bitangent.cross(n);

        return tangent * localDir.x + bitangent * localDir.y + n * localDir.z;
    }
};

// ==================== Материалы ====================

enum class MaterialType
{
    DIFFUSE,  // Диффузное отражение (Ламберт)
    SPECULAR, // Зеркальное отражение
    MIXED,    // Смешанный (диффузный + зеркальный)
    EMISSIVE  // Излучающий (источник света)
};

struct Material
{
    MaterialType type;
    Color diffuseColor;  // Коэффициент диффузного отражения
    Color specularColor; // Коэффициент зеркального отражения
    Color emission;      // Излучение (для источников света)
    float specularProb;  // Вероятность зеркального отражения для смешанного материала

    Material() : type(MaterialType::DIFFUSE),
                 diffuseColor(0.8f, 0.8f, 0.8f),
                 specularColor(1.0f, 1.0f, 1.0f),
                 emission(0, 0, 0),
                 specularProb(0.0f) {}

    // Создание диффузного материала
    static Material Diffuse(const Color &color)
    {
        Material m;
        m.type = MaterialType::DIFFUSE;
        m.diffuseColor = color;
        return m;
    }

    // Создание зеркального материала
    static Material Specular(const Color &color = Color(1, 1, 1))
    {
        Material m;
        m.type = MaterialType::SPECULAR;
        m.specularColor = color;
        return m;
    }

    // Создание смешанного материала
    static Material Mixed(const Color &diffuse, const Color &specular, float specProb)
    {
        Material m;
        m.type = MaterialType::MIXED;
        m.diffuseColor = diffuse;
        m.specularColor = specular;
        m.specularProb = std::min(std::max(specProb, 0.0f), 1.0f);
        return m;
    }

    // Создание излучающего материала
    static Material Emissive(const Color &emission, const Color &diffuse = Color(0, 0, 0))
    {
        Material m;
        m.type = MaterialType::EMISSIVE;
        m.emission = emission;
        m.diffuseColor = diffuse;
        return m;
    }

    // Проверка физичности: сумма отражений <= 1 для каждого канала
    bool isPhysical() const
    {
        Color total = diffuseColor + specularColor;
        return total.x <= 1.0f && total.y <= 1.0f && total.z <= 1.0f;
    }
};

// ==================== Треугольник ====================

struct Triangle
{
    Vec3 v0, v1, v2; // Вершины
    Vec3 normal;     // Нормаль
    int materialId;  // Индекс материала
    float area;      // Площадь

    Triangle() : materialId(0), area(0) {}

    Triangle(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2, int matId = 0)
        : v0(v0), v1(v1), v2(v2), materialId(matId)
    {
        Vec3 edge1 = v1 - v0;
        Vec3 edge2 = v2 - v0;
        Vec3 cross = edge1.cross(edge2);
        normal = cross.normalized();
        area = cross.length() * 0.5f;
    }

    // Получить случайную точку на треугольнике
    Vec3 samplePoint(RandomGenerator &rng) const
    {
        return rng.uniformTriangle(v0, v1, v2);
    }
};

// ==================== Источник света ====================

struct AreaLight
{
    Triangle triangle;
    Color emission;
    float power; // Общая мощность источника

    AreaLight(const Triangle &tri, const Color &emit)
        : triangle(tri), emission(emit)
    {
        // Мощность = интеграл по площади (emission * cos(theta) * dA)
        // Для равномерного Ламберта: power = emission * area * pi
        power = (emit.x + emit.y + emit.z) / 3.0f * tri.area * M_PI;
    }
};

// ==================== Камера ====================

struct Camera
{
    Vec3 position;
    Vec3 forward;
    Vec3 right;
    Vec3 up;
    float fov; // Угол поля зрения (в радианах)
    int width, height;

    Camera(const Vec3 &pos, const Vec3 &lookAt, const Vec3 &worldUp,
           float fovDegrees, int w, int h)
        : position(pos), width(w), height(h)
    {
        fov = fovDegrees * M_PI / 180.0f;
        forward = (lookAt - pos).normalized();
        right = forward.cross(worldUp).normalized();
        up = right.cross(forward);
    }

    // Генерация луча для пикселя (x, y) с антиалиасингом
    void generateRay(int px, int py, float &ox, float &oy, float &oz,
                     float &dx, float &dy, float &dz, RandomGenerator &rng) const
    {
        // Случайное смещение внутри пикселя для антиалиасинга
        float jitterX = rng.uniform() - 0.5f;
        float jitterY = rng.uniform() - 0.5f;

        float x = (2.0f * (px + 0.5f + jitterX) / width - 1.0f);
        float y = (1.0f - 2.0f * (py + 0.5f + jitterY) / height);

        float aspectRatio = float(width) / float(height);
        float tanHalfFov = std::tan(fov * 0.5f);

        x *= aspectRatio * tanHalfFov;
        y *= tanHalfFov;

        Vec3 dir = (forward + right * x + up * y).normalized();

        ox = position.x;
        oy = position.y;
        oz = position.z;
        dx = dir.x;
        dy = dir.y;
        dz = dir.z;
    }
};

// ==================== Сцена с Embree ====================

class Scene
{
private:
    RTCDevice device;
    RTCScene scene;
    std::vector<Triangle> triangles;
    std::vector<Material> materials;
    std::vector<AreaLight> lights;
    float totalLightPower;

public:
    Scene() : totalLightPower(0)
    {
        device = rtcNewDevice(nullptr);
        if (!device)
        {
            std::cerr << "Ошибка: не удалось создать Embree device" << std::endl;
            exit(1);
        }
        scene = rtcNewScene(device);
    }

    ~Scene()
    {
        rtcReleaseScene(scene);
        rtcReleaseDevice(device);
    }

    // Добавить материал
    int addMaterial(const Material &mat)
    {
        materials.push_back(mat);
        return materials.size() - 1;
    }

    // Добавить треугольник
    void addTriangle(const Triangle &tri)
    {
        triangles.push_back(tri);
    }

    // Добавить четырехугольник (два треугольника)
    void addQuad(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2, const Vec3 &v3, int matId)
    {
        addTriangle(Triangle(v0, v1, v2, matId));
        addTriangle(Triangle(v0, v2, v3, matId));
    }

    // Добавить источник света
    void addLight(const Triangle &tri, const Color &emission)
    {
        AreaLight light(tri, emission);
        lights.push_back(light);
        totalLightPower += light.power;
    }

    // Построить BVH после добавления всей геометрии
    void commit()
    {
        RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);

        // Выделяем буферы вершин и индексов
        float *vertices = (float *)rtcSetNewGeometryBuffer(
            geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3,
            3 * sizeof(float), triangles.size() * 3);

        unsigned *indices = (unsigned *)rtcSetNewGeometryBuffer(
            geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,
            3 * sizeof(unsigned), triangles.size());

        // Заполняем буферы
        for (size_t i = 0; i < triangles.size(); ++i)
        {
            vertices[i * 9 + 0] = triangles[i].v0.x;
            vertices[i * 9 + 1] = triangles[i].v0.y;
            vertices[i * 9 + 2] = triangles[i].v0.z;
            vertices[i * 9 + 3] = triangles[i].v1.x;
            vertices[i * 9 + 4] = triangles[i].v1.y;
            vertices[i * 9 + 5] = triangles[i].v1.z;
            vertices[i * 9 + 6] = triangles[i].v2.x;
            vertices[i * 9 + 7] = triangles[i].v2.y;
            vertices[i * 9 + 8] = triangles[i].v2.z;

            indices[i * 3 + 0] = i * 3 + 0;
            indices[i * 3 + 1] = i * 3 + 1;
            indices[i * 3 + 2] = i * 3 + 2;
        }

        rtcCommitGeometry(geom);
        rtcAttachGeometry(scene, geom);
        rtcReleaseGeometry(geom);
        rtcCommitScene(scene);
    }

    // Трассировка луча
    bool intersect(float ox, float oy, float oz,
                   float dx, float dy, float dz,
                   float &t, int &triId) const
    {
        RTCRayHit rayhit;
        rayhit.ray.org_x = ox;
        rayhit.ray.org_y = oy;
        rayhit.ray.org_z = oz;
        rayhit.ray.dir_x = dx;
        rayhit.ray.dir_y = dy;
        rayhit.ray.dir_z = dz;
        rayhit.ray.tnear = 0.001f;
        rayhit.ray.tfar = std::numeric_limits<float>::infinity();
        rayhit.ray.mask = -1;
        rayhit.ray.flags = 0;
        rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
        rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;

        rtcIntersect1(scene, &rayhit);

        if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID)
        {
            t = rayhit.ray.tfar;
            triId = rayhit.hit.primID;
            return true;
        }
        return false;
    }

    // Проверка видимости (теневой луч)
    bool isVisible(const Vec3 &from, const Vec3 &to) const
    {
        Vec3 dir = to - from;
        float dist = dir.length();
        dir = dir / dist;

        RTCRay ray;
        ray.org_x = from.x;
        ray.org_y = from.y;
        ray.org_z = from.z;
        ray.dir_x = dir.x;
        ray.dir_y = dir.y;
        ray.dir_z = dir.z;
        ray.tnear = 0.001f;
        ray.tfar = dist - 0.001f;
        ray.mask = -1;
        ray.flags = 0;

        rtcOccluded1(scene, &ray);

        return ray.tfar >= 0; // Если tfar < 0, луч был заблокирован
    }

    const Triangle &getTriangle(int id) const { return triangles[id]; }
    const Material &getMaterial(int id) const { return materials[id]; }
    const std::vector<AreaLight> &getLights() const { return lights; }
    float getTotalLightPower() const { return totalLightPower; }
    size_t getTriangleCount() const { return triangles.size(); }

    // Загрузка OBJ файла
    bool loadOBJ(const std::string &filename, int materialId,
                 const Vec3 &offset = Vec3(0, 0, 0), float scale = 1.0f)
    {
        std::ifstream file(filename);
        if (!file.is_open())
        {
            std::cerr << "Ошибка: не удалось открыть OBJ файл: " << filename << std::endl;
            return false;
        }

        std::vector<Vec3> vertices;
        std::vector<Vec3> normals;
        int triangleCount = 0;

        std::string line;
        while (std::getline(file, line))
        {
            std::istringstream iss(line);
            std::string prefix;
            iss >> prefix;

            if (prefix == "v")
            {
                // Вершина
                float x, y, z;
                iss >> x >> y >> z;
                vertices.push_back(Vec3(x * scale + offset.x,
                                        y * scale + offset.y,
                                        z * scale + offset.z));
            }
            else if (prefix == "vn")
            {
                // Нормаль
                float nx, ny, nz;
                iss >> nx >> ny >> nz;
                normals.push_back(Vec3(nx, ny, nz));
            }
            else if (prefix == "f")
            {
                // Грань (face) - поддержка форматов: v, v/vt, v/vt/vn, v//vn
                std::vector<int> faceVertices;
                std::string vertexData;

                while (iss >> vertexData)
                {
                    int vIdx = 0;
                    size_t pos = vertexData.find('/');
                    if (pos != std::string::npos)
                    {
                        vIdx = std::stoi(vertexData.substr(0, pos));
                    }
                    else
                    {
                        vIdx = std::stoi(vertexData);
                    }
                    // OBJ индексы начинаются с 1
                    faceVertices.push_back(vIdx - 1);
                }

                // Триангуляция полигона (fan triangulation)
                for (size_t i = 1; i + 1 < faceVertices.size(); ++i)
                {
                    int i0 = faceVertices[0];
                    int i1 = faceVertices[i];
                    int i2 = faceVertices[i + 1];

                    if (i0 >= 0 && i0 < (int)vertices.size() &&
                        i1 >= 0 && i1 < (int)vertices.size() &&
                        i2 >= 0 && i2 < (int)vertices.size())
                    {
                        addTriangle(Triangle(vertices[i0], vertices[i1], vertices[i2], materialId));
                        triangleCount++;
                    }
                }
            }
        }

        file.close();
        std::cout << "Загружен OBJ: " << filename << " (" << vertices.size()
                  << " вершин, " << triangleCount << " треугольников)" << std::endl;
        return true;
    }
};

// ==================== Path Tracer ====================

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

// ==================== Постобработка и сохранение ====================

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

// ==================== Создание сцены Cornell Box ====================

void createCornellBox(Scene &scene)
{
    // Материалы
    int whiteMat = scene.addMaterial(Material::Diffuse(Color(0.73f, 0.73f, 0.73f)));
    int redMat = scene.addMaterial(Material::Diffuse(Color(0.65f, 0.05f, 0.05f)));
    int greenMat = scene.addMaterial(Material::Diffuse(Color(0.12f, 0.45f, 0.15f)));
    int mirrorMat = scene.addMaterial(Material::Specular(Color(0.99f, 0.99f, 0.99f)));
    int mixedMat = scene.addMaterial(Material::Mixed(Color(0.7f, 0.7f, 0.7f), Color(0.8f, 0.8f, 0.8f), 0.3f));
    int lightMat = scene.addMaterial(Material::Emissive(Color(15.0f, 15.0f, 15.0f)));

    // Размеры Cornell Box
    float size = 5.0f;

    // Пол (белый)
    scene.addQuad(
        Vec3(-size, 0, -size), Vec3(size, 0, -size),
        Vec3(size, 0, size), Vec3(-size, 0, size),
        whiteMat);

    // Потолок (белый)
    scene.addQuad(
        Vec3(-size, size * 2, size), Vec3(size, size * 2, size),
        Vec3(size, size * 2, -size), Vec3(-size, size * 2, -size),
        whiteMat);

    // Задняя стена (белая)
    scene.addQuad(
        Vec3(-size, 0, -size), Vec3(-size, size * 2, -size),
        Vec3(size, size * 2, -size), Vec3(size, 0, -size),
        whiteMat);

    // Левая стена (красная)
    scene.addQuad(
        Vec3(-size, 0, size), Vec3(-size, size * 2, size),
        Vec3(-size, size * 2, -size), Vec3(-size, 0, -size),
        redMat);

    // Правая стена (зеленая)
    scene.addQuad(
        Vec3(size, 0, -size), Vec3(size, size * 2, -size),
        Vec3(size, size * 2, size), Vec3(size, 0, size),
        greenMat);

    // Высокий блок (зеркальный)
    float bx1 = 1.5f, bz1 = -2.0f;
    float blockHeight1 = 6.0f;
    float blockSize1 = 2.5f;

    // Верх блока
    scene.addQuad(
        Vec3(bx1 - blockSize1 / 2, blockHeight1, bz1 - blockSize1 / 2),
        Vec3(bx1 + blockSize1 / 2, blockHeight1, bz1 - blockSize1 / 2),
        Vec3(bx1 + blockSize1 / 2, blockHeight1, bz1 + blockSize1 / 2),
        Vec3(bx1 - blockSize1 / 2, blockHeight1, bz1 + blockSize1 / 2),
        mirrorMat);

    // Стороны высокого блока
    scene.addQuad(
        Vec3(bx1 - blockSize1 / 2, 0, bz1 + blockSize1 / 2),
        Vec3(bx1 - blockSize1 / 2, blockHeight1, bz1 + blockSize1 / 2),
        Vec3(bx1 - blockSize1 / 2, blockHeight1, bz1 - blockSize1 / 2),
        Vec3(bx1 - blockSize1 / 2, 0, bz1 - blockSize1 / 2),
        mirrorMat);
    scene.addQuad(
        Vec3(bx1 + blockSize1 / 2, 0, bz1 - blockSize1 / 2),
        Vec3(bx1 + blockSize1 / 2, blockHeight1, bz1 - blockSize1 / 2),
        Vec3(bx1 + blockSize1 / 2, blockHeight1, bz1 + blockSize1 / 2),
        Vec3(bx1 + blockSize1 / 2, 0, bz1 + blockSize1 / 2),
        mirrorMat);
    scene.addQuad(
        Vec3(bx1 - blockSize1 / 2, 0, bz1 - blockSize1 / 2),
        Vec3(bx1 - blockSize1 / 2, blockHeight1, bz1 - blockSize1 / 2),
        Vec3(bx1 + blockSize1 / 2, blockHeight1, bz1 - blockSize1 / 2),
        Vec3(bx1 + blockSize1 / 2, 0, bz1 - blockSize1 / 2),
        mirrorMat);
    scene.addQuad(
        Vec3(bx1 + blockSize1 / 2, 0, bz1 + blockSize1 / 2),
        Vec3(bx1 + blockSize1 / 2, blockHeight1, bz1 + blockSize1 / 2),
        Vec3(bx1 - blockSize1 / 2, blockHeight1, bz1 + blockSize1 / 2),
        Vec3(bx1 - blockSize1 / 2, 0, bz1 + blockSize1 / 2),
        mirrorMat);

    // Низкий блок (смешанный материал)
    float bx2 = -1.5f, bz2 = 1.0f;
    float blockHeight2 = 3.0f;
    float blockSize2 = 2.5f;

    // Верх блока
    scene.addQuad(
        Vec3(bx2 - blockSize2 / 2, blockHeight2, bz2 - blockSize2 / 2),
        Vec3(bx2 + blockSize2 / 2, blockHeight2, bz2 - blockSize2 / 2),
        Vec3(bx2 + blockSize2 / 2, blockHeight2, bz2 + blockSize2 / 2),
        Vec3(bx2 - blockSize2 / 2, blockHeight2, bz2 + blockSize2 / 2),
        mixedMat);

    // Стороны низкого блока
    scene.addQuad(
        Vec3(bx2 - blockSize2 / 2, 0, bz2 + blockSize2 / 2),
        Vec3(bx2 - blockSize2 / 2, blockHeight2, bz2 + blockSize2 / 2),
        Vec3(bx2 - blockSize2 / 2, blockHeight2, bz2 - blockSize2 / 2),
        Vec3(bx2 - blockSize2 / 2, 0, bz2 - blockSize2 / 2),
        mixedMat);
    scene.addQuad(
        Vec3(bx2 + blockSize2 / 2, 0, bz2 - blockSize2 / 2),
        Vec3(bx2 + blockSize2 / 2, blockHeight2, bz2 - blockSize2 / 2),
        Vec3(bx2 + blockSize2 / 2, blockHeight2, bz2 + blockSize2 / 2),
        Vec3(bx2 + blockSize2 / 2, 0, bz2 + blockSize2 / 2),
        mixedMat);
    scene.addQuad(
        Vec3(bx2 - blockSize2 / 2, 0, bz2 - blockSize2 / 2),
        Vec3(bx2 - blockSize2 / 2, blockHeight2, bz2 - blockSize2 / 2),
        Vec3(bx2 + blockSize2 / 2, blockHeight2, bz2 - blockSize2 / 2),
        Vec3(bx2 + blockSize2 / 2, 0, bz2 - blockSize2 / 2),
        mixedMat);
    scene.addQuad(
        Vec3(bx2 + blockSize2 / 2, 0, bz2 + blockSize2 / 2),
        Vec3(bx2 + blockSize2 / 2, blockHeight2, bz2 + blockSize2 / 2),
        Vec3(bx2 - blockSize2 / 2, blockHeight2, bz2 + blockSize2 / 2),
        Vec3(bx2 - blockSize2 / 2, 0, bz2 + blockSize2 / 2),
        mixedMat);

    // Источник света на потолке
    float lightSize = 2.0f;
    float lightY = size * 2 - 0.01f;

    Triangle lightTri1(
        Vec3(-lightSize, lightY, -lightSize),
        Vec3(lightSize, lightY, -lightSize),
        Vec3(lightSize, lightY, lightSize),
        lightMat);
    Triangle lightTri2(
        Vec3(-lightSize, lightY, -lightSize),
        Vec3(lightSize, lightY, lightSize),
        Vec3(-lightSize, lightY, lightSize),
        lightMat);

    scene.addTriangle(lightTri1);
    scene.addTriangle(lightTri2);
    scene.addLight(lightTri1, Color(15.0f, 15.0f, 15.0f));
    scene.addLight(lightTri2, Color(15.0f, 15.0f, 15.0f));
}

// ==================== Сцена: Комната с зеркалами ====================

void createMirrorRoom(Scene &scene)
{
    // Материалы
    int whiteMat = scene.addMaterial(Material::Diffuse(Color(0.9f, 0.9f, 0.9f)));
    int blueMat = scene.addMaterial(Material::Diffuse(Color(0.1f, 0.3f, 0.7f)));
    int mirrorMat = scene.addMaterial(Material::Specular(Color(0.95f, 0.95f, 0.95f)));
    int goldMat = scene.addMaterial(Material::Mixed(Color(0.3f, 0.2f, 0.1f), Color(0.9f, 0.7f, 0.3f), 0.6f));
    int lightMat = scene.addMaterial(Material::Emissive(Color(20.0f, 20.0f, 18.0f)));

    float size = 5.0f;

    // Пол (шахматный паттерн - упрощенно белый)
    scene.addQuad(
        Vec3(-size, 0, -size), Vec3(size, 0, -size),
        Vec3(size, 0, size), Vec3(-size, 0, size),
        whiteMat);

    // Потолок
    scene.addQuad(
        Vec3(-size, size * 2, size), Vec3(size, size * 2, size),
        Vec3(size, size * 2, -size), Vec3(-size, size * 2, -size),
        whiteMat);

    // Задняя стена (зеркало)
    scene.addQuad(
        Vec3(-size, 0, -size), Vec3(-size, size * 2, -size),
        Vec3(size, size * 2, -size), Vec3(size, 0, -size),
        mirrorMat);

    // Левая стена (синяя)
    scene.addQuad(
        Vec3(-size, 0, size), Vec3(-size, size * 2, size),
        Vec3(-size, size * 2, -size), Vec3(-size, 0, -size),
        blueMat);

    // Правая стена (зеркало)
    scene.addQuad(
        Vec3(size, 0, -size), Vec3(size, size * 2, -size),
        Vec3(size, size * 2, size), Vec3(size, 0, size),
        mirrorMat);

    // Сфера (аппроксимация октаэдром для простоты)
    float cx = 0, cy = 2.5f, cz = -1.0f;
    float r = 1.5f;

    // Верхняя половина
    scene.addTriangle(Triangle(Vec3(cx, cy + r, cz), Vec3(cx + r, cy, cz), Vec3(cx, cy, cz + r), goldMat));
    scene.addTriangle(Triangle(Vec3(cx, cy + r, cz), Vec3(cx, cy, cz + r), Vec3(cx - r, cy, cz), goldMat));
    scene.addTriangle(Triangle(Vec3(cx, cy + r, cz), Vec3(cx - r, cy, cz), Vec3(cx, cy, cz - r), goldMat));
    scene.addTriangle(Triangle(Vec3(cx, cy + r, cz), Vec3(cx, cy, cz - r), Vec3(cx + r, cy, cz), goldMat));

    // Нижняя половина
    scene.addTriangle(Triangle(Vec3(cx, cy - r, cz), Vec3(cx, cy, cz + r), Vec3(cx + r, cy, cz), goldMat));
    scene.addTriangle(Triangle(Vec3(cx, cy - r, cz), Vec3(cx - r, cy, cz), Vec3(cx, cy, cz + r), goldMat));
    scene.addTriangle(Triangle(Vec3(cx, cy - r, cz), Vec3(cx, cy, cz - r), Vec3(cx - r, cy, cz), goldMat));
    scene.addTriangle(Triangle(Vec3(cx, cy - r, cz), Vec3(cx + r, cy, cz), Vec3(cx, cy, cz - r), goldMat));

    // Источник света
    float lightSize = 1.5f;
    float lightY = size * 2 - 0.01f;

    Triangle lightTri1(Vec3(-lightSize, lightY, -lightSize), Vec3(lightSize, lightY, -lightSize), Vec3(lightSize, lightY, lightSize), lightMat);
    Triangle lightTri2(Vec3(-lightSize, lightY, -lightSize), Vec3(lightSize, lightY, lightSize), Vec3(-lightSize, lightY, lightSize), lightMat);

    scene.addTriangle(lightTri1);
    scene.addTriangle(lightTri2);
    scene.addLight(lightTri1, Color(20.0f, 20.0f, 18.0f));
    scene.addLight(lightTri2, Color(20.0f, 20.0f, 18.0f));
}

// ==================== Сцена: Три сферы ====================

void createThreeSpheres(Scene &scene)
{
    // Материалы
    int groundMat = scene.addMaterial(Material::Diffuse(Color(0.5f, 0.5f, 0.5f)));
    int diffuseMat = scene.addMaterial(Material::Diffuse(Color(0.8f, 0.2f, 0.2f)));
    int mirrorMat = scene.addMaterial(Material::Specular(Color(0.9f, 0.9f, 0.9f)));
    int glassMat = scene.addMaterial(Material::Mixed(Color(0.1f, 0.1f, 0.1f), Color(0.9f, 0.9f, 0.9f), 0.9f));
    int skyMat = scene.addMaterial(Material::Emissive(Color(0.5f, 0.7f, 1.0f)));

    // Большая плоскость (пол)
    float groundSize = 20.0f;
    scene.addQuad(
        Vec3(-groundSize, 0, -groundSize), Vec3(groundSize, 0, -groundSize),
        Vec3(groundSize, 0, groundSize), Vec3(-groundSize, 0, groundSize),
        groundMat);

    // Функция для создания икосаэдра (аппроксимация сферы)
    auto addSphere = [&scene](float cx, float cy, float cz, float r, int matId)
    {
        // Упрощенная сфера из 8 треугольников (октаэдр)
        scene.addTriangle(Triangle(Vec3(cx, cy + r, cz), Vec3(cx + r, cy, cz), Vec3(cx, cy, cz + r), matId));
        scene.addTriangle(Triangle(Vec3(cx, cy + r, cz), Vec3(cx, cy, cz + r), Vec3(cx - r, cy, cz), matId));
        scene.addTriangle(Triangle(Vec3(cx, cy + r, cz), Vec3(cx - r, cy, cz), Vec3(cx, cy, cz - r), matId));
        scene.addTriangle(Triangle(Vec3(cx, cy + r, cz), Vec3(cx, cy, cz - r), Vec3(cx + r, cy, cz), matId));
        scene.addTriangle(Triangle(Vec3(cx, cy - r, cz), Vec3(cx, cy, cz + r), Vec3(cx + r, cy, cz), matId));
        scene.addTriangle(Triangle(Vec3(cx, cy - r, cz), Vec3(cx - r, cy, cz), Vec3(cx, cy, cz + r), matId));
        scene.addTriangle(Triangle(Vec3(cx, cy - r, cz), Vec3(cx, cy, cz - r), Vec3(cx - r, cy, cz), matId));
        scene.addTriangle(Triangle(Vec3(cx, cy - r, cz), Vec3(cx + r, cy, cz), Vec3(cx, cy, cz - r), matId));
    };

    // Три сферы
    addSphere(-3.0f, 1.5f, 0, 1.5f, diffuseMat); // Красная диффузная
    addSphere(0, 1.5f, -2.0f, 1.5f, mirrorMat);  // Зеркальная
    addSphere(3.0f, 1.5f, 0, 1.5f, glassMat);    // Полузеркальная

    // Небесная сфера (источник света)
    float skyR = 100.0f;
    // Упрощенно - просто большой квад сверху
    Triangle skyTri1(Vec3(-50, 30, -50), Vec3(50, 30, -50), Vec3(50, 30, 50), skyMat);
    Triangle skyTri2(Vec3(-50, 30, -50), Vec3(50, 30, 50), Vec3(-50, 30, 50), skyMat);

    scene.addTriangle(skyTri1);
    scene.addTriangle(skyTri2);
    scene.addLight(skyTri1, Color(0.5f, 0.7f, 1.0f));
    scene.addLight(skyTri2, Color(0.5f, 0.7f, 1.0f));
}

// ==================== Сцена: Цветное освещение ====================

void createColoredLights(Scene &scene)
{
    // Материалы
    int whiteMat = scene.addMaterial(Material::Diffuse(Color(0.9f, 0.9f, 0.9f)));
    int mirrorMat = scene.addMaterial(Material::Specular(Color(0.95f, 0.95f, 0.95f)));
    int redLightMat = scene.addMaterial(Material::Emissive(Color(25.0f, 2.0f, 2.0f)));
    int greenLightMat = scene.addMaterial(Material::Emissive(Color(2.0f, 25.0f, 2.0f)));
    int blueLightMat = scene.addMaterial(Material::Emissive(Color(2.0f, 2.0f, 25.0f)));

    float size = 5.0f;

    // Белые стены
    // Пол
    scene.addQuad(Vec3(-size, 0, -size), Vec3(size, 0, -size), Vec3(size, 0, size), Vec3(-size, 0, size), whiteMat);
    // Потолок
    scene.addQuad(Vec3(-size, size * 2, size), Vec3(size, size * 2, size), Vec3(size, size * 2, -size), Vec3(-size, size * 2, -size), whiteMat);
    // Задняя стена
    scene.addQuad(Vec3(-size, 0, -size), Vec3(-size, size * 2, -size), Vec3(size, size * 2, -size), Vec3(size, 0, -size), whiteMat);
    // Левая стена
    scene.addQuad(Vec3(-size, 0, size), Vec3(-size, size * 2, size), Vec3(-size, size * 2, -size), Vec3(-size, 0, -size), whiteMat);
    // Правая стена
    scene.addQuad(Vec3(size, 0, -size), Vec3(size, size * 2, -size), Vec3(size, size * 2, size), Vec3(size, 0, size), whiteMat);

    // Зеркальный шар в центре
    float cx = 0, cy = 3.0f, cz = -1.0f, r = 2.0f;
    scene.addTriangle(Triangle(Vec3(cx, cy + r, cz), Vec3(cx + r, cy, cz), Vec3(cx, cy, cz + r), mirrorMat));
    scene.addTriangle(Triangle(Vec3(cx, cy + r, cz), Vec3(cx, cy, cz + r), Vec3(cx - r, cy, cz), mirrorMat));
    scene.addTriangle(Triangle(Vec3(cx, cy + r, cz), Vec3(cx - r, cy, cz), Vec3(cx, cy, cz - r), mirrorMat));
    scene.addTriangle(Triangle(Vec3(cx, cy + r, cz), Vec3(cx, cy, cz - r), Vec3(cx + r, cy, cz), mirrorMat));
    scene.addTriangle(Triangle(Vec3(cx, cy - r, cz), Vec3(cx, cy, cz + r), Vec3(cx + r, cy, cz), mirrorMat));
    scene.addTriangle(Triangle(Vec3(cx, cy - r, cz), Vec3(cx - r, cy, cz), Vec3(cx, cy, cz + r), mirrorMat));
    scene.addTriangle(Triangle(Vec3(cx, cy - r, cz), Vec3(cx, cy, cz - r), Vec3(cx - r, cy, cz), mirrorMat));
    scene.addTriangle(Triangle(Vec3(cx, cy - r, cz), Vec3(cx + r, cy, cz), Vec3(cx, cy, cz - r), mirrorMat));

    // Три цветных источника света
    float lightY = size * 2 - 0.01f;
    float lightSize = 1.0f;

    // Красный свет (слева)
    Triangle redLight1(Vec3(-3.5f - lightSize, lightY, -lightSize), Vec3(-3.5f + lightSize, lightY, -lightSize), Vec3(-3.5f + lightSize, lightY, lightSize), redLightMat);
    Triangle redLight2(Vec3(-3.5f - lightSize, lightY, -lightSize), Vec3(-3.5f + lightSize, lightY, lightSize), Vec3(-3.5f - lightSize, lightY, lightSize), redLightMat);
    scene.addTriangle(redLight1);
    scene.addTriangle(redLight2);
    scene.addLight(redLight1, Color(25.0f, 2.0f, 2.0f));
    scene.addLight(redLight2, Color(25.0f, 2.0f, 2.0f));

    // Зеленый свет (центр)
    Triangle greenLight1(Vec3(-lightSize, lightY, -lightSize), Vec3(lightSize, lightY, -lightSize), Vec3(lightSize, lightY, lightSize), greenLightMat);
    Triangle greenLight2(Vec3(-lightSize, lightY, -lightSize), Vec3(lightSize, lightY, lightSize), Vec3(-lightSize, lightY, lightSize), greenLightMat);
    scene.addTriangle(greenLight1);
    scene.addTriangle(greenLight2);
    scene.addLight(greenLight1, Color(2.0f, 25.0f, 2.0f));
    scene.addLight(greenLight2, Color(2.0f, 25.0f, 2.0f));

    // Синий свет (справа)
    Triangle blueLight1(Vec3(3.5f - lightSize, lightY, -lightSize), Vec3(3.5f + lightSize, lightY, -lightSize), Vec3(3.5f + lightSize, lightY, lightSize), blueLightMat);
    Triangle blueLight2(Vec3(3.5f - lightSize, lightY, -lightSize), Vec3(3.5f + lightSize, lightY, lightSize), Vec3(3.5f - lightSize, lightY, lightSize), blueLightMat);
    scene.addTriangle(blueLight1);
    scene.addTriangle(blueLight2);
    scene.addLight(blueLight1, Color(2.0f, 2.0f, 25.0f));
    scene.addLight(blueLight2, Color(2.0f, 2.0f, 25.0f));
}

// ==================== Вспомогательные функции ====================

void printUsage(const char *programName)
{
    std::cout << "Использование: " << programName << " [опции]" << std::endl;
    std::cout << "Опции:" << std::endl;
    std::cout << "  -w <ширина>      Ширина изображения (по умолчанию: 512)" << std::endl;
    std::cout << "  -h <высота>      Высота изображения (по умолчанию: 512)" << std::endl;
    std::cout << "  -s <samples>     Число лучей на пиксель (по умолчанию: 64)" << std::endl;
    std::cout << "  -d <depth>       Максимальная глубина трассировки (по умолчанию: 10)" << std::endl;
    std::cout << "  -o <файл>        Имя выходного файла (по умолчанию: output.ppm)" << std::endl;
    std::cout << "  -e <exposure>    Экспозиция (по умолчанию: авто)" << std::endl;
    std::cout << "  -g <gamma>       Гамма-коррекция (по умолчанию: 2.2)" << std::endl;
    std::cout << "  --obj <файл>     Загрузить OBJ модель" << std::endl;
    std::cout << "  --scale <число>  Масштаб OBJ модели (по умолчанию: 1.0)" << std::endl;
    std::cout << "  --help           Показать эту справку" << std::endl;
}

// ==================== Main ====================

int main(int argc, char *argv[])
{
    // Параметры по умолчанию
    int width = 512;
    int height = 512;
    int samplesPerPixel = 64;
    int maxDepth = 10;
    std::string outputFile = "output.ppm";
    std::string objFile = "";
    float objScale = 1.0f;
    float exposure = -1.0f; // -1 означает автоэкспозицию
    float gamma = 2.2f;

    // Парсинг аргументов командной строки
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];

        if (arg == "--help")
        {
            printUsage(argv[0]);
            return 0;
        }
        else if (arg == "-w" && i + 1 < argc)
        {
            width = std::atoi(argv[++i]);
        }
        else if (arg == "-h" && i + 1 < argc)
        {
            height = std::atoi(argv[++i]);
        }
        else if (arg == "-s" && i + 1 < argc)
        {
            samplesPerPixel = std::atoi(argv[++i]);
        }
        else if (arg == "-d" && i + 1 < argc)
        {
            maxDepth = std::atoi(argv[++i]);
        }
        else if (arg == "-o" && i + 1 < argc)
        {
            outputFile = argv[++i];
        }
        else if (arg == "-e" && i + 1 < argc)
        {
            exposure = std::atof(argv[++i]);
        }
        else if (arg == "-g" && i + 1 < argc)
        {
            gamma = std::atof(argv[++i]);
        }
        else if (arg == "--obj" && i + 1 < argc)
        {
            objFile = argv[++i];
        }
        else if (arg == "--scale" && i + 1 < argc)
        {
            objScale = std::atof(argv[++i]);
        }
    }

    // Проверка параметров
    width = std::max(100, std::min(width, 2000));
    height = std::max(100, std::min(height, 2000));
    samplesPerPixel = std::max(1, std::min(samplesPerPixel, 10000));
    maxDepth = std::max(1, std::min(maxDepth, 100));

    std::cout << "========================================" << std::endl;
    std::cout << "    Path Tracer - Лабораторная работа 5" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Параметры:" << std::endl;
    std::cout << "  Разрешение: " << width << "x" << height << std::endl;
    std::cout << "  Лучей на пиксель: " << samplesPerPixel << std::endl;
    std::cout << "  Макс. глубина: " << maxDepth << std::endl;
    std::cout << "  Выходной файл: " << outputFile << std::endl;
    std::cout << "  Гамма: " << gamma << std::endl;
    if (!objFile.empty())
    {
        std::cout << "  OBJ файл: " << objFile << " (масштаб: " << objScale << ")" << std::endl;
    }
    std::cout << "========================================" << std::endl;

    // Создание сцены
    Scene scene;

    if (!objFile.empty())
    {
        // Создаём простую сцену с OBJ моделью
        int whiteMat = scene.addMaterial(Material::Diffuse(Color(0.73f, 0.73f, 0.73f)));
        int objMat = scene.addMaterial(Material::Diffuse(Color(0.8f, 0.3f, 0.3f)));
        int lightMat = scene.addMaterial(Material::Emissive(Color(15.0f, 15.0f, 15.0f)));

        // Пол
        float floorSize = 10.0f;
        scene.addQuad(
            Vec3(-floorSize, 0, -floorSize), Vec3(floorSize, 0, -floorSize),
            Vec3(floorSize, 0, floorSize), Vec3(-floorSize, 0, floorSize),
            whiteMat);

        // Загружаем OBJ модель
        scene.loadOBJ(objFile, objMat, Vec3(0, 0, 0), objScale);

        // Источник света
        Triangle light1(Vec3(-3, 10, -3), Vec3(3, 10, -3), Vec3(3, 10, 3), lightMat);
        Triangle light2(Vec3(-3, 10, -3), Vec3(3, 10, 3), Vec3(-3, 10, 3), lightMat);
        scene.addTriangle(light1);
        scene.addTriangle(light2);
        scene.addLight(light1, Color(15.0f, 15.0f, 15.0f));
        scene.addLight(light2, Color(15.0f, 15.0f, 15.0f));
    }
    else
    {
        createColoredLights(scene);
    }

    scene.commit();

    std::cout << "Сцена создана: " << scene.getTriangleCount() << " треугольников" << std::endl;

    // Создание камеры (разная позиция для OBJ и Cornell Box)
    Vec3 cameraPos, cameraLookAt;
    Vec3 cameraUp(0, 1, 0);

    if (!objFile.empty())
    {
        // Для OBJ модели: камера смотрит на центр сцены
        cameraPos = Vec3(5.0f * objScale, 3.0f * objScale, 5.0f * objScale);
        cameraLookAt = Vec3(0, 1.0f * objScale, 0);
    }
    else
    {
        // Для Cornell Box
        cameraPos = Vec3(0, 5.0f, 14.0f);
        cameraLookAt = Vec3(0, 5.0f, 0);
    }

    Camera camera(cameraPos, cameraLookAt, cameraUp, 45.0f, width, height);

    // Рендеринг
    PathTracer tracer(scene, camera, maxDepth, samplesPerPixel);
    std::vector<Color> hdrImage = tracer.render();

    // Постобработка и сохранение
    bool autoExposure = (exposure < 0);
    if (exposure < 0)
        exposure = 1.0f;

    std::vector<unsigned char> ldrImage = tonemapAndGammaCorrect(hdrImage, exposure, gamma, autoExposure);
    savePPM(outputFile, ldrImage, width, height);

    std::cout << "========================================" << std::endl;
    std::cout << "    Рендеринг завершен!" << std::endl;
    std::cout << "========================================" << std::endl;

    return 0;
}
